import pandas as pd
import networkx as nx
from cobra.io import read_sbml_model
from pyvis.network import Network
from collections import defaultdict
import argparse
import os
from datetime import datetime

def get_locus_tags_for_reaction(model, reaction_id):
    """
    Obtain the locus tags associated with the genes in a reaction.

    :param model: COBRA model.
    :param reaction_id: Reaction IDs.

    :returns List of locus tags associated with the reaction.
    """

    reaction = model.reactions.get_by_id(reaction_id)

    # If no genes are associated, set "Unknown"
    if reaction.genes:
        genes = [gene.id for gene in reaction.genes]
    else:
        genes = ["Unknown"]

    return genes

def total_metabolite_frequency(model, output_file):
    """
    Calculates the number of reactions in which each metabolite participates and saves it to a file.

    :param model: COBRA model.
    :param output_file: File name to save the frequency table.

    :returns DataFrame with metabolites and their frequencies.
    """
    frequencies = defaultdict(int)
    for reaction in model.reactions:
        for metabolite in reaction.metabolites:
            frequencies[metabolite.id] += 1

    if not frequencies:
        print("Warning: No metabolites found in model. Creating empty frequency table.")
        df_frequencies = pd.DataFrame(columns=["Metabolite", "Frequency"])
    else:
        df_frequencies = pd.DataFrame(
            list(frequencies.items()), columns=["Metabolite", "Frequency"]
        )
        df_frequencies.sort_values(by="Frequency", ascending=False, inplace=True)
    
    df_frequencies.to_csv(output_file, index=False, sep = "\t")

    return df_frequencies

def filter_selected_metabolites(df_frequencies, modified_file = None):
    """
    Reads a modified file of metabolites with their frequencies and returns the list of selected metabolites.
    The user can provide a modified file with the metabolites to be selected. The file must contain a column named "Metabolite".
    If no modified file is provided, the function will return the list of metabolites with a frequency higher than 20.

    :param modified_file: Path to the modified file.

    :returns List of selected metabolites.
    """

    if modified_file is None:
        selected = df_frequencies[df_frequencies["Frequency"] > 20]["Metabolite"].tolist()
    else:
        df_modified = pd.read_csv(modified_file, sep = "\t")
        if "Metabolite" not in df_modified.columns:
            raise ValueError(f"Modified file must contain a 'Metabolite' column. Found columns: {df_modified.columns.tolist()}")
        selected = df_modified["Metabolite"].tolist()
    
    if not selected:
        print("Warning: No metabolites selected for filtering.")

    return selected

def metabolite_reactant_product_frequency(model, output_file):
    """
    Calculates how many times each metabolite appears as a product or reactant in the metabolic model.
    It considers reversible reactions (lower_bound < 0) to account for metabolites that can be both products and reactants.

    :param model: COBRA model.
    :param output_file: File name to save the frequency table.

    :returns DataFrame with metabolites, frequency as product, and as reactant.
    """
    frequencies = {}

    for reaction in model.reactions:
        for metabolite, coeff in reaction.metabolites.items():
            if coeff < 0:  # Reactant
                frequencies.setdefault(metabolite.id, {"Reactant": 0, "Product": 0})
                frequencies[metabolite.id]["Reactant"] += 1
                # If reaction is reversible, consider metabolite as potential product
                if reaction.lower_bound < 0:
                    frequencies[metabolite.id]["Product"] += 1
            if coeff > 0:  # Product
                frequencies.setdefault(metabolite.id, {"Reactant": 0, "Product": 0})
                frequencies[metabolite.id]["Product"] += 1
                # If reaction is reversible, consider metabolite as potential reactant
                if reaction.lower_bound < 0:
                    frequencies[metabolite.id]["Reactant"] += 1

    if not frequencies:
        print("Warning: No metabolite frequencies found. Creating empty table.")
        df_metabolites = pd.DataFrame(columns=["Metabolite", "Reactant", "Product"])
    else:
        df_metabolites = pd.DataFrame.from_dict(frequencies, orient="index").reset_index()
        df_metabolites.rename(columns={"index": "Metabolite"}, inplace=True)

    df_metabolites.to_csv(output_file, index=False, sep = "\t")

    return df_metabolites

def metabolite_reaction_mapper(model):
    """
    Obtain a dictionary with metabolites as keys and lists of reactions where they are substrates and products.
    It considers reversible reactions (lower_bound < 0) to account for metabolites that can be both products and reactants.

    :param model: COBRA model.

    :returns Dictionary with metabolites as keys, and values as a dictionary containing lists of reactions where they are substrates and products.
    """
    metabolite_reactions = {}

    for reaction in model.reactions:
        for metabolite, coeff in reaction.metabolites.items():
            if metabolite.id not in metabolite_reactions:
                metabolite_reactions[metabolite.id] = {"Reactant": [], "Product": []}

            if coeff < 0:  # Reactant
                metabolite_reactions[metabolite.id]["Reactant"].append(reaction.id)
                # If reaction is reversible, consider metabolite as potential product
                if reaction.lower_bound < 0:
                    metabolite_reactions[metabolite.id]["Product"].append(reaction.id)
            if coeff > 0:  # Product
                metabolite_reactions[metabolite.id]["Product"].append(reaction.id)
                # If reaction is reversible, consider metabolite as potential reactant
                if reaction.lower_bound < 0:
                    metabolite_reactions[metabolite.id]["Reactant"].append(reaction.id)

    # Ensure unique reaction IDs in the lists
    for metabolite in metabolite_reactions:
        metabolite_reactions[metabolite]["Reactant"] = list(set(metabolite_reactions[metabolite]["Reactant"]))
        metabolite_reactions[metabolite]["Product"] = list(set(metabolite_reactions[metabolite]["Product"]))

    return metabolite_reactions


def identify_chokepoint_reactions_from_df(df_metabolites, model, output_dir):
    """
    Identify chokepoint reactions in the metabolic model based on the frequency of metabolites as reactants and products.
    If a metabolite only appears as a reactant in a single reaction, the reaction is considered a consumption chokepoint.
    If a metabolite only appears as a product in a single reaction, the reaction is considered a production chokepoint.

    :param df_metabolites (pd.DataFrame): DataFrame with metabolite frequency as reactants and products. Must contain columns ["Metabolite", "Reactant", "Product"].
    :param model (cobra.Model): COBRApy metabolic model.
    :param output_dir (str): Output directory to save the results.

    :returns pd.DataFrame: DataFrame listing chokepoint reactions categorized as production or consumption.
    """

    # Validate input DataFrame
    if df_metabolites.empty:
        print("Warning: Empty metabolite DataFrame. No chokepoints can be identified.")
        return pd.DataFrame(columns=["Reaction", "Name", "Gene", "Consumption_Chokepoint", "Production_Chokepoint"])
    
    required_cols = ["Metabolite", "Reactant", "Product"]
    missing_cols = [col for col in required_cols if col not in df_metabolites.columns]
    if missing_cols:
        raise ValueError(f"DataFrame missing required columns: {missing_cols}")

    # Filter metabolites that only participate as reactants in a single reaction
    unique_consumers = set(df_metabolites[df_metabolites["Reactant"] == 1]["Metabolite"].tolist())
    # Filter metabolites that only participate as products in a single reaction
    unique_producers = set(df_metabolites[df_metabolites["Product"] == 1]["Metabolite"].tolist())

    chokepoint_reactions = []

    # Iterate through reactions in the model

    for reaction in model.reactions:
        # Identify reactants and products in the reaction
        reactant_ids = {met.id for met in reaction.reactants}
        product_ids = {met.id for met in reaction.products}

        # Check if the reaction is a chokepoint for consumption or production
        is_consumption_chokepoint = bool(reactant_ids & set(unique_consumers))
        is_production_chokepoint = bool(product_ids & set(unique_producers))

        if is_consumption_chokepoint or is_production_chokepoint:
            chokepoint_reactions.append({
                "Reaction": reaction.id,
                "Name": reaction.name,
                "Gene": [gene.id for gene in reaction.genes] if reaction.genes else ["Unknown"],
                "Consumption_Chokepoint": is_consumption_chokepoint,
                "Production_Chokepoint": is_production_chokepoint
            })

    # Convert results to a DataFrame
    chokepoint_df = pd.DataFrame(chokepoint_reactions)
    output_file = os.path.join(output_dir, "chokepoint_reactions.tsv")
    chokepoint_df.to_csv(output_file, index=False, sep = "\t")

    print(f"Found {len(chokepoint_df)} chokepoint reactions.")
    print(f"Consumption chokepoints: {chokepoint_df['Consumption_Chokepoint'].sum()}")
    print(f"Production chokepoints: {chokepoint_df['Production_Chokepoint'].sum()}")
    
    return chokepoint_df

def map_genes_to_chokepoints(chokepoint_df, model, output_dir):
    """
    Map genes to chokepoint reactions in the metabolic model.

    :param chokepoint_df: DataFrame with chokepoint reactions.
    :param model: COBRA model.
    :param output_dir: Output directory to save the results.

    :returns DataFrame with genes associated with chokepoint reactions.
    """
    data = []

    for _, row in chokepoint_df.iterrows():
        try:
            reaction = model.reactions.get_by_id(row["Reaction"])
            for locus_tag in row["Gene"]:
                data.append({
                    "Gene": locus_tag,
                    "Reaction": row["Reaction"],
                    "EC_Number": reaction.annotation.get("ec-code", None),
                    "KEGG_id": reaction.annotation.get("kegg.reaction", None),
                    "Rhea_id": reaction.annotation.get("rhea", None),
                    "Biocyc_id": reaction.annotation.get("biocyc", None),
                    "Consumption_Chokepoint": row["Consumption_Chokepoint"],
                    "Production_Chokepoint": row["Production_Chokepoint"]
                })
        except KeyError:
            print(f"Warning: Reaction {row['Reaction']} not found in model. Skipping.")
            continue

    if not data:
        print("Warning: No gene data to save for chokepoints.")
        gene_df = pd.DataFrame(columns=["Gene", "Reaction", "EC_Number", "KEGG_id", "Rhea_id", "Biocyc_id", "Consumption_Chokepoint", "Production_Chokepoint"])
    else:
        gene_df = pd.DataFrame(data)
    
    output_file = os.path.join(output_dir, "chokepoint_genes.tsv")
    gene_df.to_csv(output_file, index=False, sep = "\t")

    return gene_df

def metabolism_to_graph(model, excluded_metabolites, metabolite_reactions):
    """
    Transforms the metabolic model into a graph where nodes are reactions and edges are established if the product 
    of one reaction is the substrate of another, excluding specific metabolites.

    :param model: COBRA model.
    :param excluded_metabolites: List of metabolites to exclude in connections.

    :returns NetworkX graph representing the metabolic network.
    """

    graph = nx.DiGraph()

    for reaction in model.reactions:
        graph.add_node(reaction.id)

    for reaction in model.reactions:
        if reaction.lower_bound < 0 and reaction.upper_bound > 0:
            products = {met.id for met, coeff in reaction.metabolites.items()}
        else:
            products = {met.id for met, coeff in reaction.metabolites.items() if coeff > 0}

        for product in products:
            if product not in excluded_metabolites:
                if product in metabolite_reactions:
                    reactant_rxn = metabolite_reactions[product]['Reactant']
                    for rxn in reactant_rxn:
                        if rxn != reaction.id:
                            graph.add_edge(reaction.id, rxn)
                else:
                    print(f"Warning: Metabolite {product} not found in metabolite_reactions mapping. Skipping.")

    return graph

def calculate_betweenness_centrality(graph, model, output_dir):
    """
    Calculate the betweenness centrality of reactions in the largest connected component of the metabolic network.
    The betweenness centrality is normalized by the maximum value in the network.

    :param graph: NetworkX graph representing the metabolic network.
    :param model: COBRA model.
    :param output_dir: Output directory to save the results.

    :returns Tuple with the largest connected component of the graph and a dictionary with normalized betweenness centrality values.
    """

    # Check if the graph is empty or has no nodes
    if graph.number_of_nodes() == 0:
        print("Warning: Graph is empty. Cannot calculate betweenness centrality.")
        return None, {}
    
    # Check if there are any weakly connected components
    components = list(nx.weakly_connected_components(graph))
    if not components:
        print("Warning: Graph has no connected components. Cannot calculate betweenness centrality.")
        return None, {}
    
    largest_component = max(components, key=len)
    subgraph = graph.subgraph(largest_component).copy()
    
    # Check if the subgraph has at least one node
    if subgraph.number_of_nodes() == 0:
        print("Warning: Largest connected component is empty. Cannot calculate betweenness centrality.")
        return None, {}

    # Calculate centrality measures
    all_betweenness_centrality = nx.betweenness_centrality(subgraph)
    
    # Handle case where all centrality values are zero
    max_value = max(all_betweenness_centrality.values())
    if max_value == 0:
        print("Warning: All betweenness centrality values are zero. Using raw values.")
        normalized_betweenness_centrality = all_betweenness_centrality
    else:
        normalized_betweenness_centrality = {key: (value / max_value) for key, value in all_betweenness_centrality.items()}

    # Calculate node degrees
    node_degrees = dict(subgraph.degree())

    # Create an empty list to store the data for DataFrame
    data = []

    # Iterate over reactions and their associated genes
    for reaction_id in normalized_betweenness_centrality:
        try:
            reaction = model.reactions.get_by_id(reaction_id)
            genes = get_locus_tags_for_reaction(model, reaction_id)
            for gene in genes:
                # Append a row for each gene associated with the reaction
                data.append({
                    'Gene': gene,
                    'Reaction': reaction_id,
                    'EC_Number': reaction.annotation.get('ec-code', None),
                    'KEGG_id': reaction.annotation.get('kegg.reaction', None),
                    'Rhea_id': reaction.annotation.get('rhea', None),
                    'Biocyc_id': reaction.annotation.get('biocyc', None),
                    'Betweenness_Centrality': normalized_betweenness_centrality[reaction_id],
                    'Degree': node_degrees[reaction_id]
                })
        except KeyError:
            print(f"Warning: Reaction {reaction_id} not found in model. Skipping.")
            continue

    # Create the DataFrame from the list of dictionaries
    if not data:
        print("Warning: No data to save for betweenness centrality.")
        return subgraph, normalized_betweenness_centrality
    
    df = pd.DataFrame(data)

    # Save the DataFrame to a file
    output_file = os.path.join(output_dir, "betweenness_centrality.tsv")
    df.to_csv(output_file, index=False, sep = "\t")

    return subgraph, normalized_betweenness_centrality


def write_sif(graph, file_path):
    """
    Write a networkx graph to a SIF file.

    :param graph: Networkx graph.
    :param file_path: Path to the SIF file.
    """
    try:
        with open(file_path, 'w') as file:
            for edge in graph.edges():
                file.write(f"{edge[0]}\tinteracts\t{edge[1]}\n")
    except IOError as e:
        print(f"Error writing SIF file {file_path}: {e}")

def write_html(graph, file_path):
    """
    Write a networkx graph to a HTML file.

    :param graph: Networkx graph.
    :param file_path: Path to the HTML file.
    """
    try:
        net = Network()
        net.from_nx(graph)
        net.save_graph(file_path)
    except Exception as e:
        print(f"Error writing HTML file {file_path}: {e}")

def print_box_dynamic(messages):
    width = max(len(line) for line in messages) + 4

    print('+' + '-' * (width - 2) + '+')
    for line in messages:
        print(f'| {line.ljust(width - 4)} |')
    print('+' + '-' * (width - 2) + '+')

def main():

    parser = argparse.ArgumentParser(
            description="Analyze metabolic models and generate visualizations."
        )
    parser.add_argument(
        "--model", type=str, required=True, 
        help="Path to the SBML file representing the metabolic model."
    )
    parser.add_argument(
        "--chokepoints", action="store_true", 
        help="Identify chokepoint reactions and save results to a file."
    )
    parser.add_argument(
        "--graph", action="store_true", 
        help="Generate a metabolic network graph."
    )
    parser.add_argument(
        "--frequency_filter_file", type=str, 
        help="Optional modified metabolite file to filter frequencies. Only if --graph is selected."
    )
    parser.add_argument(
        "--output", type=str, 
        help="Output directory to save the folder with the results."
    )

    args = parser.parse_args()
    
    # Validate arguments
    if args.frequency_filter_file and not args.graph:
        parser.error("--frequency_filter_file can only be used with --graph")
    
    if not args.chokepoints and not args.graph:
        parser.error("At least one of --chokepoints or --graph must be specified")

    if args.output:
        output_dir = args.output
    else:
        output_dir = os.getcwd()
    
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            print(f"Error: Cannot create output directory {output_dir}: {e}")
            return
    
    if not os.access(output_dir, os.W_OK):
        print(f"Error: Output directory {output_dir} is not writable.")
        return

    today = datetime.today().strftime('%Y-%m-%d_%H-%M-%S')
    results_dir = os.path.join(output_dir, f"MGT_results_{today}")
    try:
        os.makedirs(results_dir, exist_ok=True)
    except OSError as e:
        print(f"Error: Cannot create results directory {results_dir}: {e}")
        return

    # Load the model
    print("Loading the metabolic model...")
    try:
        if not os.path.exists(args.model):
            raise FileNotFoundError(f"Model file not found: {args.model}")
        model = read_sbml_model(args.model)
        if len(model.reactions) == 0:
            raise ValueError("Model has no reactions. Cannot proceed with analysis.")
        print(f"Model loaded successfully with {len(model.reactions)} reactions and {len(model.metabolites)} metabolites.")
    except Exception as e:
        print(f"Error loading model: {e}")
        return

    if args.chokepoints:
        print("Identifying chokepoint reactions...")
        metabolite_reactant_product_file = os.path.join(results_dir, "metabolite_reactant_product.tsv")
        reactant_product_frequencies = metabolite_reactant_product_frequency(model, metabolite_reactant_product_file)
        chokepoints = identify_chokepoint_reactions_from_df(reactant_product_frequencies, model, results_dir)
        if not chokepoints.empty:
            mapped_genes = map_genes_to_chokepoints(chokepoints, model, results_dir)
            print("Chokepoint files saved.")
        else:
            print("No chokepoint reactions found. Skipping gene mapping.")

    if args.graph:

        print("Calculating metabolite frequencies...")
        metabolite_file = os.path.join(results_dir, "all_metabolite_frequencies.tsv")
        metabolite_frequencies = total_metabolite_frequency(model, metabolite_file)
        print("Metabolite frequencies saved.")

        metabolite_reactions = metabolite_reaction_mapper(model)
        
        # Filter selected metabolites if modified_file is provided
        print("Filtering selected metabolites...")
        if args.frequency_filter_file:
            if not os.path.exists(args.frequency_filter_file):
                print(f"Warning: Frequency filter file not found: {args.frequency_filter_file}")
                print("Using default frequency filter (20 reactions) instead.")
                selected_metabolites = filter_selected_metabolites(metabolite_frequencies)
            else:
                try:
                    selected_metabolites = filter_selected_metabolites(metabolite_frequencies, args.frequency_filter_file)
                except Exception as e:
                    print(f"Error reading frequency filter file: {e}")
                    print("Using default frequency filter (20 reactions) instead.")
                    selected_metabolites = filter_selected_metabolites(metabolite_frequencies)
        else:
            selected_metabolites = filter_selected_metabolites(metabolite_frequencies)
            note0 = ["Note:"]
            note1 = ["No curated metabolite file provided. Using default frequency filter (20 reactions)."]
            note2 = ["If you want to select certain metabolites, modify the 'all_metabolite_frequencies.tsv' file to keep only the desired metabolites."]
            note3 = ["Once modified, provide the file as an argument to the --frequency_filter_file option."]
            print_box_dynamic(note0+note1+note2+note3)

        print("Creating metabolic network graph...")
        G = metabolism_to_graph(model, selected_metabolites, metabolite_reactions)

        # Save the graph as a visualization or file
        network_file = os.path.join(results_dir, "metabolic_network")
        write_sif(G, network_file + ".sif")
        write_html(G, network_file + ".html")
        print("Complete metabolic network graph saved.")

        print("Calculating betweenness centrality...")
        result = calculate_betweenness_centrality(G, model, results_dir)
        
        if result is not None:
            subgraph, betweenness = result
            # Save the largest connected component of the graph as a visualization or file
            network_filter_file = os.path.join(results_dir, "metabolic_network_filter")
            write_sif(subgraph, network_filter_file + ".sif")
            write_html(subgraph, network_filter_file + ".html")
            print("Largest component of metabolic network graph filtered and saved.")
        else:
            print("Skipping filtered network output due to empty graph.")

if __name__ == "__main__":
    main()