import pandas as pd
import networkx as nx
import cobra.io
from cobra.io import read_sbml_model
from pyvis.network import Network
import zipfile
from collections import defaultdict

def get_locus_tags_for_reaction(model, reaction_id):
    """
    Obtain the locus tags associated with the genes in a reaction.

    :param model: COBRA model.
    :param reaction_id: Reaction IDs.

    :returns List of locus tags associated with the reaction.
    """
    locus_tags = {}
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
        selected = df_modified["Metabolite"].tolist()

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


def identify_chokepoint_reactions_from_df(df_metabolites, model):
    """
    Identify chokepoint reactions in the metabolic model based on the frequency of metabolites as reactants and products.
    If a metabolite only appears as a reactant in a single reaction, the reaction is considered a consumption chokepoint.
    If a metabolite only appears as a product in a single reaction, the reaction is considered a production chokepoint.

    :param df_metabolites (pd.DataFrame): DataFrame with metabolite frequency as reactants and products. Must contain columns ["Metabolite", "Reactant", "Product"].
    :param model (cobra.Model): COBRApy metabolic model.

    :returns pd.DataFrame: DataFrame listing chokepoint reactions categorized as production or consumption.
    """

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
                "Gene": reaction.genes,
                "Consumption_Chokepoint": is_consumption_chokepoint,
                "Production_Chokepoint": is_production_chokepoint
            })

    # Convert results to a DataFrame
    chokepoint_df = pd.DataFrame(chokepoint_reactions)
    chokepoint_df.to_csv("chokepoint_reactions.csv", index=False, sep = "\t")

    print(f"Found {len(chokepoint_df)} chokepoint reactions.")
    print(f"Consumption chokepoints: {chokepoint_df['Consumption_Chokepoint'].sum()}")
    print(f"Production chokepoints: {chokepoint_df['Production_Chokepoint'].sum()}")
    
    return chokepoint_df

def map_genes_to_chokepoints(chokepoint_df, model):
    """
    Map genes to chokepoint reactions in the metabolic model.

    :param chokepoint_df: DataFrame with chokepoint reactions.
    :param model: COBRA model.

    :returns DataFrame with genes associated with chokepoint reactions.
    """
    data = []

    for _, row in chokepoint_df.iterrows():
        for locus_tag in row["Gene"]:
            data.append({
                "Gene": locus_tag,
                "Reaction": row["Reaction"],
                "EC_Number": model.reactions.get_by_id(row["Reaction"]).annotation.get("ec-code", None),
                "KEGG_id": model.reactions.get_by_id(row["Reaction"]).annotation.get("kegg.reaction", None),
                "Rhea_id": model.reactions.get_by_id(row["Reaction"]).annotation.get("rhea", None),
                "Biocyc_id": model.reactions.get_by_id(row["Reaction"]).annotation.get("biocyc", None),
                "Consumption_Chokepoint": row["Consumption_Chokepoint"],
                "Production_Chokepoint": row["Production_Chokepoint"]
            })

    gene_df = pd.DataFrame(data)
    gene_df.to_csv("chokepoint_genes.csv", index=False, sep = "\t")

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
        if reaction.lower_bound < 0:
            products = {met.id for met, coeff in reaction.metabolites.items()}
        else:
            products = {met.id for met, coeff in reaction.metabolites.items() if coeff > 0}

        for product in products:
            if product not in excluded_metabolites:
                reactant_rxn = metabolite_reactions[product]['Reactant']
                for rxn in reactant_rxn:
                    if rxn != reaction.id:
                        graph.add_edge(reaction.id, rxn)

    return graph

def betweenness_centrality(graph, model):
    """
    Calculate the betweenness centrality of reactions in the largest connected component of the metabolic network.
    The betweenness centrality is normalized by the maximum value in the network.

    :param graph: NetworkX graph representing the metabolic network.
    :param model: COBRA model.

    :returns Tuple with the largest connected component of the graph and a dictionary with normalized betweenness centrality values.
    """

    largest_component = max(nx.weakly_connected_components(G), key=len)
    subgraph = G.subgraph(largest_component).copy()

    # Calculate centrality measures
    betweenness_centrality = nx.betweenness_centrality(subgraph)
    max_value = max(value for value in betweenness_centrality.values())
    normalized_betweenness_centrality = {key: (value / max_value) for key, value in betweenness_centrality.items()}

    # Calculate node degrees
    node_degrees = dict(subgraph.degree())

    # Create an empty list to store the data for DataFrame
    data = []

    # Iterate over reactions and their associated genes
    for reaction_id in betweenness_centrality:
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
                'Betweenness_Centrality': betweenness_centrality[reaction_id],
                'Degree': node_degrees[reaction_id]
            })

    # Create the DataFrame from the list of dictionaries
    df = pd.DataFrame(data)

    # Save the DataFrame to a file
    df.to_csv("betweenness_centrality.csv", index=False, sep = "\t")

    return subgraph, normalized_betweenness_centrality


def write_sif(graph, file_path):
    """
    Write a networkx graph to a SIF file.

    :param graph: Networkx graph.
    :param file_path: Path to the SIF file.
    """
    with open(file_path, 'w') as file:
        for edge in graph.edges():
            file.write(f"{edge[0]}\tinteracts\t{edge[1]}\n")

def write_html(graph, file_path):
    """
    Write a networkx graph to a HTML file.

    :param graph: Networkx graph.
    :param file_path: Path to the HTML file.
    """
    net = Network("500px", "500px")
    net.from_nx(graph)
    net.save_graph(file_path)


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
    args = parser.parse_args()

    # Load the model
    print("Loading the metabolic model...")
    model = read_sbml_model(args.model)

    if args.chokepoints:
        print("Identifying chokepoint reactions...")
        reactant_product_frequencies = metabolite_reactant_product_frequency(model, f"{args.output}_reactant_product.csv")
        chokepoints = identify_chokepoint_reactions_from_df(reactant_product_frequencies, model)
        mapped_genes = map_genes_to_chokepoints(chokepoints, model)
        print("Chokepoint files saved.")

    if args.graph:

        print("Calculating metabolite frequencies...")
        metabolite_frequencies = total_metabolite_frequency(model, f"all_metabolite_frequencies.csv")
        print("Metabolite frequencies saved.")

        metabolite_reactions = metabolite_reaction_mapper(model)
        
        # Filter selected metabolites if modified_file is provided
        print("Filtering selected metabolites...")
        if args.modified_file:
            selected_metabolites = filter_selected_metabolites(metabolite_frequencies, args.modified_file)
        else:
            selected_metabolites = filter_selected_metabolites(metabolite_frequencies)
            print("No curated metabolite file provided. Using default frequency filter (20 reactions).")
            print("If you want to select certain metabolites, modify the 'all_metabolite_frequencies.csv' file to keep only the desired metabolites.")
            print("Once modified, provide the file as an argument to the --frequency_filter_file option.")

        print("Creating metabolic network graph...")
        G = metabolism_to_graph(model, selected_metabolites, metabolite_reactions)

        # Save the graph as a visualization or file
        write_sif(G, "metabolic_network.sif")
        write_html(G, "metabolic_network.html")
        print("Complete metabolic network graph saved.")

        print("Calculating betweenness centrality...")
        subgraph, betweenness = betweenness_centrality(G, model)

        write_sif(subgraph, "metabolic_network_filter.sif")
        write_html(subgraph, "metabolic_network_filter.html")
        print("Largest component of metabolic network graph filtered and saved.")