import pandas as pd
import networkx as nx
import cobra.io
from cobra.io import read_sbml_model
from pyvis.network import Network
import zipfile

def calculate_metabolite_frequency(model, output_file):
    """
    Calculates the frequency of each metabolite (as a reactant or product) in the model
    and saves it to a file.

    :param model: COBRA model.
    :param output_file: File name to save the frequency table.

    :returns DataFrame with metabolites and their frequencies.
    """
    frequencies = {}
    for reaction in model.reactions:
        for metabolite in reaction.metabolites:
            frequencies[metabolite.id] = frequencies.get(metabolite.id, 0) + 1

    df_frequencies = pd.DataFrame(
        list(frequencies.items()), columns=["Metabolite", "Frequency"]
    )
    df_frequencies.sort_values(by="Frequency", ascending=False, inplace=True)
    df_frequencies.to_csv(output_file, index=False)
    return df_frequencies

def filter_selected_metabolites(modified_file):
    """
    Reads a modified file by the user to extract the list of selected metabolites.

    :param modified_file: Path to the modified file.

    :returns List of selected metabolites.
    """
    df_modified = pd.read_csv(modified_file)
    selected = df_modified["Metabolite"].tolist()
    return selected

def metabolite_reactant_product_frequency(model):
    """
    Calculates how many times each metabolite appears as a product and as a reactant, considering reversible reactions.

    :param model: COBRA model.

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

    df_frequencies = pd.DataFrame.from_dict(frequencies, orient="index").reset_index()
    df_frequencies.rename(columns={"index": "Metabolite"}, inplace=True)
    return df_frequencies

def metabolite_reaction_lists(model):
    """
    Identifies the reactions where each metabolite appears as a substrate (reactant) and as a product,
    considering reversible reactions.

    :param model: COBRA model.

    :returns Dictionary with metabolites as keys, and values as a dictionary containing lists of reactions
        where they are substrates and products.
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


def identify_chokepoint_reactions_from_df(df, model):
    """
    Identify chokepoint reactions based on the frequency DataFrame.

    :param df (pd.DataFrame): DataFrame with metabolite frequency as reactants and products.
                         Must contain columns ["Metabolite", "Reactant", "Product"].
    :param model (cobra.Model): COBRApy metabolic model.

    :returns pd.DataFrame: DataFrame listing chokepoint reactions categorized as production or consumption.
    """
    # Filter metabolites that only participate as reactants in a single reaction
    unique_consumers = df[df["Reactant"] == 1]["Metabolite"].tolist()
    # Filter metabolites that only participate as products in a single reaction
    unique_producers = df[df["Product"] == 1]["Metabolite"].tolist()

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
    return chokepoint_df

def metabolism_to_graph(model, excluded_metabolites, metabolite_reactions):
    """
    Transforms the metabolic model into a graph where nodes are reactions and edges are
    established if the product of one reaction is the substrate of another, excluding
    specific metabolites.

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

def write_sif(graph, file_path):
    """
    Write a networkx graph to a SIF file.

    :param graph: Networkx graph.
    :param file_path: Path to the SIF file.
    """
    with open(file_path, 'w') as file:
        for edge in graph.edges():
            file.write(f"{edge[0]}\tinteracts\t{edge[1]}\n")

