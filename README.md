# MetaGraphTools

# Metabolic Model Analyzer

This Python tool is designed to analyze and visualize metabolic networks by processing SBML-format models. It calculates **chokepoint reactions**, constructs a reaction graph, and computes **betweenness centrality** to identify key reactions with high metabolic impact. The analysis helps identify potential metabolic targets with significant implications for metabolic engineering or therapeutic interventions.

## Key Features

- **Chokepoint Detection**: Identifies unique reactions that exclusively consume or produce a metabolite. This is crucial for pinpointing key reactions in the network.
  
- **Graph Construction**: Transforms the metabolic network into a graph where:
  - Each **node** represents a reaction.
  - An **edge** between nodes is created if the product of one reaction is the substrate of another.
  
- **Filtering of Ubiquitous Compounds**: Compounds that are involved in more than 20 reactions (by default) are considered ubiquitous and filtered out to avoid generating artificial interactions in the graph. Alternatively, users can provide their own list of metabolites to filter.
  
- **SIF and HTML Output**: The reaction graph is saved in SIF (Simple Interaction Format) for further analysis and visualization, and also as an HTML file using `pyvis` to allow interactive graph exploration.
  
- **Component Analysis**: Computes the **largest connected component** of the metabolic network to focus on the most significant interactions.
  
- **Betweenness Centrality Calculation**: Calculates the **normalized betweenness centrality** for each reaction, which quantifies the importance of a reaction in facilitating the flow of metabolites through the network. Reactions with high centrality are often crucial for metabolic regulation and can be key targets for intervention.
  
- **Gene-Level Analysis**: Centrality calculations are also performed for genes, providing insight into the impact of gene deletions or modifications on network dynamics.

## Usage

To run the tool, use the following command:

```bash
python3 main.py --model MODEL [--chokepoints] [--graph] [--frequency_filter_file FREQUENCY_FILTER_FILE] [--output OUTPUT]
```

# Arguments:
- `-h, --help`: Show help message and exit.
- `--model MODEL`: Path to the SBML file representing the metabolic model (required).
- `--chokepoints`: If specified, the tool identifies chokepoint reactions and saves results to a file.
- `--graph`: If specified, generates the reaction graph.
- `--frequency_filter_file FREQUENCY_FILTER_FILE`: Optional file containing metabolites to filter based on their frequency in the network (used only with `--graph`).
- `--output OUTPUT`: Path to the output directory to save the results.


# Example Commands:
1. **Analyze metabolic model without graph generation**:
```bash
python3 main.py --model /path/to/your/model.sbml
```

2. Identify Chokepoint Reactions and Save the Results:
```bash
python3 main.py --model /path/to/your/model.sbml --chokepoints --output /path/to/output/
```

3. Generate a Graph of the Metabolic Network:
```bash
python3 main.py --model /path/to/your/model.sbml --graph --output /path/to/output/
```

4. Generate a Graph with Custom Frequency Filtering:
```bash
python3 main.py --model /path/to/your/model.sbml --graph --frequency_filter_file /path/to/freq
```

# Output Files

The following files will be generated:

## Chokepoint Calculator
- **chokepoint_genes.tsv**: A list of genes associated with chokepoint reactions.
- **chokepoint_reactions.tsv**: A list of reactions identified as chokepoints.
- **metabolite_reactant_product.tsv**: A table detailing, for each metabolite, the number of reactions in which it acts as a reactant or product.

## Graph Analysis
- **all_metabolite_frequencies.tsv**: A table listing metabolites and their frequencies within the model. This table can be curated to filter out specific metabolites when generating the graph.
- **betweenness_centrality.tsv**: A table containing the normalized betweenness centrality for each reaction and gene.

## Network Files
- **metabolic_network.sif**: The original SIF file representing the complete metabolic network.
- **metabolic_network_filter.sif**: A SIF file representing the largest connected component of the metabolic network.

## Visualizations
- **metabolic_network.html**: An interactive HTML visualization of the complete metabolic network.
- **metabolic_network_filter.html**: An interactive HTML visualization of the largest connected component of the metabolic network.


