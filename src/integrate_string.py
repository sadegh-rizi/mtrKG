import requests
import time
from rdflib import Literal, RDF, URIRef
from rdflib.namespace import RDFS, XSD
# Import your specific namespaces
from src.schema_definition import MTR, BIOLINK, PROV, HGNC

STRING_API_URL = "https://string-db.org/api/json/network"


def enrich_graph_with_string(g, min_score=0.4):
    """
    Scans the graph for causal Genes, queries the STRING database for
    protein-protein interactions among them, and adds the network to the graph.

    Args:
        g: The rdflib Graph
        min_score: Minimum STRING confidence score (0.0 to 1.0). 0.4 is medium confidence.
    """
    print("\n--- Starting STRING Database Integration ---")

    # 1. Gather all unique Gene symbols currently in your graph
    gene_nodes = list(g.subjects(RDF.type, MTR.Gene))
    gene_symbols = [str(node).split('/')[-1] for node in gene_nodes]

    if not gene_symbols:
        print("No genes found in the graph. Skipping STRING integration.")
        return g

    print(f"Found {len(gene_symbols)} genes. Fetching interaction network...")

    # 2. Prepare the STRING API request
    # STRING allows you to query multiple genes at once!
    params = {
        "identifiers": "%0d".join(gene_symbols),  # Carriage return separated list
        "species": 9606,  # NCBI taxonomy ID for Homo sapiens
        "required_score": int(min_score * 1000),  # STRING API expects 0-1000
        "caller_identity": "unil_rqtl_knowledge_graph_project"  # Good practice for APIs
    }

    try:
        response = requests.post(STRING_API_URL, data=params)
        if response.status_code == 200:
            interactions = response.json()
            print(f"Successfully retrieved {len(interactions)} interactions from STRING.")

            # 3. Parse the JSON and add triples to the graph
            for interaction in interactions:
                # Get the standard gene names (e.g., FADS1, THEM4)
                gene_a_symbol = interaction['preferredName_A']
                gene_b_symbol = interaction['preferredName_B']
                score = interaction['score']  # Comes back as a float 0.0 - 1.0

                # Create the node URIs
                gene_a_node = URIRef(f"http://identifiers.org/hgnc.symbol/{gene_a_symbol}")
                gene_b_node = URIRef(f"http://identifiers.org/hgnc.symbol/{gene_b_symbol}")

                # Create a unique Association node for the interaction
                interaction_node = MTR[f"PPI_{gene_a_symbol}_{gene_b_symbol}"]

                # Build the OBAN/Biolink Association structure
                g.add((interaction_node, RDF.type, BIOLINK.Association))
                g.add((interaction_node, BIOLINK.has_subject, gene_a_node))
                g.add((interaction_node, BIOLINK.has_object, gene_b_node))
                g.add((interaction_node, BIOLINK.category, BIOLINK.ProteinProteinInteraction))

                # Add the STRING confidence score and provenance
                g.add((interaction_node, MTR.interaction_score, Literal(score, datatype=XSD.float)))
                g.add((interaction_node, PROV.wasDerivedFrom, URIRef("https://string-db.org/")))

                # For direct network traversal, add a direct edge too
                g.add((gene_a_node, BIOLINK.interacts_with, gene_b_node))
                g.add((gene_b_node, BIOLINK.interacts_with, gene_a_node))  # PPIs are bidirectional

        else:
            print(f"STRING API returned an error: {response.status_code}")

    except Exception as e:
        print(f"STRING API Request failed: {e}")

    print("--- STRING Integration Complete ---\n")
    return g