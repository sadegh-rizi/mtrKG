import pandas as pd
import numpy as np
import torch
import logging
from scipy.spatial.distance import cdist
from rdflib import Graph, URIRef, BNode
from pykeen.triples import TriplesFactory
from pykeen.pipeline import pipeline

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def prepare_kg_for_ml(ttl_file_path):
    """
    Strips out Literals and prepares a clean edge-list for the ML model.
    """
    logger.info("Loading Knowledge Graph for Machine Learning...")
    g = Graph()
    g.parse(ttl_file_path, format="turtle")

    triples = []
    drugs = set()
    diseases = set()

    BIOLINK_DRUG = URIRef("https://w3id.org/biolink/vocab/Drug")
    BIOLINK_DISEASE = URIRef("https://w3id.org/biolink/vocab/Disease")
    RDF_TYPE = URIRef("http://www.w3.org/1999/02/22-rdf-syntax-ns#type")

    for s, p, o in g:
        if isinstance(o, (URIRef, BNode)) and isinstance(s, (URIRef, BNode)):
            triples.append([str(s), str(p), str(o)])

            if o == BIOLINK_DRUG and p == RDF_TYPE:
                drugs.add(str(s))
            if o == BIOLINK_DISEASE and p == RDF_TYPE:
                diseases.add(str(s))

    df_triples = pd.DataFrame(triples, columns=["head", "relation", "tail"])
    logger.info(f"Extracted {len(df_triples)} structural edges.")
    logger.info(f"Found {len(drugs)} Drugs and {len(diseases)} Diseases.")

    return df_triples, list(drugs), list(diseases)


def run_drug_repurposing(ttl_file_path="../data/mtrKG.ttl"):
    # 1. Prepare Data
    df_triples, drugs_list, diseases_list = prepare_kg_for_ml(ttl_file_path)

    if not drugs_list or not diseases_list:
        logger.error("You need both Drugs and Diseases in your graph to predict links!")
        return

    tf = TriplesFactory.from_labeled_triples(df_triples.values)
    training, testing = tf.split(random_state=42)

    # ==========================================
    # 2. TRAIN THE A.I. MODEL (TransE)
    # ==========================================
    logger.info("Training TransE Embedding Model...")

    result = pipeline(
        training=training,
        testing=testing,
        model='TransE',
        epochs=200,  # Increased to 200 for presentation-quality embeddings
        random_seed=42,  # GUARANTEES REPRODUCIBILITY!
        device='cpu'
    )

    logger.info("Model Training & Evaluation Complete!")

    # ==========================================
    # 3. DRUG REPURPOSING (VECTORIZED PROXIMITY)
    # ==========================================
    logger.info("Predicting novel indications using Vector Space Proximity...")

    # Extract all embeddings
    entity_embeddings = result.model.entity_representations[0](indices=None).detach().cpu().numpy()
    entity_to_id = result.training.entity_to_id

    # Filter out entities that didn't make it into the training set
    valid_diseases = [uri for uri in diseases_list if uri in entity_to_id]
    valid_drugs = [uri for uri in drugs_list if uri in entity_to_id]

    if not valid_diseases or not valid_drugs:
        logger.error("Not enough valid entities found in the trained model.")
        return

    # Map URIs to their specific matrix indices
    disease_indices = [entity_to_id[uri] for uri in valid_diseases]
    drug_indices = [entity_to_id[uri] for uri in valid_drugs]

    # Extract sub-matrices for diseases and drugs
    disease_matrix = entity_embeddings[disease_indices]
    drug_matrix = entity_embeddings[drug_indices]

    # MATHEMATICAL OPTIMIZATION: Calculate distance of all pairs simultaneously
    # cdist returns a matrix of shape (len(diseases), len(drugs))
    distance_matrix = cdist(disease_matrix, drug_matrix, metric='euclidean')

    # For each disease, find the index of the drug with the minimum distance
    best_drug_indices = np.argmin(distance_matrix, axis=1)
    min_distances = np.min(distance_matrix, axis=1)

    # Build the final list of predictions instantly
    predictions = [
        {
            "Disease": valid_diseases[i].split('/')[-1],
            "Predicted_Drug": valid_drugs[best_drug_indices[i]].split('/')[-1],
            "Confidence_Score": -min_distances[i]  # Negative so closest = highest score
        }
        for i in range(len(valid_diseases))
    ]

    # 4. Clean and Save Data
    df_predictions = pd.DataFrame(predictions)
    df_predictions = df_predictions.sort_values(by="Confidence_Score", ascending=False)

    df_predictions.to_csv("../output/drug_repurposing_predictions.csv", index=False)

    logger.info("--- Top 5 Novel Drug Repurposing Candidates ---")
    print(df_predictions.head(5).to_markdown())

    return df_predictions

# Execute the script
# run_drug_repurposing()