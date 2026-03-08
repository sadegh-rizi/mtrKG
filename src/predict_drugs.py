import pandas as pd
import torch
from rdflib import Graph, URIRef, BNode
from pykeen.triples import TriplesFactory
from pykeen.pipeline import pipeline
from pykeen.predict import predict_target

import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def prepare_kg_for_ml(ttl_file_path):
    """
    Machine Learning models cannot read text strings or numbers (Literals).
    They only understand structural links between entities (URIs).
    This function strips out the Literals and prepares a clean edge-list.
    """
    logger.info("Loading Knowledge Graph for Machine Learning...")
    g = Graph()
    g.parse(ttl_file_path, format="turtle")

    # 1. Extract only Entity-to-Entity triples (ignore strings, p-values, etc.)
    triples = []
    drugs = set()
    diseases = set()

    # Namespaces for filtering later
    BIOLINK_DRUG = URIRef("https://w3id.org/biolink/vocab/Drug")
    BIOLINK_DISEASE = URIRef("https://w3id.org/biolink/vocab/Disease")

    for s, p, o in g:
        # Skip literal values and Blank Nodes
        if isinstance(o, (URIRef, BNode)) and isinstance(s, (URIRef, BNode)):
            triples.append([str(s), str(p), str(o)])

            # Catalog our Drugs and Diseases for the prediction phase
            if o == BIOLINK_DRUG and p == URIRef("http://www.w3.org/1999/02/22-rdf-syntax-ns#type"):
                drugs.add(str(s))
            if o == BIOLINK_DISEASE and p == URIRef("http://www.w3.org/1999/02/22-rdf-syntax-ns#type"):
                diseases.add(str(s))

    df_triples = pd.DataFrame(triples, columns=["head", "relation", "tail"])
    logger.info(f"Extracted {len(df_triples)} structural edges.")
    logger.info(f"Found {len(drugs)} Drugs and {len(diseases)} Diseases.")

    return df_triples, list(drugs), list(diseases)


def run_drug_repurposing(ttl_file_path="../data/mtrKG.ttl"):
    # 1. Prepare Data
    df_triples, drugs_list, diseases_list = prepare_kg_for_ml(ttl_file_path)

    if len(drugs_list) == 0 or len(diseases_list) == 0:
        logger.error("You need both Drugs and Diseases in your graph to predict links between them!")
        return

    # Convert to PyKEEN TriplesFactory format
    tf = TriplesFactory.from_labeled_triples(df_triples.values)

    # Split into Training and Testing sets (80% train, 20% test)
    training, testing = tf.split(random_state=42)

    # ==========================================
    # 2. TRAIN THE A.I. MODEL (TransE)
    # ==========================================
    logger.info("Training TransE Embedding Model (This may take a few minutes)...")

    # Note: epochs=30 is just for a quick test.
    # For your final presentation, change epochs to 200 for accurate predictions!
    result = pipeline(
        training=training,
        testing=testing,
        model='TransE',  # The classic translation embedding model
        epochs=30,  # Increase this for better accuracy!
        device='cpu'  # Change to 'cuda' if you have an Nvidia GPU
    )

    # Print the model's accuracy on the hidden 20% of edges
    # Print the model's accuracy safely
    logger.info("Model Training & Evaluation Complete!")
    try:
        # Safely fetch MRR using PyKEEN's native built-in getter
        mrr = result.metric_results.get_metric('mean_reciprocal_rank')
        logger.info(f"Mean Reciprocal Rank (MRR): {mrr:.4f}")
    except Exception as e:
        logger.warning("Could not print MRR, but moving directly to predictions anyway!")
    # ==========================================
    # 3. DRUG REPURPOSING (PREDICT NOVEL LINKS)
    # ==========================================
    logger.info("Predicting novel Drug-to-Disease indications...")
    model = result.model

    # We want to predict: (Drug, biolink:treats, Disease)
    # Even if biolink:treats is rare in your graph, the model infers it from topology!
    TREATS_RELATION = "https://w3id.org/biolink/vocab/treats"

    predictions = []

    # For every disease, what drug is most likely to treat it based on the network structure?
    for disease_uri in diseases_list:
        try:
            # Ask the model to rank all possible head nodes (Drugs) for this disease
            pred_df = predict_target(
                model=model,
                relation=TREATS_RELATION,
                tail=disease_uri,
                triples_factory=result.training,
            ).df

            # Filter the top predictions to ONLY include our known Drugs
            drug_preds = pred_df[pred_df['head_label'].isin(drugs_list)].copy()

            # Get the #1 most probable drug for this disease
            if not drug_preds.empty:
                top_drug = drug_preds.iloc[0]
                predictions.append({
                    "Disease": disease_uri.split('/')[-1],
                    "Predicted_Drug": top_drug['head_label'].split('/')[-1],
                    "Confidence_Score": top_drug['score']
                })
        except Exception as e:
            # This happens if a node was entirely isolated and the model couldn't learn it
            continue

    # Save the predictions to a CSV
    df_predictions = pd.DataFrame(predictions)

    # Sort by highest confidence
    df_predictions = df_predictions.sort_values(by="Confidence_Score", ascending=False)
    df_predictions.to_csv("../output/drug_repurposing_predictions.csv", index=False)

    logger.info("Top 5 Novel Drug Repurposing Candidates:")
    print(df_predictions.head(5))

    return df_predictions

