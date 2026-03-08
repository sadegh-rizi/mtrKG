import requests
import time
import logging
import pandas as pd
import concurrent.futures
from rdflib import Literal, RDF, URIRef, Namespace
from rdflib.namespace import RDFS
from src.schema_definition import *  # Ensure BIOLINK, MTR, PROV are defined

# ==========================================
# 1. SETUP LOGGING
# ==========================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("../output/rhea_integration.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Define the Rhea Namespace
RHEA = Namespace("http://rdf.rhea-db.org/")
RHEA_SPARQL_ENDPOINT = "https://sparql.rhea-db.org/sparql"


# ==========================================
# 2. API FETCH FUNCTION (SPARQL to JSON)
# ==========================================
def fetch_rhea_reactions(chebi_num):
    """
    Queries the official Rhea SPARQL endpoint to find all biochemical reactions
    that involve a specific ChEBI molecule as a substrate or product.
    """
    # Rhea's RDF structure is beautiful: Reaction -> Side -> Participant -> Compound -> ChEBI
    query = f"""
    PREFIX rh: <http://rdf.rhea-db.org/>
    PREFIX chebi: <http://purl.obolibrary.org/obo/CHEBI_>

    SELECT DISTINCT ?reaction ?equation
    WHERE {{
      ?reaction rh:equation ?equation .
      ?reaction rh:side ?side .
      ?side rh:contains ?participant .
      ?participant rh:compound ?compound .
      ?compound rh:chebi chebi:{chebi_num} .
    }}
    """

    headers = {
        "Accept": "application/sparql-results+json",
        "User-Agent": "Metabolite-Ratio-KG-Builder/1.0"
    }

    try:
        response = requests.get(RHEA_SPARQL_ENDPOINT, params={"query": query}, headers=headers, timeout=15)
        if response.status_code == 200:
            data = response.json()
            return data["results"]["bindings"]
        else:
            logger.warning(f"Rhea SPARQL returned status {response.status_code} for ChEBI:{chebi_num}")
            return []
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to fetch Rhea data for ChEBI:{chebi_num}: {e}")
        return []


# ==========================================
# 3. MULTI-THREADED GRAPH ENRICHMENT
# ==========================================
def enrich_graph_with_rhea(g, max_workers=5):
    logger.info("--- Starting Multi-Threaded Rhea Biochemistry Integration ---")

    metabolite_nodes = list(g.subjects(RDF.type, MTR.Metabolite))

    # Deduplication
    unique_chebis = {}
    for met_node in metabolite_nodes:
        uri_str = str(met_node)
        if "CHEBI_" in uri_str:
            chebi_num = uri_str.split("CHEBI_")[-1]
            if chebi_num not in unique_chebis:
                unique_chebis[chebi_num] = []
            unique_chebis[chebi_num].append(met_node)

    logger.info(f"Querying Rhea for {len(unique_chebis)} unique metabolites...")

    integration_report = []

    def fetch_task(chebi_num):
        try:
            reactions = fetch_rhea_reactions(chebi_num)
            if not reactions:
                return chebi_num, [], "Unmapped"
            return chebi_num, reactions, "Mapped"
        except Exception as e:
            return chebi_num, [], f"Failed: {e}"

    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_chebi = {executor.submit(fetch_task, chebi): chebi for chebi in unique_chebis.keys()}

        completed = 0
        for future in concurrent.futures.as_completed(future_to_chebi):
            completed += 1
            chebi_num = future_to_chebi[future]
            chebi_num_returned, reactions, status = future.result()

            results[chebi_num_returned] = {"reactions": reactions, "status": status}
            logger.info(
                f"[{completed}/{len(unique_chebis)}] Processed ChEBI:{chebi_num_returned} -> {status} ({len(reactions)} reactions)")
            time.sleep(0.1)  # Polite delay

    # Safely write to the Graph sequentially
    logger.info("Writing threaded results to the Knowledge Graph...")
    for chebi_num, data in results.items():
        status = data["status"]
        reactions = data["reactions"]

        integration_report.append({
            "Metabolite": f"CHEBI:{chebi_num}",
            "Status": status,
            "Reactions_Found": len(reactions)
        })

        if reactions:
            for rxn in reactions:
                reaction_uri_str = rxn["reaction"]["value"]
                equation = rxn["equation"]["value"]

                # Use the exact Rhea URI so it links seamlessly with external databases
                reaction_node = URIRef(reaction_uri_str)

                g.add((reaction_node, RDF.type, BIOLINK.MolecularActivity))
                g.add((reaction_node, RDFS.label, Literal(equation)))
                g.add((reaction_node, PROV.wasDerivedFrom, URIRef("https://www.rhea-db.org/")))

                # Link the Metabolite to the Reaction
                for met_node in unique_chebis[chebi_num]:
                    g.add((met_node, BIOLINK.participates_in, reaction_node))

    logger.info("--- Rhea Integration Complete ---")
    df_report = pd.DataFrame(integration_report)
    df_report.to_csv("../output/rhea_mapping_report.csv", index=False)

    return g