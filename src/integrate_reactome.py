import requests
import time
import logging
import pandas as pd
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
        logging.FileHandler("../output/reactome_integration.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Define the Reactome Namespace for URIs
REACTOME = Namespace("https://reactome.org/content/detail/")


# ==========================================
# 2. API FETCH FUNCTIONS
# ==========================================

def fetch_reactome_pathways(gene_symbol):
    """
    Hits the Reactome AnalysisService API to find pathways for a given GENE identifier.
    """
    # FIXED the URL string formatting here
    url = f"https://reactome.org/AnalysisService/identifier/{gene_symbol}"
    headers = {"accept": "application/json"}

    params = {
        "interactors": "false",
        "pageSize": 20,
        "page": 1,
        "sortBy": "ENTITIES_PVALUE",
        "order": "ASC",
        "resource": "TOTAL",
        "pValue": 1,
        "includeDisease": "true",
        "importableOnly": "false"
    }

    try:
        response = requests.get(url, headers=headers, params=params)
        if response.status_code == 200:
            return response.json()
        elif response.status_code == 404:
            return None  # Gene not found in Reactome
        else:
            logger.warning(f"Unexpected status {response.status_code} for {gene_symbol}")
            return None
    except requests.exceptions.RequestException as e:
        logger.error(f"Request exception for {gene_symbol}: {e}")
        return None


def fetch_reactome_pathways_for_metabolite(chebi_num):
    """
    Hits the Reactome ContentService API to find pathways for a given ChEBI ID.
    """
    # Uses the specific ContentService mapping endpoint you provided
    url = f"https://reactome.org/ContentService/data/mapping/ChEBI/{chebi_num}/pathways"
    headers = {"accept": "application/json"}

    # 9606 is the taxonomy ID for Homo sapiens
    params = {"species": "9606"}

    try:
        response = requests.get(url, headers=headers, params=params)
        if response.status_code == 200:
            return response.json()
        elif response.status_code == 404:
            return None  # Metabolite not mapped in Reactome
        else:
            logger.warning(f"Unexpected status {response.status_code} for ChEBI {chebi_num}")
            return None
    except requests.exceptions.RequestException as e:
        logger.error(f"Request exception for ChEBI {chebi_num}: {e}")
        return None


# ==========================================
# 3. GRAPH ENRICHMENT FUNCTIONS
# ==========================================

def enrich_graph_with_reactome(g):
    logger.info("--- Starting Reactome AnalysisService Integration for GENES ---")

    gene_nodes = list(g.subjects(RDF.type, MTR.Gene))
    integration_report = []

    for i, gene_node in enumerate(gene_nodes):
        gene_symbol = str(gene_node).split('/')[-1].strip()
        logger.info(f"[{i + 1}/{len(gene_nodes)}] Fetching pathways for {gene_symbol}...")

        try:
            pathways_data = fetch_reactome_pathways(gene_symbol)

            if isinstance(pathways_data, dict) and 'pathways' in pathways_data:
                pathways_list = pathways_data['pathways']
            elif isinstance(pathways_data, list):
                pathways_list = pathways_data
            else:
                pathways_list = []

            if not pathways_list:
                logger.info("   -> No pathways found.")
                integration_report.append({"Gene": gene_symbol, "Status": "Unmapped", "Pathways_Found": 0})
                time.sleep(0.1)
                continue

            logger.info(f"   -> Found {len(pathways_list)} pathways!")
            integration_report.append({"Gene": gene_symbol, "Status": "Mapped", "Pathways_Found": len(pathways_list)})

            for pw in pathways_list:
                pathway_stId = pw.get('stId')
                pathway_name = pw.get('name')  # AnalysisService uses 'name'

                if not pathway_stId or not pathway_name:
                    continue

                pathway_node = REACTOME[pathway_stId]
                g.add((pathway_node, RDF.type, BIOLINK.Pathway))
                g.add((pathway_node, RDFS.label, Literal(pathway_name)))
                g.add((gene_node, BIOLINK.participates_in, pathway_node))
                g.add((pathway_node, PROV.wasDerivedFrom, URIRef("https://reactome.org/")))

        except Exception as e:
            logger.error(f"   -> ERROR processing {gene_symbol}: {e}")
            integration_report.append({"Gene": gene_symbol, "Status": "Failed", "Pathways_Found": 0})

        time.sleep(0.15)

    logger.info("--- Reactome Gene Integration Complete ---")
    df_report = pd.DataFrame(integration_report)
    df_report.to_csv("../output/reactome_gene_mapping_report.csv", index=False)

    return g


def enrich_metabolites_with_reactome(g):
    logger.info("--- Starting Reactome ContentService Integration for METABOLITES ---")

    metabolite_nodes = list(g.subjects(RDF.type, MTR.Metabolite))
    integration_report = []

    for i, met_node in enumerate(metabolite_nodes):
        uri_str = str(met_node)

        # We need JUST the number for the ContentService API URL
        if "CHEBI_" in uri_str:
            chebi_num = uri_str.split("CHEBI_")[-1]
        else:
            logger.info(f"   -> Skipping {uri_str}, not a standard ChEBI format.")
            continue

        logger.info(f"[{i + 1}/{len(metabolite_nodes)}] Fetching pathways for ChEBI:{chebi_num}...")

        try:
            # Use the new Metabolite-specific function
            pathways_list = fetch_reactome_pathways_for_metabolite(chebi_num)

            # ContentService natively returns a list, so we just check if it's valid
            if not pathways_list or not isinstance(pathways_list, list):
                logger.info("   -> No pathways found.")
                integration_report.append(
                    {"Metabolite": f"CHEBI:{chebi_num}", "Status": "Unmapped", "Pathways_Found": 0})
                time.sleep(0.15)
                continue

            logger.info(f"   -> Found {len(pathways_list)} pathways!")
            integration_report.append(
                {"Metabolite": f"CHEBI:{chebi_num}", "Status": "Mapped", "Pathways_Found": len(pathways_list)})

            for pw in pathways_list:
                pathway_stId = pw.get('stId')
                # IMPORTANT: ContentService uses 'displayName' instead of 'name'!
                pathway_name = pw.get('displayName')

                if not pathway_stId or not pathway_name:
                    continue

                # Create the Pathway Node
                pathway_node = REACTOME[pathway_stId]
                g.add((pathway_node, RDF.type, BIOLINK.Pathway))
                g.add((pathway_node, RDFS.label, Literal(pathway_name)))

                # Link the METABOLITE to the Pathway
                g.add((met_node, BIOLINK.participates_in, pathway_node))
                g.add((pathway_node, PROV.wasDerivedFrom, URIRef("https://reactome.org/")))

        except Exception as e:
            logger.error(f"   -> ERROR processing ChEBI:{chebi_num}: {e}")
            integration_report.append({"Metabolite": f"CHEBI:{chebi_num}", "Status": "Failed", "Pathways_Found": 0})

        time.sleep(0.1)

    logger.info("--- Reactome Metabolite Integration Complete ---")
    df_report = pd.DataFrame(integration_report)
    df_report.to_csv("../output/reactome_metabolite_mapping_report.csv", index=False)

    return g