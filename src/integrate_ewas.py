import requests
import time
import urllib.parse
import logging
import pandas as pd
from rdflib import Graph, Literal, RDF, URIRef
from rdflib.namespace import RDFS, XSD
from src.schema_definition import *  # Ensure BIOLINK, MTR, EFO, PROV, PUBMED are defined

# ==========================================
# 1. SETUP LOGGING
# ==========================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("../output/ewas_integration.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

EWAS_API_URL = "https://www.ewascatalog.org/api/"


def ewas_api_request(gene_symbol):
    """
    Queries the EWAS Catalog API for a specific gene symbol.
    """
    params = {"gene": gene_symbol}
    try:
        response = requests.get(EWAS_API_URL, params=params, headers={"Accept": "application/json"})
        if response.status_code == 200:
            return response.json()
        return None
    except requests.exceptions.RequestException as e:
        logger.error(f"API Request failed for {gene_symbol}: {e}")
        return None


def enrich_graph_with_ewas(g):
    logger.info("--- Starting EWAS Catalog Integration ---")

    # 1. Gather all Gene nodes from the current graph
    gene_nodes = list(g.subjects(RDF.type, MTR.Gene))

    integration_report = []

    for i, gene_node in enumerate(gene_nodes):
        # Extract the pure HGNC symbol (e.g., "FADS1")
        gene_symbol = str(gene_node).split('/')[-1].strip()
        logger.info(f"[{i + 1}/{len(gene_nodes)}] Fetching EWAS associations for Gene: {gene_symbol}...")

        try:
            ewas_data = ewas_api_request(gene_symbol)

            # The EWAS API returns a list of association dictionaries
            if not ewas_data or len(ewas_data) == 0:
                logger.info("   -> No epigenetic associations found.")
                integration_report.append({"Gene": gene_symbol, "Status": "Unmapped", "CpG_Found": 0})
                time.sleep(0.2)
                continue

            logger.info(f"   -> Found {len(ewas_data)} CpG associations!")
            integration_report.append({"Gene": gene_symbol, "Status": "Mapped", "CpG_Found": len(ewas_data)})

            for assoc in ewas_data:
                cpg_id = assoc.get('cpg')
                trait_name = assoc.get('trait')
                p_value = assoc.get('p')
                beta = assoc.get('beta')
                pmid = assoc.get('pmid')

                if not cpg_id or not trait_name:
                    continue

                # ==========================================
                # 2. GRAPH CONSTRUCTION: The CpG Site
                # ==========================================
                # Treat the CpG methylation site as a Genomic Entity
                cpg_node = MTR[cpg_id]
                g.add((cpg_node, RDF.type, BIOLINK.GenomicEntity))
                g.add((cpg_node, RDFS.label, Literal(f"CpG site {cpg_id}")))

                # Link the Gene to the CpG site (Epigenetic regulation)
                g.add((cpg_node, BIOLINK.regulates, gene_node))

                # ==========================================
                # 3. GRAPH CONSTRUCTION: The Trait
                # ==========================================
                safe_trait_id = urllib.parse.quote(trait_name.strip().replace(" ", "_").replace('"', ''))
                trait_node = EFO[f"Custom_EWAS_{safe_trait_id}"]

                g.add((trait_node, RDF.type, BIOLINK.PhenotypicFeature))
                g.add((trait_node, RDFS.label, Literal(trait_name)))

                # ==========================================
                # 4. GRAPH CONSTRUCTION: The Association
                # ==========================================
                # Create a unique Association connecting the CpG site to the Trait
                assoc_node = MTR[f"EWAS_Assoc_{cpg_id}_{safe_trait_id}"]
                g.add((assoc_node, RDF.type, BIOLINK.Association))
                g.add((assoc_node, BIOLINK.has_subject, cpg_node))
                g.add((assoc_node, BIOLINK.has_object, trait_node))

                # Add Stats
                if p_value is not None:
                    g.add((assoc_node, MTR.p_value, Literal(float(p_value), datatype=XSD.float)))
                if beta is not None:
                    g.add((assoc_node, MTR.beta, Literal(float(beta), datatype=XSD.float)))

                # Add Provenance
                g.add((assoc_node, PROV.wasGeneratedBy, Literal("EWAS_Catalog_API")))
                if pmid:
                    g.add((assoc_node, BIOLINK.publications, PUBMED[str(pmid)]))

        except Exception as e:
            logger.error(f"   -> ERROR processing {gene_symbol}: {e}")
            integration_report.append({"Gene": gene_symbol, "Status": "Failed", "CpG_Found": 0})

        # Be polite to the EWAS server
        time.sleep(0.2)

    logger.info("--- EWAS Integration Complete ---")

    # Save the tracking report
    df_report = pd.DataFrame(integration_report)
    df_report.to_csv("../output/ewas_mapping_report.csv", index=False)

    return g