import requests
import time
import logging
import pandas as pd
from rdflib import Literal, RDF, URIRef
from rdflib.namespace import RDFS
from src.schema_definition import *  # Ensure BIOLINK, MTR, PROV are defined

# ==========================================
# 1. SETUP LOGGING (Ensuring Trustworthiness)
# ==========================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("../output/encode_integration.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

ENSEMBL_SERVER = "https://rest.ensembl.org"


def get_snp_location(rsid):
    """Hits Ensembl to find the exact genomic coordinates of a SNP."""
    ext = f"/variation/human/{rsid}"
    try:
        r = requests.get(ENSEMBL_SERVER + ext, headers={"Content-Type": "application/json"})
        if r.status_code == 200:
            data = r.json()
            # Extract the first mapping location (e.g., "6:123456-123456")
            mappings = data.get('mappings', [])
            if mappings:
                return mappings[0].get('location')
    except requests.exceptions.RequestException:
        pass
    return None


def get_regulatory_overlap(location):
    """Hits Ensembl to find ENCODE regulatory features at a specific location."""
    ext = f"/overlap/region/human/{location}?feature=regulatory"
    try:
        r = requests.get(ENSEMBL_SERVER + ext, headers={"Content-Type": "application/json"})
        if r.status_code == 200:
            return r.json()
    except requests.exceptions.RequestException:
        pass
    return None


def enrich_graph_with_encode(g, max_snps=None):
    logger.info("--- Starting ENCODE Epigenetic Integration (via Ensembl) ---")

    snp_nodes = list(g.subjects(RDF.type, MTR.SNP))
    if max_snps:
        snp_nodes = snp_nodes[:max_snps]

    integration_report = []

    for i, snp_node in enumerate(snp_nodes):
        rsid = str(snp_node).split('/')[-1].strip()
        logger.info(f"[{i + 1}/{len(snp_nodes)}] Checking ENCODE regulatory overlaps for {rsid}...")

        try:
            # Step 1: Find where the SNP is on the genome
            location = get_snp_location(rsid)
            if not location:
                logger.info(f"   -> Could not find genomic coordinates for {rsid}.")
                integration_report.append({"SNP": rsid, "Status": "Unmapped_Location", "Features_Found": 0})
                time.sleep(0.1)
                continue

            # Step 2: Ask if any ENCODE regulatory features overlap this exact spot
            reg_features = get_regulatory_overlap(location)

            if not reg_features:
                logger.info("   -> No overlapping epigenetic features found.")
                integration_report.append({"SNP": rsid, "Status": "No_Overlap", "Features_Found": 0})
                time.sleep(0.1)
                continue

            logger.info(f"   -> Found {len(reg_features)} ENCODE regulatory features!")
            integration_report.append({"SNP": rsid, "Status": "Mapped", "Features_Found": len(reg_features)})

            # Step 3: Add the Epigenetic Data to the Knowledge Graph
            for feature in reg_features:
                reg_id = feature.get('id')  # e.g., ENSR0000012345
                feature_type = feature.get('description')  # e.g., "Enhancer", "Promoter"

                if not reg_id or not feature_type:
                    continue

                # Create the Regulatory Feature Node
                reg_node = MTR[reg_id]
                g.add((reg_node, RDF.type, BIOLINK.GenomicEntity))
                g.add((reg_node, RDFS.label, Literal(f"{feature_type} ({reg_id})")))

                # Link the SNP to the Regulatory Feature
                # Using biolink:overlaps ensures semantic accuracy
                g.add((snp_node, BIOLINK.overlaps, reg_node))

                # Add Provenance to ensure Trustworthiness
                g.add((reg_node, PROV.wasDerivedFrom, URIRef("https://www.ensembl.org/info/genome/funcgen/index.html")))

        except Exception as e:
            logger.error(f"   -> ERROR processing {rsid}: {e}")
            integration_report.append({"SNP": rsid, "Status": "Failed", "Features_Found": 0})

        # Ensembl allows 15 requests per second, but a small sleep is polite
        time.sleep(0.15)

    logger.info("--- ENCODE Integration Complete ---")

    df_report = pd.DataFrame(integration_report)
    df_report.to_csv("../output/encode_mapping_report.csv", index=False)

    return g