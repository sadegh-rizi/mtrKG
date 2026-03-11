import requests
import time
import urllib.parse
import re
import logging
from pathlib import Path
import pandas as pd
from rdflib import Graph, Literal, RDF, URIRef, Namespace
from rdflib.namespace import RDFS, XSD

try:
    from src.schema_definition import *  # Ensure BIOLINK, MTR, EFO, PROV, SKOS, PUBMED defined
except ModuleNotFoundError:
    from schema_definition import *

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ==========================================
# 1. SETUP LOGGING
# ==========================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(OUTPUT_DIR / "gwas_integration.log", mode='w'),  # Saves to file
        logging.StreamHandler()  # Prints to Jupyter cell
    ]
)
logger = logging.getLogger(__name__)

BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"


def gwas_api_request(endpoint, params=None):
    full_url = f"{BASE_URL.rstrip('/')}/{endpoint.lstrip('/')}"
    try:
        response = requests.get(full_url, params=params, headers={"Accept": "application/json"})
        if response.status_code == 200:
            return response.json()
        elif response.status_code == 404:
            return None
        else:
            return None
    except requests.exceptions.RequestException:
        return None


def parse_gwas_beta(beta_string):
    """
    Parses a GWAS Catalog beta string (e.g., "0.61413 unit decrease")
    Returns a tuple: (numeric_float, unit_description_string)
    """
    if not beta_string:
        return None, None

    beta_str = str(beta_string).strip()
    match = re.search(r"[-+]?\d*\.?\d+", beta_str)

    if match:
        numeric_val = float(match.group())
        description = beta_str.replace(match.group(), "").strip()

        if "decrease" in description.lower() and numeric_val > 0:
            numeric_val = -numeric_val

        return numeric_val, description

    return None, beta_str


def enrich_graph_with_gwas(g, max_snps=5):
    logger.info("--- Starting GWAS Catalog Integration (v2 API) ---")

    snp_nodes = list(g.subjects(RDF.type, MTR.SNP))
    snps_to_process = snp_nodes[:max_snps]

    # ==========================================
    # 2. INITIALIZE TRACKING LIST
    # ==========================================
    integration_report = []

    # Keywords to distinguish molecular traits from clinical diseases
    molecular_keywords = [
        'measurement', 'level', 'amount', 'ratio', 'percentage',
        'acid', 'protein', 'metabolite', 'lipid', 'cholesterol',
        ' to ', 'vitamin', 'peptide', 'urate', 'carnitine'
    ]

    for i, snp_node in enumerate(snps_to_process):
        rsid = str(snp_node).split('/')[-1].strip()
        logger.info(f"[{i + 1}/{len(snps_to_process)}] Fetching v2 associations for {rsid}...")

        endpoint = "v2/associations"
        params = {
            "rs_id": rsid,
            "sort": "p_value",
            "direction": "asc",
            "size": 10
        }

        try:
            associations_data = gwas_api_request(endpoint, params)

            # Check if data exists
            if not associations_data or '_embedded' not in associations_data or 'associations' not in associations_data[
                '_embedded']:
                logger.info("   -> No associated traits found in the GWAS Catalog for this SNP.")
                integration_report.append(
                    {"SNP": rsid, "Status": "Unmapped", "Associations_Found": 0, "Error_Message": "No data returned"})
                time.sleep(0.1)
                continue

            associations_list = associations_data['_embedded']['associations']

            if len(associations_list) == 0:
                logger.info("   -> No associated traits found in the GWAS Catalog for this SNP.")
                integration_report.append(
                    {"SNP": rsid, "Status": "Unmapped", "Associations_Found": 0, "Error_Message": "Empty list"})
                time.sleep(0.1)
                continue

            logger.info(f"   -> Found {len(associations_list)} significant associations!")
            integration_report.append(
                {"SNP": rsid, "Status": "Mapped", "Associations_Found": len(associations_list), "Error_Message": None})

            for assoc in associations_list:
                p_value = assoc.get('p_value')
                beta = assoc.get('beta')
                ci_lower = assoc.get('ci_lower')
                ci_upper = assoc.get('ci_upper')
                effect_allele = assoc.get('snp_effect_allele')

                pubmed_id = assoc.get('pubmed_id')
                accession_id = assoc.get('accession_id')
                reported_trait = assoc.get('reported_trait')
                efo_traits = assoc.get('efo_traits', [])

                for trait_dict in efo_traits:
                    trait_name = trait_dict.get('efo_trait')
                    if not trait_name:
                        continue

                    safe_trait_id = urllib.parse.quote(trait_name.strip().replace(" ", "_").replace('"', ''))
                    trait_node = EFO[f"Custom_{safe_trait_id}"]

                    # ==========================================
                    # 3. SEMANTIC TYPING (Disease vs Molecular)
                    # ==========================================
                    trait_label_lower = trait_name.lower()
                    if any(kw in trait_label_lower for kw in molecular_keywords):
                        g.add((trait_node, RDF.type, BIOLINK.ClinicalMeasurement))
                    else:
                        g.add((trait_node, RDF.type, BIOLINK.Disease))

                    g.add((trait_node, RDFS.label, Literal(trait_name)))

                    if reported_trait and reported_trait != trait_name:
                        g.add((trait_node, SKOS.altLabel, Literal(reported_trait)))

                    unique_assoc_id = accession_id if accession_id else f"Assoc_{rsid}_{safe_trait_id}"
                    assoc_node = MTR[f"GWAS_{unique_assoc_id}"]

                    g.add((assoc_node, RDF.type, BIOLINK.Association))
                    g.add((assoc_node, BIOLINK.has_subject, snp_node))
                    g.add((assoc_node, BIOLINK.has_object, trait_node))

                    if p_value is not None:
                        g.add((assoc_node, MTR.p_value, Literal(p_value, datatype=XSD.float)))

                    if beta is not None:
                        numeric_beta, beta_unit = parse_gwas_beta(beta)
                        if numeric_beta is not None:
                            g.add((assoc_node, MTR.beta, Literal(numeric_beta, datatype=XSD.float)))
                        if beta_unit:
                            g.add((assoc_node, MTR.beta_unit, Literal(beta_unit)))

                    if ci_lower is not None and ci_upper is not None:
                        g.add((assoc_node, MTR.ci_lower, Literal(ci_lower, datatype=XSD.float)))
                        g.add((assoc_node, MTR.ci_upper, Literal(ci_upper, datatype=XSD.float)))

                    if effect_allele:
                        g.add((assoc_node, MTR.effect_allele, Literal(effect_allele)))

                    g.add((assoc_node, PROV.wasGeneratedBy, Literal("EBI_GWAS_REST_API_v2")))

                    if pubmed_id:
                        paper_node = PUBMED[str(pubmed_id)]
                        g.add((assoc_node, BIOLINK.publications, paper_node))

                    if accession_id:
                        study_node = URIRef(f"https://www.ebi.ac.uk/gwas/studies/{accession_id}")
                        g.add((assoc_node, PROV.wasDerivedFrom, study_node))

        except Exception as e:
            logger.error(f"   -> ERROR processing {rsid}: {e}")
            integration_report.append(
                {"SNP": rsid, "Status": "Failed", "Associations_Found": 0, "Error_Message": str(e)})

        time.sleep(0.1)

    logger.info("--- GWAS Integration Complete ---")

    # ==========================================
    # 4. SAVE THE REPORT
    # ==========================================
    df_report = pd.DataFrame(integration_report)
    df_report.to_csv(OUTPUT_DIR / "gwas_mapping_report.csv", index=False)
    logger.info(f"Saved mapping report to {OUTPUT_DIR / 'gwas_mapping_report.csv'}")

    return g
