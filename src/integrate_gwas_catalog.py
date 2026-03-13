import logging
import re
import time
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from rdflib import Graph, Literal, RDF, URIRef
from rdflib.namespace import RDFS, XSD, SKOS
from urllib3.util.retry import Retry

try:
    from src.schema_definition import *
except ModuleNotFoundError:
    from schema_definition import *

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(OUTPUT_DIR / "gwas_integration.log", mode="w"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"


def build_http_session() -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=1.0,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset(["GET"]),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update(
        {
            "Accept": "application/json",
            "User-Agent": "rqtl-kg-integrator/1.0 (+https://metabolite-ratio-app.unil.ch/)",
        }
    )
    return session


SESSION = build_http_session()


def gwas_api_request(endpoint: str, params=None, timeout: int = 30):
    full_url = f"{BASE_URL.rstrip('/')}/{endpoint.lstrip('/')}"
    try:
        response = SESSION.get(full_url, params=params, timeout=timeout)
        if response.status_code == 404:
            return None
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as exc:
        logger.warning("GWAS API request failed for %s: %s", full_url, exc)
        return None


def parse_gwas_beta(beta_value):
    """
    Parses a GWAS Catalog beta string like:
      '0.61413 unit decrease'
    Returns:
      (numeric_beta: float|None, beta_unit: str|None)
    """
    if beta_value in [None, "", "NA", "NaN"]:
        return None, None

    if isinstance(beta_value, (int, float)):
        return float(beta_value), None

    beta_str = str(beta_value).strip()
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", beta_str)

    if not match:
        return None, beta_str

    numeric_val = float(match.group())
    unit_text = beta_str.replace(match.group(), "").strip()

    if "decrease" in unit_text.lower() and numeric_val > 0:
        numeric_val = -numeric_val

    return numeric_val, unit_text or None


def safe_float(value):
    if value in [None, "", "NA", "NaN"]:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def safe_literal(graph, subject, predicate, value, datatype=None):
    if value in [None, "", "NA", "NaN"]:
        return
    try:
        if datatype is not None:
            graph.add((subject, predicate, Literal(value, datatype=datatype)))
        else:
            graph.add((subject, predicate, Literal(value)))
    except Exception:
        logger.debug("Skipped invalid literal for %s %s %r", subject, predicate, value)


def rsid_from_uri(node) -> str:
    return str(node).rstrip("/").split("/")[-1].strip()


def mint_gwas_run_node():
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return MTR[f"activity/gwas_catalog_enrichment/{ts}"]


def mint_local_trait_node(trait_name: str):
    safe_trait = urllib.parse.quote(trait_name.strip().replace(" ", "_").replace('"', ""))
    return MTR[f"gwas_trait/{safe_trait}"]


def mint_ontology_term(term_id: str | None, fallback_label: str):
    """
    Mints the correct URI for ontology-backed traits.
    """
    if not term_id:
        return mint_local_trait_node(fallback_label)

    term_id = str(term_id).strip()

    # Native EFO IDs in the form EFO_0000001
    if term_id.startswith("EFO_"):
        return URIRef(f"http://www.ebi.ac.uk/efo/{term_id}")

    # CURIE style IDs like MONDO:0005148
    if ":" in term_id:
        normalized = term_id.replace(":", "_")
        return URIRef(f"http://purl.obolibrary.org/obo/{normalized}")

    # OBO-style underscore IDs like OBA_123, HP_123, MONDO_123
    if "_" in term_id:
        prefix = term_id.split("_", 1)[0]
        if prefix in {"OBA", "HP", "MONDO", "GO", "UBERON", "CHEBI", "CL", "NCIT", "ORPHA", "DOID"}:
            return URIRef(f"http://purl.obolibrary.org/obo/{term_id}")

    return mint_local_trait_node(fallback_label)


def make_trait_node(trait_dict: dict, fallback_label: str):
    """
    Returns:
      trait_node, trait_label, term_id
    """
    if not isinstance(trait_dict, dict):
        return mint_local_trait_node(fallback_label), fallback_label, None

    label = trait_dict.get("efo_trait") or fallback_label

    # Some GWAS responses include a direct URI
    uri_val = trait_dict.get("uri")
    if uri_val:
        return URIRef(uri_val), label, uri_val

    # Some responses include an ontology short form / ID
    term_id = trait_dict.get("efo_id") or trait_dict.get("short_form")
    if term_id:
        return mint_ontology_term(term_id, label), label, str(term_id)

    return mint_local_trait_node(label), label, None


def classify_trait(term_id: str | None, label: str | None):
    """
    Returns:
      (specific_type, parent_type, assoc_type)
    """
    tid = (term_id or "").upper()
    lbl = (label or "").lower()

    disease_keywords = [
        "disease", "syndrome", "disorder", "cancer", "diabetes",
        "arthritis", "infection", "stroke", "hypertension",
    ]
    phenotype_keywords = [
        "measurement", "level", "levels", "amount", "ratio", "percentage",
        "concentration", "trait", "biomarker", "cholesterol", "lipid",
        "metabolite", "protein", "urate", "glucose", "vitamin",
    ]

    if tid.startswith(("MONDO", "DOID", "ORPHA")) or any(k in lbl for k in disease_keywords):
        return (
            BIOLINK.Disease,
            BIOLINK.DiseaseOrPhenotypicFeature,
            BIOLINK.VariantToDiseaseAssociation,
        )

    if tid.startswith(("OBA", "HP", "EFO")) or any(k in lbl for k in phenotype_keywords):
        return (
            BIOLINK.PhenotypicFeature,
            BIOLINK.DiseaseOrPhenotypicFeature,
            BIOLINK.VariantToPhenotypicFeatureAssociation,
        )

    return (
        BIOLINK.DiseaseOrPhenotypicFeature,
        BIOLINK.DiseaseOrPhenotypicFeature,
        BIOLINK.Association,
    )


def ensure_node(graph, node, rdf_type=None, label=None):
    if rdf_type is not None:
        graph.add((node, RDF.type, rdf_type))
    if label:
        graph.set((node, RDFS.label, Literal(label, datatype=XSD.string)))


def add_publication(graph, assoc_node, pubmed_id):
    if not pubmed_id:
        return

    paper_node = PUBMED[str(pubmed_id)]
    graph.add((paper_node, RDF.type, BIOLINK.JournalArticle))
    graph.add((paper_node, RDFS.label, Literal(f"PubMed {pubmed_id}", datatype=XSD.string)))
    graph.add((assoc_node, BIOLINK.publications, paper_node))


def add_study(graph, assoc_node, accession_id):
    if not accession_id:
        return

    study_node = URIRef(f"https://www.ebi.ac.uk/gwas/studies/{accession_id}")
    graph.add((study_node, RDF.type, BIOLINK.Study))
    graph.add((study_node, RDFS.label, Literal(f"GWAS Catalog Study {accession_id}", datatype=XSD.string)))
    graph.add((assoc_node, PROV.wasDerivedFrom, study_node))


def enrich_graph_with_gwas(g: Graph, max_snps: int = 5, page_size: int = 10, sleep_s: float = 0.15) -> Graph:
    logger.info("--- Starting GWAS Catalog Integration (v2 API) ---")

    snp_nodes = set(g.subjects(RDF.type, MTR.SNP)) | set(g.subjects(RDF.type, BIOLINK.SequenceVariant))
    snps_to_process = list(snp_nodes)[:max_snps]

    integration_report = []

    run_node = mint_gwas_run_node()
    g.add((run_node, RDF.type, PROV.Activity))
    g.add((run_node, RDFS.label, Literal("GWAS Catalog enrichment run", datatype=XSD.string)))

    for i, snp_node in enumerate(snps_to_process):
        rsid = rsid_from_uri(snp_node)
        logger.info("[%d/%d] Fetching GWAS associations for %s...", i + 1, len(snps_to_process), rsid)

        params = {
            "rs_id": rsid,
            "sort": "p_value",
            "direction": "asc",
            "size": page_size,
        }

        try:
            associations_data = gwas_api_request("v2/associations", params=params)

            if not associations_data:
                logger.info("   -> No data returned.")
                integration_report.append(
                    {"SNP": rsid, "Status": "Unmapped", "Associations_Found": 0, "Error_Message": "No data returned"}
                )
                time.sleep(sleep_s)
                continue

            associations_list = (
                associations_data.get("_embedded", {}).get("associations", [])
                if isinstance(associations_data, dict)
                else []
            )

            if not associations_list:
                logger.info("   -> No associated traits found.")
                integration_report.append(
                    {"SNP": rsid, "Status": "Unmapped", "Associations_Found": 0, "Error_Message": "Empty list"}
                )
                time.sleep(sleep_s)
                continue

            logger.info("   -> Found %d associations.", len(associations_list))
            integration_report.append(
                {"SNP": rsid, "Status": "Mapped", "Associations_Found": len(associations_list), "Error_Message": None}
            )

            for assoc in associations_list:
                p_value = safe_float(assoc.get("p_value"))
                beta_val = assoc.get("beta")
                ci_lower = safe_float(assoc.get("ci_lower"))
                ci_upper = safe_float(assoc.get("ci_upper"))
                effect_allele = assoc.get("snp_effect_allele")

                pubmed_id = assoc.get("pubmed_id")
                accession_id = assoc.get("accession_id")
                reported_trait = assoc.get("reported_trait")
                efo_traits = assoc.get("efo_traits") or []

                if not efo_traits and reported_trait:
                    efo_traits = [{"efo_trait": reported_trait}]

                for trait_dict in efo_traits:
                    trait_node, trait_label, term_id = make_trait_node(
                        trait_dict,
                        reported_trait or "Unknown trait"
                    )

                    specific_type, parent_type, assoc_type = classify_trait(term_id, trait_label)

                    ensure_node(g, trait_node, specific_type, trait_label)
                    g.add((trait_node, RDF.type, parent_type))

                    if reported_trait and reported_trait != trait_label:
                        g.add((trait_node, SKOS.altLabel, Literal(reported_trait, datatype=XSD.string)))

                    assoc_local_id = urllib.parse.quote(
                        f"{rsid}__{accession_id or 'NA'}__{trait_label}".replace(" ", "_").replace('"', "")
                    )
                    assoc_node = MTR[f"gwas_association/{assoc_local_id}"]

                    g.add((assoc_node, RDF.type, assoc_type))
                    g.add((assoc_node, BIOLINK.subject, snp_node))
                    g.add((assoc_node, BIOLINK.object, trait_node))
                    g.add((assoc_node, BIOLINK.predicate, BIOLINK.genetically_associated_with))

                    if p_value is not None:
                        g.add((assoc_node, BIOLINK.p_value, Literal(p_value, datatype=XSD.double)))

                    numeric_beta, beta_unit = parse_gwas_beta(beta_val)
                    if numeric_beta is not None:
                        g.add((assoc_node, MTR.beta, Literal(numeric_beta, datatype=XSD.double)))
                    if beta_unit:
                        g.add((assoc_node, MTR.beta_unit, Literal(beta_unit, datatype=XSD.string)))

                    if ci_lower is not None:
                        g.add((assoc_node, MTR.ci_lower, Literal(ci_lower, datatype=XSD.double)))
                    if ci_upper is not None:
                        g.add((assoc_node, MTR.ci_upper, Literal(ci_upper, datatype=XSD.double)))

                    if effect_allele:
                        g.add((assoc_node, MTR.effect_allele, Literal(effect_allele, datatype=XSD.string)))

                    g.add((assoc_node, PROV.wasGeneratedBy, run_node))
                    add_publication(g, assoc_node, pubmed_id)
                    add_study(g, assoc_node, accession_id)

        except Exception as exc:
            logger.exception("   -> ERROR processing %s: %s", rsid, exc)
            integration_report.append(
                {"SNP": rsid, "Status": "Failed", "Associations_Found": 0, "Error_Message": str(exc)}
            )

        time.sleep(sleep_s)

    logger.info("--- GWAS Integration Complete ---")

    df_report = pd.DataFrame(integration_report)
    report_path = OUTPUT_DIR / "gwas_mapping_report.csv"
    df_report.to_csv(report_path, index=False)
    logger.info("Saved mapping report to %s", report_path)

    return g