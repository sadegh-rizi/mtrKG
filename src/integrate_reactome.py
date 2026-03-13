import logging
import re
import time
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from rdflib import Graph, Literal, RDF, URIRef, Namespace
from rdflib.namespace import RDFS, XSD
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
        logging.FileHandler(OUTPUT_DIR / "reactome_integration.log", mode="w"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

REACTOME = Namespace("https://reactome.org/content/detail/")
REACTOME_SITE = URIRef("https://reactome.org/")
REACTOME_ANALYSIS_BASE = "https://reactome.org/AnalysisService"
REACTOME_CONTENT_BASE = "https://reactome.org/ContentService"

CHEBI_URI_RE = re.compile(r"CHEBI_(\d+)$")


def build_http_session() -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=1.0,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset(["GET", "POST"]),
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


def safe_literal(graph: Graph, subject, predicate, value, datatype=None):
    if value in [None, "", "NA", "NaN"]:
        return
    try:
        if datatype is not None:
            graph.add((subject, predicate, Literal(value, datatype=datatype)))
        else:
            graph.add((subject, predicate, Literal(value)))
    except Exception:
        logger.debug("Skipped invalid literal for %s %s %r", subject, predicate, value)


def ensure_node(graph: Graph, node, rdf_type=None, label=None):
    if rdf_type is not None:
        graph.add((node, RDF.type, rdf_type))
    if label:
        graph.set((node, RDFS.label, Literal(str(label), datatype=XSD.string)))


def get_preferred_label(graph: Graph, node):
    label = next(graph.objects(node, RDFS.label), None)
    if label:
        return str(label)
    return urllib.parse.unquote(str(node).rstrip("/").split("/")[-1])


def mint_run_node():
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return MTR[f"activity/reactome_enrichment/{ts}"]


def mint_assoc_node(prefix: str, *parts: str):
    safe = "__".join(
        urllib.parse.quote(str(p).strip().replace(" ", "_").replace('"', ""))
        for p in parts
    )
    return MTR[f"{prefix}/{safe}"]


def reactome_analysis_post(identifier_text: str, page_size: int = 20):
    """
    Current documented AnalysisService usage is POST to /identifiers/projection
    with text/plain identifier input.
    """
    url = f"{REACTOME_ANALYSIS_BASE}/identifiers/projection"
    params = {
        "pageSize": page_size,
        "page": 1,
        "sortBy": "ENTITIES_PVALUE",
        "order": "ASC",
        "resource": "TOTAL",
        "pValue": 1,
        "includeDisease": "true",
        "importableOnly": "false",
    }
    headers = {"Content-Type": "text/plain", "Accept": "application/json"}

    payload = f"#Identifiers\n{identifier_text}\n"

    try:
        response = SESSION.post(url, params=params, headers=headers, data=payload, timeout=45)
        if response.status_code == 404:
            return None
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as exc:
        logger.warning("Reactome AnalysisService request failed for %s: %s", identifier_text, exc)
        return None


def fetch_reactome_pathways_for_metabolite(chebi_num: str):
    """
    ContentService direct curated mapping for ChEBI -> pathways.
    """
    url = f"{REACTOME_CONTENT_BASE}/data/mapping/ChEBI/{chebi_num}/pathways"
    params = {"species": "9606"}  # Homo sapiens
    headers = {"Accept": "application/json"}

    try:
        response = SESSION.get(url, params=params, headers=headers, timeout=30)
        if response.status_code == 404:
            return None
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as exc:
        logger.warning("Reactome ContentService request failed for ChEBI %s: %s", chebi_num, exc)
        return None


def extract_chebi_number(node) -> str | None:
    match = CHEBI_URI_RE.search(str(node))
    if match:
        return match.group(1)
    return None


def add_reactome_provenance(graph: Graph):
    run_node = mint_run_node()
    ensure_node(graph, run_node, PROV.Activity, "Reactome enrichment run")

    analysis_source = URIRef(f"{REACTOME_ANALYSIS_BASE}/identifiers/projection")
    content_source = URIRef(f"{REACTOME_CONTENT_BASE}/data/mapping/ChEBI")
    ensure_node(graph, analysis_source, PROV.Entity, "Reactome AnalysisService identifiers/projection")
    ensure_node(graph, content_source, PROV.Entity, "Reactome ContentService ChEBI mapping")

    return run_node, analysis_source, content_source


def enrich_graph_with_reactome_genes(g: Graph, max_genes: int | None = None, sleep_s: float = 0.15) -> Graph:
    logger.info("--- Starting Reactome gene integration ---")

    gene_nodes = list(g.subjects(RDF.type, MTR.Gene))
    if max_genes is not None:
        gene_nodes = gene_nodes[:max_genes]

    run_node, analysis_source, _ = add_reactome_provenance(g)
    integration_report = []

    for i, gene_node in enumerate(gene_nodes):
        gene_symbol = get_preferred_label(g, gene_node)
        logger.info("[%d/%d] Fetching Reactome analysis pathways for %s...", i + 1, len(gene_nodes), gene_symbol)

        try:
            pathways_data = reactome_analysis_post(gene_symbol)

            pathways_list = pathways_data.get("pathways", []) if isinstance(pathways_data, dict) else []
            summary = pathways_data.get("summary", {}) if isinstance(pathways_data, dict) else {}
            analysis_token = summary.get("token")

            if not pathways_list:
                logger.info("   -> No pathways found.")
                integration_report.append(
                    {"Entity": gene_symbol, "Entity_Type": "Gene", "Status": "Unmapped", "Pathways_Found": 0}
                )
                time.sleep(sleep_s)
                continue

            logger.info("   -> Found %d pathways.", len(pathways_list))
            integration_report.append(
                {"Entity": gene_symbol, "Entity_Type": "Gene", "Status": "Mapped", "Pathways_Found": len(pathways_list)}
            )

            ensure_node(g, gene_node, MTR.Gene, gene_symbol)

            for pw in pathways_list:
                pathway_stid = pw.get("stId")
                pathway_name = pw.get("name")
                if not pathway_stid or not pathway_name:
                    continue

                pathway_node = REACTOME[pathway_stid]
                ensure_node(g, pathway_node, BIOLINK.Pathway, pathway_name)

                assoc_node = mint_assoc_node("reactome_gene_pathway", gene_symbol, pathway_stid)
                g.add((assoc_node, RDF.type, BIOLINK.GeneToPathwayAssociation))
                g.add((assoc_node, BIOLINK.subject, gene_node))
                g.add((assoc_node, BIOLINK.object, pathway_node))
                g.add((assoc_node, BIOLINK.predicate, BIOLINK.participates_in))
                g.add((assoc_node, PROV.wasGeneratedBy, run_node))
                g.add((assoc_node, PROV.wasDerivedFrom, analysis_source))

                safe_literal(g, assoc_node, MTR.reactome_source_type, "AnalysisService", XSD.string)
                safe_literal(g, assoc_node, MTR.reactome_analysis_token, analysis_token, XSD.string)
                safe_literal(g, assoc_node, MTR.reactome_species, pw.get("speciesName"), XSD.string)

                entities = pw.get("entities", {}) or {}
                reactions = pw.get("reactions", {}) or {}

                # AnalysisService statistics if present
                safe_literal(g, assoc_node, MTR.reactome_entities_p_value, entities.get("pValue"), XSD.double)
                safe_literal(g, assoc_node, MTR.reactome_entities_fdr, entities.get("fdr"), XSD.double)
                safe_literal(g, assoc_node, MTR.reactome_entities_found, entities.get("found"), XSD.integer)
                safe_literal(g, assoc_node, MTR.reactome_entities_total, entities.get("total"), XSD.integer)
                safe_literal(g, assoc_node, MTR.reactome_reactions_found, reactions.get("found"), XSD.integer)
                safe_literal(g, assoc_node, MTR.reactome_reactions_total, reactions.get("total"), XSD.integer)

        except Exception as exc:
            logger.exception("   -> ERROR processing gene %s: %s", gene_symbol, exc)
            integration_report.append(
                {"Entity": gene_symbol, "Entity_Type": "Gene", "Status": "Failed", "Pathways_Found": 0}
            )

        time.sleep(sleep_s)

    report_path = OUTPUT_DIR / "reactome_gene_mapping_report.csv"
    pd.DataFrame(integration_report).to_csv(report_path, index=False)
    logger.info("Saved gene mapping report to %s", report_path)
    logger.info("--- Reactome gene integration complete ---")
    return g


def enrich_graph_with_reactome_metabolites(g: Graph, max_metabolites: int | None = None, sleep_s: float = 0.10) -> Graph:
    logger.info("--- Starting Reactome metabolite integration ---")

    metabolite_nodes = list(g.subjects(RDF.type, MTR.Metabolite))
    if max_metabolites is not None:
        metabolite_nodes = metabolite_nodes[:max_metabolites]

    run_node, _, content_source = add_reactome_provenance(g)
    integration_report = []

    for i, met_node in enumerate(metabolite_nodes):
        chebi_num = extract_chebi_number(met_node)
        met_label = get_preferred_label(g, met_node)

        if not chebi_num:
            logger.info("[%d/%d] Skipping %s, no ChEBI numeric ID in URI.", i + 1, len(metabolite_nodes), met_label)
            integration_report.append(
                {"Entity": met_label, "Entity_Type": "Metabolite", "Status": "Skipped", "Pathways_Found": 0}
            )
            continue

        logger.info("[%d/%d] Fetching Reactome pathways for ChEBI:%s...", i + 1, len(metabolite_nodes), chebi_num)

        try:
            pathways_list = fetch_reactome_pathways_for_metabolite(chebi_num)

            if not pathways_list or not isinstance(pathways_list, list):
                logger.info("   -> No pathways found.")
                integration_report.append(
                    {"Entity": f"CHEBI:{chebi_num}", "Entity_Type": "Metabolite", "Status": "Unmapped", "Pathways_Found": 0}
                )
                time.sleep(sleep_s)
                continue

            logger.info("   -> Found %d pathways.", len(pathways_list))
            integration_report.append(
                {"Entity": f"CHEBI:{chebi_num}", "Entity_Type": "Metabolite", "Status": "Mapped", "Pathways_Found": len(pathways_list)}
            )

            for pw in pathways_list:
                pathway_stid = pw.get("stId")
                pathway_name = pw.get("displayName") or pw.get("name")
                if not pathway_stid or not pathway_name:
                    continue

                pathway_node = REACTOME[pathway_stid]
                ensure_node(g, pathway_node, BIOLINK.Pathway, pathway_name)

                assoc_node = mint_assoc_node("reactome_metabolite_pathway", chebi_num, pathway_stid)
                g.add((assoc_node, RDF.type, BIOLINK.ChemicalEntityToPathwayAssociation))
                g.add((assoc_node, BIOLINK.subject, met_node))
                g.add((assoc_node, BIOLINK.object, pathway_node))
                g.add((assoc_node, BIOLINK.predicate, BIOLINK.participates_in))
                g.add((assoc_node, PROV.wasGeneratedBy, run_node))
                g.add((assoc_node, PROV.wasDerivedFrom, content_source))

                safe_literal(g, assoc_node, MTR.reactome_source_type, "ContentService", XSD.string)
                safe_literal(g, assoc_node, MTR.reactome_species, "Homo sapiens", XSD.string)

        except Exception as exc:
            logger.exception("   -> ERROR processing ChEBI:%s: %s", chebi_num, exc)
            integration_report.append(
                {"Entity": f"CHEBI:{chebi_num}", "Entity_Type": "Metabolite", "Status": "Failed", "Pathways_Found": 0}
            )

        time.sleep(sleep_s)

    report_path = OUTPUT_DIR / "reactome_metabolite_mapping_report.csv"
    pd.DataFrame(integration_report).to_csv(report_path, index=False)
    logger.info("Saved metabolite mapping report to %s", report_path)
    logger.info("--- Reactome metabolite integration complete ---")
    return g


def enrich_graph_with_reactome(g: Graph, max_genes: int | None = None, max_metabolites: int | None = None) -> Graph:
    g = enrich_graph_with_reactome_genes(g, max_genes=max_genes)
    g = enrich_graph_with_reactome_metabolites(g, max_metabolites=max_metabolites)
    return g