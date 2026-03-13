import logging
import time
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from rdflib import Graph, Literal, RDF, URIRef
from rdflib.namespace import RDFS, XSD
from urllib3.util.retry import Retry

from schema_definition import *

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(OUTPUT_DIR / "opentargets_integration.log", mode="w"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

OPENTARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"


def build_http_session() -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=1.0,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset(["POST"]),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update(
        {
            "Accept": "application/json",
            "Content-Type": "application/json",
            "User-Agent": "rqtl-kg-integrator/1.0 (+https://metabolite-ratio-app.unil.ch/)",
        }
    )
    return session


SESSION = build_http_session()


def graphql_post(query: str, variables: dict, timeout: int = 45):
    try:
        response = SESSION.post(
            OPENTARGETS_API,
            json={"query": query, "variables": variables},
            timeout=timeout,
        )
        response.raise_for_status()
        payload = response.json()

        if payload.get("errors"):
            logger.warning("Open Targets GraphQL errors: %s", payload["errors"])
            return None

        return payload.get("data")
    except requests.exceptions.RequestException as exc:
        logger.warning("Open Targets request failed: %s", exc)
        return None


def safe_id(text: str) -> str:
    return urllib.parse.quote(str(text).strip().replace(" ", "_").replace('"', ""))


def ensure_node(graph: Graph, node, rdf_type=None, label=None):
    if rdf_type is not None:
        graph.add((node, RDF.type, rdf_type))
    if label:
        graph.set((node, RDFS.label, Literal(str(label), datatype=XSD.string)))


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


def get_gene_symbol(g: Graph, gene_node) -> str:
    label = next(g.objects(gene_node, RDFS.label), None)
    if label:
        return str(label)
    return urllib.parse.unquote(str(gene_node).rstrip("/").split("/")[-1].strip())


def mint_run_node():
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return MTR[f"activity/opentargets_enrichment/{ts}"]


def mint_assoc_node(prefix: str, *parts: str):
    safe = "__".join(
        urllib.parse.quote(str(p).strip().replace(" ", "_").replace('"', ""))
        for p in parts
    )
    return MTR[f"{prefix}/{safe}"]


def obo_or_local_term(term_id: str | None, label: str):
    if not term_id:
        return MTR[f"ot_term/{safe_id(label)}"]

    term_id = str(term_id).strip()

    if term_id.startswith("EFO_"):
        return URIRef(f"http://www.ebi.ac.uk/efo/{term_id}")

    if ":" in term_id:
        normalized = term_id.replace(":", "_")
        return URIRef(f"http://purl.obolibrary.org/obo/{normalized}")

    if "_" in term_id:
        prefix = term_id.split("_", 1)[0]
        if prefix in {"MONDO", "HP", "OBA", "GO", "UBERON", "CHEBI", "CL", "NCIT", "ORPHA"}:
            return URIRef(f"http://purl.obolibrary.org/obo/{term_id}")

    return MTR[f"ot_term/{safe_id(term_id)}"]


def classify_disease_or_phenotype(term_id: str | None, label: str | None):
    tid = (term_id or "").upper()
    lbl = (label or "").lower()

    if tid.startswith(("MONDO", "DOID", "ORPHA")):
        return BIOLINK.Disease

    if tid.startswith(("HP", "OBA")):
        return BIOLINK.PhenotypicFeature

    disease_keywords = [
        "disease", "syndrome", "disorder", "cancer", "diabetes",
        "asthma", "arthritis", "infection", "stroke", "hypertension",
    ]
    phenotype_keywords = [
        "phenotype", "measurement", "level", "concentration",
        "ratio", "amount", "trait", "abnormality",
    ]

    if any(k in lbl for k in disease_keywords):
        return BIOLINK.Disease
    if any(k in lbl for k in phenotype_keywords):
        return BIOLINK.PhenotypicFeature

    return BIOLINK.DiseaseOrPhenotypicFeature


def search_target_by_symbol_exact(gene_symbol: str):
    query = """
    query resolveTarget($symbol: String!) {
      search(queryString: $symbol, entityNames: ["target"], page: {index: 0, size: 10}) {
        hits {
          id
          object {
            ... on Target {
              id
              approvedSymbol
            }
          }
        }
      }
    }
    """
    data = graphql_post(query, {"symbol": gene_symbol})
    if not data:
        return None

    hits = data.get("search", {}).get("hits", [])
    if not hits:
        return None

    for hit in hits:
        obj = hit.get("object", {})
        approved = obj.get("approvedSymbol")
        if approved and approved.upper() == gene_symbol.upper():
            return {"id": obj.get("id"), "approvedSymbol": approved}

    obj = hits[0].get("object", {})
    if obj.get("id"):
        return {"id": obj.get("id"), "approvedSymbol": obj.get("approvedSymbol")}
    return None


def fetch_target_bundle(ensembl_id: str):
    query = """
    query targetBundle($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        id
        approvedSymbol
        tractability {
          modality
          label
          value
        }
        safetyLiabilities {
          event
          eventId
          biosamples {
            tissueLabel
          }
        }
        knownDrugs(size: 10) {
          rows {
            phase
            drug {
              id
              name
            }
          }
        }
        associatedDiseases(page: {index: 0, size: 10}) {
          rows {
            disease {
              id
              name
            }
            datasourceScores {
              id
              score
            }
          }
        }
      }
    }
    """
    data = graphql_post(query, {"ensemblId": ensembl_id})
    if not data:
        return None
    return data.get("target")


def enrich_graph_with_opentargets(g: Graph, max_genes: int = 5, sleep_s: float = 0.15) -> Graph:
    logger.info("--- Starting Open Targets Integration ---")

    gene_nodes = list(g.subjects(RDF.type, MTR.Gene))
    logger.info("Found %d genes in the local graph.", len(gene_nodes))
    genes_to_process = gene_nodes[:max_genes]

    run_node = mint_run_node()
    source_node = URIRef("https://platform.opentargets.org/")
    g.add((run_node, RDF.type, PROV.Activity))
    g.add((run_node, RDFS.label, Literal("Open Targets enrichment run", datatype=XSD.string)))
    g.add((source_node, RDF.type, PROV.Entity))
    g.add((source_node, RDFS.label, Literal("Open Targets Platform", datatype=XSD.string)))

    for i, gene_node in enumerate(genes_to_process):
        gene_symbol = get_gene_symbol(g, gene_node)
        logger.info("[%d/%d] Fetching OT data for %s...", i + 1, len(genes_to_process), gene_symbol)

        resolved = search_target_by_symbol_exact(gene_symbol)
        if not resolved or not resolved.get("id"):
            logger.warning("   -> %s not found in Open Targets database.", gene_symbol)
            time.sleep(sleep_s)
            continue

        ensembl_id = resolved["id"]
        approved_symbol = resolved.get("approvedSymbol") or gene_symbol

        ensure_node(g, gene_node, MTR.Gene, approved_symbol)
        safe_literal(g, gene_node, MTR.open_targets_target_id, ensembl_id, XSD.string)

        target_obj = fetch_target_bundle(ensembl_id)
        if not target_obj:
            logger.warning("   -> Failed to fetch target bundle for %s (%s).", approved_symbol, ensembl_id)
            time.sleep(sleep_s)
            continue

        # A. Tractability
        for item in target_obj.get("tractability", []) or []:
            if item.get("value") is not True:
                continue

            modality = item.get("modality")
            if modality == "SM":
                safe_literal(g, gene_node, MTR.small_molecule_tractable, True, XSD.boolean)
            elif modality == "AB":
                safe_literal(g, gene_node, MTR.antibody_tractable, True, XSD.boolean)
            elif modality == "PR":
                safe_literal(g, gene_node, MTR.protac_tractable, True, XSD.boolean)
            elif modality == "OC":
                safe_literal(g, gene_node, MTR.other_clinical_tractable, True, XSD.boolean)

        # B. Safety liabilities
        safety_liabilities = target_obj.get("safetyLiabilities") or []
        if safety_liabilities:
            logger.info("   -> Found %d safety liabilities.", len(safety_liabilities))

        for liability in safety_liabilities:
            event_name = liability.get("event")
            event_id = liability.get("eventId")
            if not event_name:
                continue

            event_node = obo_or_local_term(event_id, event_name)
            ensure_node(g, event_node, classify_disease_or_phenotype(event_id, event_name), event_name)

            assoc_node = mint_assoc_node("ot_safety_liability", approved_symbol, event_id or event_name)
            g.add((assoc_node, RDF.type, MTR.TargetSafetyLiabilityAssociation))
            g.add((assoc_node, BIOLINK.subject, gene_node))
            g.add((assoc_node, BIOLINK.object, event_node))
            g.add((assoc_node, BIOLINK.predicate, MTR.has_safety_liability))
            g.add((assoc_node, PROV.wasGeneratedBy, run_node))
            g.add((assoc_node, PROV.wasDerivedFrom, source_node))

            if event_id:
                safe_literal(g, assoc_node, MTR.external_event_id, event_id, XSD.string)

            for bio in liability.get("biosamples", []) or []:
                tissue = bio.get("tissueLabel")
                if tissue:
                    safe_literal(g, assoc_node, MTR.tissue_context, tissue, XSD.string)

        # C. Known drugs
        drug_rows = (target_obj.get("knownDrugs") or {}).get("rows", [])
        if drug_rows:
            logger.info("   -> Found %d known drugs.", len(drug_rows))

        for row in drug_rows:
            drug_info = row.get("drug") or {}
            drug_id = drug_info.get("id")
            drug_name = drug_info.get("name")
            phase = row.get("phase")

            if not drug_id or not drug_name:
                continue

            drug_node = CHEMBL[drug_id]
            ensure_node(g, drug_node, BIOLINK.Drug, drug_name)
            if phase is not None:
                safe_literal(g, drug_node, MTR.max_clinical_phase, int(phase), XSD.integer)

            assoc_node = mint_assoc_node("ot_drug_target", drug_id, approved_symbol)
            g.add((assoc_node, RDF.type, BIOLINK.DrugToGeneAssociation))
            g.add((assoc_node, BIOLINK.subject, drug_node))
            g.add((assoc_node, BIOLINK.object, gene_node))
            g.add((assoc_node, BIOLINK.predicate, BIOLINK.interacts_with))
            g.add((assoc_node, PROV.wasGeneratedBy, run_node))
            g.add((assoc_node, PROV.wasDerivedFrom, source_node))

            if phase is not None:
                safe_literal(g, assoc_node, MTR.max_clinical_phase, int(phase), XSD.integer)

        # D. Associated diseases / phenotypes
        disease_rows = (target_obj.get("associatedDiseases") or {}).get("rows", [])
        if disease_rows:
            logger.info("   -> Found %d associated disease/phenotype rows.", len(disease_rows))

        for row in disease_rows:
            disease_info = row.get("disease") or {}
            term_id = disease_info.get("id")
            term_name = disease_info.get("name")

            if not term_name:
                continue

            term_node = obo_or_local_term(term_id, term_name)
            term_type = classify_disease_or_phenotype(term_id, term_name)
            ensure_node(g, term_node, term_type, term_name)
            g.add((term_node, RDF.type, BIOLINK.DiseaseOrPhenotypicFeature))

            assoc_node = mint_assoc_node("ot_target_condition", approved_symbol, term_id or term_name)
            g.add((assoc_node, RDF.type, BIOLINK.GeneToDiseaseOrPhenotypicFeatureAssociation))
            g.add((assoc_node, BIOLINK.subject, gene_node))
            g.add((assoc_node, BIOLINK.object, term_node))
            g.add((assoc_node, BIOLINK.predicate, BIOLINK.gene_associated_with_condition))
            g.add((assoc_node, PROV.wasGeneratedBy, run_node))
            g.add((assoc_node, PROV.wasDerivedFrom, source_node))

            for ds in row.get("datasourceScores", []) or []:
                ds_id = ds.get("id")
                ds_score = ds.get("score")
                if ds_id is None or ds_score is None:
                    continue

                if ds_id == "genetics":
                    safe_literal(g, assoc_node, MTR.ot_genetics_score, float(ds_score), XSD.double)

                score_node = mint_assoc_node("ot_datasource_score", approved_symbol, term_id or term_name, ds_id)
                g.add((score_node, RDF.type, MTR.OpenTargetsDatasourceScore))
                g.add((score_node, MTR.about_association, assoc_node))
                safe_literal(g, score_node, MTR.datasource_id, ds_id, XSD.string)
                safe_literal(g, score_node, MTR.datasource_score, float(ds_score), XSD.double)

        time.sleep(sleep_s)

    logger.info("--- Open Targets Integration Complete ---")
    return g