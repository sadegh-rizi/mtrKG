import logging
import math
import re
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
        logging.FileHandler(OUTPUT_DIR / "rhea_integration.log", mode="w"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

RHEA = Namespace("http://rdf.rhea-db.org/")
RHEA_SPARQL_ENDPOINT = "https://sparql.rhea-db.org/sparql"

CHEBI_URI_RE = re.compile(r"CHEBI_(\d+)$")


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
            "Accept": "application/sparql-results+json",
            "User-Agent": "Metabolite-Ratio-KG-Builder/1.0",
        }
    )
    return session


SESSION = build_http_session()


def chunked(items, chunk_size):
    for i in range(0, len(items), chunk_size):
        yield items[i:i + chunk_size]


def safe_literal(g: Graph, subject, predicate, value, datatype=None):
    if value in [None, "", "NA", "NaN"]:
        return
    try:
        if datatype is not None:
            g.add((subject, predicate, Literal(value, datatype=datatype)))
        else:
            g.add((subject, predicate, Literal(value)))
    except Exception:
        logger.debug("Skipped invalid literal for %s %s %r", subject, predicate, value)


def ensure_node(g: Graph, node, rdf_type=None, label=None):
    if rdf_type is not None:
        g.add((node, RDF.type, rdf_type))
    if label:
        g.set((node, RDFS.label, Literal(str(label), datatype=XSD.string)))


def extract_chebi_number(node) -> str | None:
    match = CHEBI_URI_RE.search(str(node))
    if match:
        return match.group(1)
    return None


def mint_run_node():
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return MTR[f"activity/rhea_enrichment/{ts}"]


def mint_assoc_node(chebi_num: str, reaction_uri: str, met_node) -> URIRef:
    safe = urllib.parse.quote(
        f"{chebi_num}__{reaction_uri}__{str(met_node)}".replace(" ", "_").replace('"', "")
    )
    return MTR[f"rhea_reaction_participant/{safe}"]


def build_rhea_batch_query(chebi_nums):
    values = " ".join(f"CHEBI:{num}" for num in chebi_nums)
    return f"""
    PREFIX rh: <http://rdf.rhea-db.org/>
    PREFIX CHEBI: <http://purl.obolibrary.org/obo/CHEBI_>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

    SELECT DISTINCT ?matchedChebi ?reaction ?equation ?ecNumber
    WHERE {{
      VALUES ?matchedChebi {{ {values} }}

      ?reaction rdfs:subClassOf rh:Reaction .
      ?reaction rh:status rh:Approved .
      ?reaction rh:equation ?equation .
      ?reaction rh:side ?reactionSide .
      ?reactionSide rh:contains ?participant .
      ?participant rh:compound ?compound .

      ?compound (rh:chebi|(rh:reactivePart/rh:chebi)|rh:underlyingChebi) ?matchedChebi .

      OPTIONAL {{ ?reaction rh:ec ?ecNumber . }}
    }}
    """


def fetch_rhea_reactions_batch(chebi_nums):
    query = build_rhea_batch_query(chebi_nums)
    try:
        response = SESSION.get(
            RHEA_SPARQL_ENDPOINT,
            params={"query": query},
            timeout=45,
        )
        response.raise_for_status()
        data = response.json()
        return data.get("results", {}).get("bindings", [])
    except requests.exceptions.RequestException as exc:
        logger.warning("Rhea SPARQL batch failed for %s: %s", chebi_nums, exc)
        return []


def enrich_graph_with_rhea(g: Graph, chunk_size: int = 50) -> Graph:
    logger.info("--- Starting batched Rhea integration ---")

    metabolite_nodes = list(g.subjects(RDF.type, MTR.Metabolite))

    # Deduplicate by ChEBI
    unique_chebis = {}
    for met_node in metabolite_nodes:
        chebi_num = extract_chebi_number(met_node)
        if not chebi_num:
            continue
        unique_chebis.setdefault(chebi_num, []).append(met_node)

    logger.info("Querying Rhea for %d unique ChEBI metabolites...", len(unique_chebis))

    integration_report = []
    run_node = mint_run_node()
    source_node = URIRef("https://www.rhea-db.org/")

    ensure_node(g, run_node, PROV.Activity, "Rhea enrichment run")
    ensure_node(g, source_node, PROV.Entity, "Rhea database")

    results_by_chebi = {chebi: [] for chebi in unique_chebis.keys()}

    chebi_list = list(unique_chebis.keys())
    total_chunks = math.ceil(len(chebi_list) / chunk_size) if chebi_list else 0

    for idx, chebi_chunk in enumerate(chunked(chebi_list, chunk_size), start=1):
        logger.info("[%d/%d] Querying Rhea for %d ChEBI IDs...", idx, total_chunks, len(chebi_chunk))
        rows = fetch_rhea_reactions_batch(chebi_chunk)

        for row in rows:
            matched_chebi_uri = row.get("matchedChebi", {}).get("value")
            reaction_uri = row.get("reaction", {}).get("value")
            equation = row.get("equation", {}).get("value")
            ec_number = row.get("ecNumber", {}).get("value")

            if not matched_chebi_uri or not reaction_uri:
                continue

            chebi_num = matched_chebi_uri.rstrip("/").split("CHEBI_")[-1]
            results_by_chebi.setdefault(chebi_num, []).append(
                {
                    "reaction_uri": reaction_uri,
                    "equation": equation,
                    "ec_number": ec_number,
                }
            )

    logger.info("Writing Rhea results to the Knowledge Graph...")

    for chebi_num, reaction_rows in results_by_chebi.items():
        integration_report.append(
            {
                "Metabolite": f"CHEBI:{chebi_num}",
                "Status": "Mapped" if reaction_rows else "Unmapped",
                "Reactions_Found": len(reaction_rows),
            }
        )

        if not reaction_rows:
            continue

        # deduplicate repeated rows per ChEBI
        seen = set()
        for rxn in reaction_rows:
            key = (rxn["reaction_uri"], rxn.get("equation"), rxn.get("ec_number"))
            if key in seen:
                continue
            seen.add(key)

            reaction_node = URIRef(rxn["reaction_uri"])
            ensure_node(g, reaction_node, MTR.BiochemicalReaction, rxn.get("equation") or rxn["reaction_uri"])
            g.add((reaction_node, PROV.wasDerivedFrom, source_node))

            if rxn.get("equation"):
                safe_literal(g, reaction_node, MTR.rhea_equation, rxn["equation"], XSD.string)

            if rxn.get("ec_number"):
                safe_literal(g, reaction_node, MTR.rhea_ec_number, rxn["ec_number"], XSD.string)

            safe_literal(g, reaction_node, MTR.rhea_status, "Approved", XSD.string)

            # connect every local metabolite node using this ChEBI
            for met_node in unique_chebis[chebi_num]:
                assoc_node = mint_assoc_node(chebi_num, rxn["reaction_uri"], met_node)

                g.add((assoc_node, RDF.type, BIOLINK.ReactionToParticipantAssociation))
                g.add((assoc_node, BIOLINK.subject, reaction_node))
                g.add((assoc_node, BIOLINK.object, met_node))
                g.add((assoc_node, BIOLINK.predicate, BIOLINK.has_participant))
                g.add((assoc_node, PROV.wasGeneratedBy, run_node))
                g.add((assoc_node, PROV.wasDerivedFrom, source_node))

    report_path = OUTPUT_DIR / "rhea_mapping_report.csv"
    pd.DataFrame(integration_report).to_csv(report_path, index=False)
    logger.info("Saved mapping report to %s", report_path)
    logger.info("--- Rhea integration complete ---")

    return g