import logging
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from rdflib import Literal, RDF, URIRef
from rdflib.namespace import RDFS, XSD
from urllib3.util.retry import Retry

from schema_definition import MTR, BIOLINK, PROV, HGNC

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(OUTPUT_DIR / "string_integration.log", mode="w"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Keep configurable. For reproducible production runs, use a versioned STRING base.
STRING_API_BASE = "https://string-db.org/api"
STRING_SPECIES = 9606
CALLER_IDENTITY = "unil_rqtl_knowledge_graph_project"


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
    session.headers.update({"User-Agent": CALLER_IDENTITY})
    return session


SESSION = build_http_session()


def safe_id(text: str) -> str:
    return urllib.parse.quote(str(text).strip().replace(" ", "_").replace('"', ""))


def ensure_node(g, node, rdf_type=None, label=None):
    if rdf_type is not None:
        g.add((node, RDF.type, rdf_type))
    if label:
        g.set((node, RDFS.label, Literal(str(label), datatype=XSD.string)))


def get_gene_symbol(g, gene_node) -> str:
    label = next(g.objects(gene_node, RDFS.label), None)
    if label:
        return str(label).strip()
    return urllib.parse.unquote(str(gene_node).rstrip("/").split("/")[-1].strip())


def post_string_api(method: str, params: dict, timeout: int = 30):
    url = f"{STRING_API_BASE}/json/{method}"
    try:
        response = SESSION.post(url, data=params, timeout=timeout)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as exc:
        logger.warning("STRING API request failed for %s: %s", method, exc)
        return None


def resolve_string_ids(symbols, species=STRING_SPECIES):
    """
    Resolve gene symbols to STRING internal identifiers.
    Returns a list of mapping rows.
    """
    params = {
        "identifiers": "\r".join(symbols),
        "species": species,
        "echo_query": 1,
        "caller_identity": CALLER_IDENTITY,
    }
    return post_string_api("get_string_ids", params) or []


def fetch_string_network(string_ids, species=STRING_SPECIES, min_score=0.4, network_type="physical"):
    """
    Fetch STRING network among the resolved STRING IDs.
    network_type:
      - physical   => closer to protein-protein interaction semantics
      - functional => broader functional association network
    """
    required_score = max(0, min(1000, int(float(min_score) * 1000)))
    params = {
        "identifiers": "\r".join(string_ids),
        "species": species,
        "required_score": required_score,
        "network_type": network_type,
        "caller_identity": CALLER_IDENTITY,
    }
    return post_string_api("network", params) or []


def mint_run_node():
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return MTR[f"activity/string_enrichment/{ts}"]


def mint_interaction_node(id_a: str, id_b: str, network_type: str):
    pair = sorted([id_a, id_b])
    return MTR[f"string_interaction/{safe_id(pair[0])}__{safe_id(pair[1])}__{safe_id(network_type)}"]


def enrich_graph_with_string(g, min_score=0.4, network_type="physical", add_direct_edges=False):
    """
    Enrich the graph with STRING interactions among genes already present.

    Args:
        g: rdflib Graph
        min_score: 0.0-1.0 STRING confidence threshold
        network_type: 'physical' or 'functional'
        add_direct_edges: if True, also add direct biolink:interacts_with edges
    """
    logger.info("--- Starting STRING integration ---")

    gene_nodes = list(g.subjects(RDF.type, MTR.Gene))
    if not gene_nodes:
        logger.warning("No genes found in the graph. Skipping STRING integration.")
        return g

    # Map existing local gene symbols to nodes
    symbol_to_gene_node = {}
    for gene_node in gene_nodes:
        symbol = get_gene_symbol(g, gene_node)
        if symbol and symbol not in symbol_to_gene_node:
            symbol_to_gene_node[symbol] = gene_node
            ensure_node(g, gene_node, MTR.Gene, symbol)

    gene_symbols = sorted(symbol_to_gene_node.keys())
    logger.info("Found %d genes. Resolving identifiers in STRING...", len(gene_symbols))

    resolved_rows = resolve_string_ids(gene_symbols)
    if not resolved_rows:
        logger.warning("STRING identifier resolution returned no mappings.")
        return g

    # Build symbol -> STRING ID and STRING ID -> local node maps
    stringid_to_gene_node = {}
    resolution_report = []

    for row in resolved_rows:
        query_item = row.get("queryItem")
        string_id = row.get("stringId")
        preferred_name = row.get("preferredName")

        if not query_item or not string_id:
            continue

        gene_node = symbol_to_gene_node.get(query_item)
        if gene_node is None:
            # fallback on preferred name if available
            gene_node = symbol_to_gene_node.get(preferred_name)

        if gene_node is None:
            continue

        stringid_to_gene_node[string_id] = gene_node
        ensure_node(g, gene_node, MTR.Gene, preferred_name or query_item)
        g.add((gene_node, MTR.string_identifier, Literal(string_id, datatype=XSD.string)))

        resolution_report.append(
            {
                "Input_Symbol": query_item,
                "Preferred_Name": preferred_name,
                "STRING_ID": string_id,
                "Status": "Mapped",
            }
        )

    if not stringid_to_gene_node:
        logger.warning("No local genes could be mapped to STRING IDs.")
        return g

    logger.info("Resolved %d genes. Fetching STRING %s network...", len(stringid_to_gene_node), network_type)

    interactions = fetch_string_network(
        list(stringid_to_gene_node.keys()),
        min_score=min_score,
        network_type=network_type,
    )

    if not interactions:
        logger.info("No STRING interactions returned.")
        pd.DataFrame(resolution_report).to_csv(OUTPUT_DIR / "string_resolution_report.csv", index=False)
        return g

    logger.info("Retrieved %d interactions from STRING.", len(interactions))

    run_node = mint_run_node()
    source_node = URIRef("https://string-db.org/")
    g.add((run_node, RDF.type, PROV.Activity))
    g.add((run_node, RDFS.label, Literal("STRING enrichment run", datatype=XSD.string)))
    g.add((source_node, RDF.type, PROV.Entity))
    g.add((source_node, RDFS.label, Literal("STRING database", datatype=XSD.string)))

    interaction_report = []
    seen_pairs = set()

    for row in interactions:
        string_id_a = row.get("stringId_A")
        string_id_b = row.get("stringId_B")
        preferred_a = row.get("preferredName_A")
        preferred_b = row.get("preferredName_B")
        score = row.get("score")

        if not string_id_a or not string_id_b:
            continue

        gene_a_node = stringid_to_gene_node.get(string_id_a)
        gene_b_node = stringid_to_gene_node.get(string_id_b)

        # network endpoint among multiple inputs should return only query proteins,
        # but keep a fallback
        if gene_a_node is None:
            gene_a_node = HGNC[safe_id(preferred_a or string_id_a)]
            ensure_node(g, gene_a_node, MTR.Gene, preferred_a or string_id_a)

        if gene_b_node is None:
            gene_b_node = HGNC[safe_id(preferred_b or string_id_b)]
            ensure_node(g, gene_b_node, MTR.Gene, preferred_b or string_id_b)

        pair_key = tuple(sorted([string_id_a, string_id_b]))
        if pair_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)

        interaction_node = mint_interaction_node(string_id_a, string_id_b, network_type)

        g.add((interaction_node, RDF.type, BIOLINK.PairwiseGeneToGeneInteraction))
        g.add((interaction_node, BIOLINK.subject, gene_a_node))
        g.add((interaction_node, BIOLINK.object, gene_b_node))
        g.add((interaction_node, BIOLINK.predicate, BIOLINK.interacts_with))

        if score is not None:
            g.add((interaction_node, MTR.string_combined_score, Literal(float(score), datatype=XSD.double)))

        # Optional evidence channels from STRING response
        for field, predicate in [
            ("escore", MTR.string_experimental_score),
            ("dscore", MTR.string_database_score),
            ("tscore", MTR.string_textmining_score),
            ("ascore", MTR.string_coexpression_score),
            ("nscore", MTR.string_neighborhood_score),
            ("fscore", MTR.string_fusion_score),
            ("pscore", MTR.string_phylogenetic_score),
        ]:
            value = row.get(field)
            if value is not None:
                g.add((interaction_node, predicate, Literal(float(value), datatype=XSD.double)))

        g.add((interaction_node, PROV.wasGeneratedBy, run_node))
        g.add((interaction_node, PROV.wasDerivedFrom, source_node))
        g.add((interaction_node, MTR.string_network_type, Literal(network_type, datatype=XSD.string)))

        if add_direct_edges:
            g.add((gene_a_node, BIOLINK.interacts_with, gene_b_node))
            g.add((gene_b_node, BIOLINK.interacts_with, gene_a_node))

        interaction_report.append(
            {
                "Gene_A": preferred_a or string_id_a,
                "Gene_B": preferred_b or string_id_b,
                "STRING_ID_A": string_id_a,
                "STRING_ID_B": string_id_b,
                "Score": score,
                "Network_Type": network_type,
            }
        )

    pd.DataFrame(resolution_report).to_csv(OUTPUT_DIR / "string_resolution_report.csv", index=False)
    pd.DataFrame(interaction_report).to_csv(OUTPUT_DIR / "string_interaction_report.csv", index=False)

    logger.info("--- STRING integration complete ---")
    return g