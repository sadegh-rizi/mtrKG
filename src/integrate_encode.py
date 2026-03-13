import logging
import time
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from rdflib import Graph, Literal, RDF, URIRef
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
        logging.FileHandler(OUTPUT_DIR / "ensembl_regulatory_integration.log", mode="w"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

ENSEMBL_SERVER = "https://rest.ensembl.org"


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
            "Content-Type": "application/json",
            "User-Agent": "rqtl-kg-integrator/1.0"
        }
    )
    return session


SESSION = build_http_session()


def safe_id(text: str) -> str:
    return urllib.parse.quote(str(text).strip().replace(" ", "_").replace('"', ""))


def ensure_node(g: Graph, node, rdf_type=None, label=None):
    if rdf_type is not None:
        g.add((node, RDF.type, rdf_type))
    if label:
        g.set((node, RDFS.label, Literal(str(label), datatype=XSD.string)))


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


def rsid_from_node(snp_node) -> str:
    return str(snp_node).rstrip("/").split("/")[-1].strip()


def mint_run_node():
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return MTR[f"activity/ensembl_regulatory_enrichment/{ts}"]


def classify_regulatory_feature(feature: dict):
    desc = (feature.get("description") or feature.get("feature_type") or "").lower()

    if "ctcf" in desc or "binding site" in desc or "tf binding" in desc:
        return BIOLINK.TranscriptionFactorBindingSite

    if "open chromatin" in desc or "accessible" in desc:
        return BIOLINK.AccessibleDnaRegion

    return BIOLINK.RegulatoryRegion


def classify_motif_feature(feature: dict):
    return BIOLINK.TranscriptionFactorBindingSite


def get_snp_location(rsid: str):
    """
    Use Ensembl variation endpoint to resolve rsID to genomic coordinates.
    """
    ext = f"/variation/human/{rsid}"
    try:
        r = SESSION.get(ENSEMBL_SERVER + ext, timeout=30)
        if r.status_code == 404:
            return None, None

        r.raise_for_status()
        data = r.json()

        mappings = data.get("mappings", [])
        if not mappings:
            return None, None

        # keep first usable mapping
        for m in mappings:
            seq_region = m.get("seq_region_name")
            start = m.get("start")
            end = m.get("end")
            if seq_region and start and end:
                return f"{seq_region}:{start}-{end}", m.get("assembly_name")

        return None, None

    except requests.exceptions.RequestException as exc:
        logger.warning("Variation lookup failed for %s: %s", rsid, exc)
        return None, None


def get_overlap_by_region(location: str, feature_type: str):
    """
    Ensembl overlap/region endpoint for regulatory or motif features.
    """
    ext = f"/overlap/region/human/{location}"
    params = {"feature": feature_type}
    try:
        r = SESSION.get(ENSEMBL_SERVER + ext, params=params, timeout=30)
        if r.status_code == 404:
            return []

        r.raise_for_status()
        data = r.json()
        return data if isinstance(data, list) else []

    except requests.exceptions.RequestException as exc:
        logger.warning("Overlap lookup failed for %s (%s): %s", location, feature_type, exc)
        return []


def add_feature_overlap(
    g: Graph,
    snp_node,
    run_node,
    source_node,
    rsid: str,
    location: str,
    assembly_name: str | None,
    feature: dict,
    feature_kind: str,
):
    feature_id = feature.get("id")
    if not feature_id:
        return

    feature_class = classify_motif_feature(feature) if feature_kind == "motif" else classify_regulatory_feature(feature)

    feature_desc = (
        feature.get("description")
        or feature.get("feature_type")
        or feature.get("logic_name")
        or ("Motif feature" if feature_kind == "motif" else "Regulatory feature")
    )

    reg_node = MTR[f"{feature_kind}_feature/{safe_id(feature_id)}"]
    ensure_node(g, reg_node, feature_class, f"{feature_desc} ({feature_id})")

    safe_literal(g, reg_node, MTR.external_feature_id, feature_id, XSD.string)
    safe_literal(g, reg_node, MTR.feature_description, feature_desc, XSD.string)
    safe_literal(g, reg_node, MTR.chromosome, feature.get("seq_region_name"), XSD.string)
    safe_literal(g, reg_node, MTR.start_position, feature.get("start"), XSD.integer)
    safe_literal(g, reg_node, MTR.end_position, feature.get("end"), XSD.integer)
    safe_literal(g, reg_node, MTR.assembly_name, assembly_name, XSD.string)

    overlap_assoc = MTR[f"variant_{feature_kind}_overlap/{safe_id(rsid)}__{safe_id(feature_id)}"]
    g.add((overlap_assoc, RDF.type, MTR.VariantToRegulatoryFeatureOverlapAssociation))
    g.add((overlap_assoc, BIOLINK.subject, snp_node))
    g.add((overlap_assoc, BIOLINK.object, reg_node))
    g.add((overlap_assoc, BIOLINK.predicate, BIOLINK.overlaps))
    g.add((overlap_assoc, PROV.wasGeneratedBy, run_node))
    g.add((overlap_assoc, PROV.wasDerivedFrom, source_node))

    safe_literal(g, overlap_assoc, MTR.assembly_name, assembly_name, XSD.string)
    safe_literal(g, overlap_assoc, MTR.query_region, location, XSD.string)
    safe_literal(g, overlap_assoc, MTR.feature_category, feature_kind, XSD.string)

    # optional denormalized edge
    g.add((snp_node, BIOLINK.overlaps, reg_node))


def enrich_graph_with_ensembl_regulatory(g: Graph, max_snps=None, sleep_s: float = 0.15) -> Graph:
    logger.info("--- Starting Ensembl Regulatory Build integration ---")

    snp_nodes = list(g.subjects(RDF.type, MTR.SNP))
    if max_snps is not None:
        snp_nodes = snp_nodes[:max_snps]

    if not snp_nodes:
        logger.warning("No SNPs found in the graph.")
        return g

    run_node = mint_run_node()
    source_node = URIRef("https://www.ensembl.org/info/genome/funcgen/index.html")
    ensure_node(g, run_node, PROV.Activity, "Ensembl regulatory enrichment run")
    ensure_node(g, source_node, PROV.Entity, "Ensembl Regulatory Build")

    integration_report = []

    for i, snp_node in enumerate(snp_nodes):
        rsid = rsid_from_node(snp_node)
        logger.info("[%d/%d] Checking regulatory overlaps for %s...", i + 1, len(snp_nodes), rsid)

        try:
            location, assembly_name = get_snp_location(rsid)
            if not location:
                logger.info("   -> Could not find genomic coordinates for %s.", rsid)
                integration_report.append(
                    {"SNP": rsid, "Status": "Unmapped_Location", "Features_Found": 0}
                )
                time.sleep(sleep_s)
                continue

            regulatory_features = get_overlap_by_region(location, "regulatory")
            motif_features = get_overlap_by_region(location, "motif")

            total_features = len(regulatory_features) + len(motif_features)
            if total_features == 0:
                logger.info("   -> No overlapping regulatory or motif features found.")
                integration_report.append(
                    {"SNP": rsid, "Status": "No_Overlap", "Features_Found": 0}
                )
                time.sleep(sleep_s)
                continue

            logger.info(
                "   -> Found %d features (%d regulatory, %d motif).",
                total_features,
                len(regulatory_features),
                len(motif_features),
            )
            integration_report.append(
                {"SNP": rsid, "Status": "Mapped", "Features_Found": total_features}
            )

            for feature in regulatory_features:
                add_feature_overlap(
                    g=g,
                    snp_node=snp_node,
                    run_node=run_node,
                    source_node=source_node,
                    rsid=rsid,
                    location=location,
                    assembly_name=assembly_name,
                    feature=feature,
                    feature_kind="regulatory",
                )

            for feature in motif_features:
                add_feature_overlap(
                    g=g,
                    snp_node=snp_node,
                    run_node=run_node,
                    source_node=source_node,
                    rsid=rsid,
                    location=location,
                    assembly_name=assembly_name,
                    feature=feature,
                    feature_kind="motif",
                )

        except Exception as exc:
            logger.exception("   -> ERROR processing %s: %s", rsid, exc)
            integration_report.append(
                {"SNP": rsid, "Status": "Failed", "Features_Found": 0}
            )

        time.sleep(sleep_s)

    report_path = OUTPUT_DIR / "ensembl_regulatory_mapping_report.csv"
    pd.DataFrame(integration_report).to_csv(report_path, index=False)

    logger.info("Saved mapping report to %s", report_path)
    logger.info("--- Ensembl Regulatory Build integration complete ---")
    return g