import logging
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
from rdflib import Graph, Literal, RDF, URIRef
from rdflib.namespace import RDFS, SKOS, XSD

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
        logging.FileHandler(OUTPUT_DIR / "hmdb_bulk_integration.log", mode="w"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def normalize_label(text: str) -> str:
    return " ".join(str(text).strip().split())


def safe_id(text: str) -> str:
    return urllib.parse.quote(normalize_label(text).replace(" ", "_"))


def ensure_node(g: Graph, node, rdf_type=None, label=None):
    if rdf_type is not None:
        g.add((node, RDF.type, rdf_type))
    if label:
        g.set((node, RDFS.label, Literal(normalize_label(label), datatype=XSD.string)))


def add_location_association(g: Graph, met_node, loc_node, source_node, assoc_prefix: str):
    """
    Association-centric modeling with provenance.
    """
    assoc_node = MTR[f"{assoc_prefix}/{safe_id(str(met_node))}__{safe_id(str(loc_node))}"]
    g.add((assoc_node, RDF.type, MTR.MetaboliteLocationAssociation))
    g.add((assoc_node, BIOLINK.subject, met_node))
    g.add((assoc_node, BIOLINK.object, loc_node))
    g.add((assoc_node, BIOLINK.predicate, BIOLINK.located_in))
    g.add((assoc_node, PROV.wasDerivedFrom, source_node))


def extract_first_text(elem, child_name: str):
    for child in elem:
        if child.tag.split("}")[-1] == child_name and child.text:
            return normalize_label(child.text)
    return None


def enrich_graph_with_hmdb(
    g: Graph,
    xml_file_path="../data/hmdb_metabolites.xml",
    use_association_nodes: bool = True,
):
    logger.info("--- Starting Local HMDB Bulk XML Integration ---")

    xml_path = Path(xml_file_path)
    hmdb_source_node = URIRef("https://hmdb.ca/")
    g.add((hmdb_source_node, RDF.type, PROV.Entity))
    g.add((hmdb_source_node, RDFS.label, Literal("HMDB bulk metabolite XML", datatype=XSD.string)))

    # 1. Find HMDB-linked metabolites already in the graph
    target_hmdb_nodes = {}
    for met_node, _, hmdb_node in g.triples((None, SKOS.exactMatch, None)):
        hmdb_uri = str(hmdb_node)
        if hmdb_uri.startswith(str(HMDB)):
            hmdb_id = hmdb_uri.rstrip("/").split("/")[-1]
            target_hmdb_nodes[hmdb_id] = met_node

    if not target_hmdb_nodes:
        logger.warning("No HMDB cross-references found in the graph. Skipping.")
        return g

    logger.info("Looking for %d specific metabolites in the bulk file...", len(target_hmdb_nodes))

    integration_report = []
    found_ids = set()

    # caches to avoid duplicate node creation
    cellular_cache = {}
    tissue_cache = {}
    biospecimen_cache = {}

    try:
        context = ET.iterparse(xml_path, events=("end",))

        for _, elem in context:
            tag_name = elem.tag.split("}")[-1]

            if tag_name != "metabolite":
                continue

            accession = extract_first_text(elem, "accession")
            if not accession:
                elem.clear()
                continue

            if accession not in target_hmdb_nodes:
                elem.clear()
                continue

            met_node = target_hmdb_nodes[accession]
            found_ids.add(accession)
            logger.info("   -> Found match in XML: %s", accession)

            cellular = set()
            biospecimen = set()
            tissue = set()

            # HMDB location fields often appear under biological properties sections.
            for desc in elem.iter():
                child_tag = desc.tag.split("}")[-1]
                if not desc.text:
                    continue

                value = normalize_label(desc.text)
                if not value or value.lower() == "not available":
                    continue

                if child_tag == "cellular_location":
                    cellular.add(value)
                elif child_tag == "biospecimen":
                    biospecimen.add(value)
                elif child_tag == "tissue":
                    tissue.add(value)

            total_locs = len(cellular) + len(biospecimen) + len(tissue)

            if total_locs == 0:
                integration_report.append(
                    {"Metabolite": accession, "Status": "Unmapped_NoData", "Locations_Found": 0}
                )
                logger.info("      No location data available for %s", accession)
                elem.clear()
                if len(found_ids) >= len(target_hmdb_nodes):
                    break
                continue

            # 2. Add cellular locations
            for cell_loc in sorted(cellular):
                if cell_loc not in cellular_cache:
                    loc_node = MTR[f"cellular_component/{safe_id(cell_loc)}"]
                    ensure_node(g, loc_node, BIOLINK.CellularComponent, cell_loc)
                    cellular_cache[cell_loc] = loc_node
                loc_node = cellular_cache[cell_loc]

                if use_association_nodes:
                    add_location_association(g, met_node, loc_node, hmdb_source_node, "hmdb_metabolite_location")
                else:
                    g.add((met_node, BIOLINK.located_in, loc_node))

            # 3. Add tissue locations
            for tissue_loc in sorted(tissue):
                if tissue_loc not in tissue_cache:
                    loc_node = MTR[f"anatomical_entity/{safe_id(tissue_loc)}"]
                    ensure_node(g, loc_node, BIOLINK.GrossAnatomicalStructure, tissue_loc)
                    tissue_cache[tissue_loc] = loc_node
                loc_node = tissue_cache[tissue_loc]

                if use_association_nodes:
                    add_location_association(g, met_node, loc_node, hmdb_source_node, "hmdb_metabolite_location")
                else:
                    g.add((met_node, BIOLINK.located_in, loc_node))

            # 4. Add biospecimen locations
            for biospec in sorted(biospecimen):
                if biospec not in biospecimen_cache:
                    loc_node = MTR[f"biospecimen/{safe_id(biospec)}"]
                    ensure_node(g, loc_node, BIOLINK.MaterialSample, biospec)
                    biospecimen_cache[biospec] = loc_node
                loc_node = biospecimen_cache[biospec]

                if use_association_nodes:
                    add_location_association(g, met_node, loc_node, hmdb_source_node, "hmdb_metabolite_location")
                else:
                    g.add((met_node, BIOLINK.located_in, loc_node))

            integration_report.append(
                {"Metabolite": accession, "Status": "Mapped", "Locations_Found": total_locs}
            )

            elem.clear()

            if len(found_ids) >= len(target_hmdb_nodes):
                logger.info("All %d HMDB targets found. Halting XML scan early.", len(target_hmdb_nodes))
                break

    except ET.ParseError as exc:
        logger.error("Failed to parse XML file: %s", exc)
    except FileNotFoundError:
        logger.error("Could not find the file at %s", xml_path)

    for target_id in target_hmdb_nodes.keys():
        if target_id not in found_ids:
            integration_report.append(
                {"Metabolite": target_id, "Status": "Not_Found_in_XML", "Locations_Found": 0}
            )

    df_report = pd.DataFrame(integration_report)
    report_path = OUTPUT_DIR / "hmdb_bulk_mapping_report.csv"
    df_report.to_csv(report_path, index=False)

    logger.info("Saved HMDB mapping report to %s", report_path)
    logger.info("--- HMDB Bulk Integration Complete ---")
    return g