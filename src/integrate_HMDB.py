import logging
import urllib.parse
import pandas as pd
import xml.etree.ElementTree as ET
from rdflib import Literal, RDF
from rdflib.namespace import RDFS, SKOS
from src.schema_definition import *  # Ensure BIOLINK, MTR, PROV are defined

# ==========================================
# 1. SETUP LOGGING
# ==========================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("../output/hmdb_bulk_integration.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def enrich_graph_with_hmdb(g, xml_file_path="../data/hmdb_metabolites.xml"):
    logger.info("--- Starting Local HMDB Bulk XML Integration ---")

    # 1. Find the target metabolites in your graph first
    # This prevents us from processing all 250,000+ HMDB metabolites!
    target_hmdb_nodes = {}
    for met_node, _, hmdb_node in g.triples((None, SKOS.exactMatch, None)):
        hmdb_uri = str(hmdb_node)
        if "HMDB" in hmdb_uri:
            hmdb_id = hmdb_uri.split('/')[-1]  # Extracts "HMDB0000168"
            target_hmdb_nodes[hmdb_id] = met_node

    if not target_hmdb_nodes:
        logger.warning("No HMDB cross-references found in the graph. Skipping.")
        return g

    logger.info(f"Looking for {len(target_hmdb_nodes)} specific metabolites in the bulk file...")

    integration_report = []
    found_count = 0

    # 2. Use iterparse to stream the XML file without crashing RAM
    try:
        context = ET.iterparse(xml_file_path, events=('end',))

        for event, elem in context:
            # HMDB XML tags include namespaces (e.g., {http://www.hmdb.ca}metabolite)
            # We strip the namespace to make it easy to read
            tag_name = elem.tag.split('}')[-1]

            if tag_name == 'metabolite':
                # Grab the main accession ID for this metabolite block
                accession_elem = None
                for child in elem:
                    if child.tag.split('}')[-1] == 'accession':
                        accession_elem = child
                        break

                if accession_elem is not None and accession_elem.text in target_hmdb_nodes:
                    hmdb_id = accession_elem.text
                    met_node = target_hmdb_nodes[hmdb_id]

                    logger.info(f"   -> Found match in XML: {hmdb_id}")

                    # Extract the nested location data
                    cellular = []
                    biospecimen = []
                    tissue = []

                    # Deep search within this specific metabolite block
                    for desc in elem.iter():
                        child_tag = desc.tag.split('}')[-1]
                        if child_tag == 'cellular_location' and desc.text:
                            cellular.append(desc.text)
                        elif child_tag == 'biospecimen' and desc.text:
                            biospecimen.append(desc.text)
                        elif child_tag == 'tissue' and desc.text:
                            tissue.append(desc.text)

                    total_locs = len(cellular) + len(biospecimen) + len(tissue)

                    if total_locs > 0:
                        integration_report.append(
                            {"Metabolite": hmdb_id, "Status": "Mapped", "Locations_Found": total_locs})

                        # 3. Add to Graph
                        for cell_loc in cellular:
                            safe_id = urllib.parse.quote(cell_loc.strip().replace(" ", "_"))
                            loc_node = MTR[f"CellularComponent_{safe_id}"]
                            g.add((loc_node, RDF.type, BIOLINK.CellularComponent))
                            g.add((loc_node, RDFS.label, Literal(cell_loc)))
                            g.add((met_node, BIOLINK.located_in, loc_node))

                        for macro_loc in tissue + biospecimen:
                            safe_id = urllib.parse.quote(macro_loc.strip().replace(" ", "_"))
                            loc_node = MTR[f"Anatomy_{safe_id}"]
                            g.add((loc_node, RDF.type, BIOLINK.GrossAnatomicalStructure))
                            g.add((loc_node, RDFS.label, Literal(macro_loc)))
                            g.add((met_node, BIOLINK.located_in, loc_node))
                    else:
                        integration_report.append(
                            {"Metabolite": hmdb_id, "Status": "Unmapped_NoData", "Locations_Found": 0})
                        logger.info(f"      No location data available for {hmdb_id}")

                    found_count += 1

                # 4. CRITICAL MEMORY MANAGEMENT: Clear the parsed element from RAM!
                elem.clear()

                # 5. Optimization: Stop parsing if we've found all our targets
                if found_count >= len(target_hmdb_nodes):
                    logger.info(f"All {len(target_hmdb_nodes)} targets found! Halting XML scan early to save time.")
                    break

    except ET.ParseError as e:
        logger.error(f"Failed to parse XML file: {e}")
    except FileNotFoundError:
        logger.error(f"Could not find the file at {xml_file_path}. Please check the path.")

    logger.info("--- HMDB Bulk Integration Complete ---")

    # Identify any metabolites we failed to find entirely
    mapped_ids = [row["Metabolite"] for row in integration_report]
    for target_id in target_hmdb_nodes.keys():
        if target_id not in mapped_ids:
            integration_report.append({"Metabolite": target_id, "Status": "Not_Found_in_XML", "Locations_Found": 0})

    df_report = pd.DataFrame(integration_report)
    df_report.to_csv("../output/hmdb_bulk_mapping_report.csv", index=False)

    return g