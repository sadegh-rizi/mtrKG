import requests
import time
import urllib.parse
import logging
from rdflib import Graph, Literal, RDF, URIRef
from rdflib.namespace import RDFS, XSD
# Assuming your schema defines MTR, BIOLINK, EFO, CHEMBL, PROV, OBAN
from schema_definition import *

# ==========================================
# 1. SETUP LOGGING
# ==========================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("../output/opentargets_integration.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

OPENTARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"


def fetch_opentargets_data(gene_symbol):
    """
    Executes a GraphQL query against Open Targets using the updated schema to find
    tractability, safety liabilities, known drugs, and associated diseases.
    """
    query = """
    query searchDrugTargets($symbol: String!) {
      search(queryString: $symbol, entityNames: ["target"], page: {index: 0, size: 1}) {
        hits {
          id
          object {
            ... on Target {
              id
              approvedSymbol

              # 1. Updated Tractability Schema
              tractability {
                modality
                label
                value
              }

              # 2. Updated Safety Liabilities Schema
              safetyLiabilities {
                event
                eventId
                biosamples {
                  tissueLabel
                }
              }

              # 3. Existing Known Drugs
              knownDrugs(size: 10) {
                rows {
                  phase
                  drug {
                    id
                    name
                  }
                }
              }

              # 4. NEW: Associated Diseases
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
        }
      }
    }
    """

    variables = {"symbol": gene_symbol}

    try:
        response = requests.post(OPENTARGETS_API, json={"query": query, "variables": variables})
        if response.status_code == 200:
            return response.json()
        return None
    except Exception as e:
        logger.error(f"Open Targets API Error: {e}")
        return None


def enrich_graph_with_opentargets(g, max_genes=5):
    """
    Scans the graph for causal Genes, queries Open Targets,
    and adds the updated schema data to the graph.
    """
    logger.info("--- Starting Open Targets Integration ---")

    # Find all Genes currently in the graph
    gene_nodes = list(g.subjects(RDF.type, MTR.Gene))
    logger.info(f"Found {len(gene_nodes)} causal Genes in the local graph.")

    genes_to_process = gene_nodes[:max_genes]

    for i, gene_node in enumerate(genes_to_process):
        gene_symbol = str(gene_node).split('/')[-1].strip()

        logger.info(f"[{i + 1}/{len(genes_to_process)}] Fetching OT data for {gene_symbol}...")

        data = fetch_opentargets_data(gene_symbol)

        if data and 'data' in data and data['data']['search']['hits']:
            target_obj = data['data']['search']['hits'][0]['object']

            # ==========================================
            # A. PARSE UPDATED TRACTABILITY
            # ==========================================
            tractability = target_obj.get('tractability', [])
            if tractability:
                for t in tractability:
                    # We look for Small Molecule (SM) or Antibody (AB) tractability where value is True
                    if t.get('value') is True:
                        if t.get('modality') == 'SM':
                            g.add((gene_node, MTR.small_molecule_tractable, Literal(True, datatype=XSD.boolean)))
                        elif t.get('modality') == 'AB':
                            g.add((gene_node, MTR.antibody_tractable, Literal(True, datatype=XSD.boolean)))

            # ==========================================
            # B. PARSE UPDATED SAFETY LIABILITIES
            # ==========================================
            safety_liabilities = target_obj.get('safetyLiabilities')
            if safety_liabilities:
                logger.info(f"   -> Found {len(safety_liabilities)} safety liabilities!")
                for liability in safety_liabilities:
                    event_name = liability.get('event')
                    event_id = liability.get('eventId')
                    if not event_name: continue

                    # Clean the ID (Use official EFO if provided, else create a custom one)
                    if event_id and event_id.startswith('EFO_'):
                        event_node = EFO[event_id]
                    else:
                        safe_event_id = urllib.parse.quote(event_name.strip().replace(" ", "_").replace('"', ''))
                        event_node = EFO[f"Custom_Event_{safe_event_id}"]

                    g.add((event_node, RDF.type, BIOLINK.PhenotypicFeature))
                    g.add((event_node, RDFS.label, Literal(event_name)))
                    g.add((gene_node, BIOLINK.has_adverse_event, event_node))

                    # Process the new 'biosamples' array
                    biosamples = liability.get('biosamples', [])
                    for bio in biosamples:
                        tissue = bio.get('tissueLabel')
                        if tissue:
                            g.add((event_node, MTR.tissue_context, Literal(tissue)))

            # ==========================================
            # C. PARSE KNOWN DRUGS
            # ==========================================
            if 'knownDrugs' in target_obj and target_obj['knownDrugs']:
                drug_rows = target_obj['knownDrugs'].get('rows', [])
                if drug_rows:
                    logger.info(f"   -> Found {len(drug_rows)} known drugs!")
                    for row in drug_rows:
                        drug_info = row['drug']
                        drug_node = CHEMBL[drug_info['id']]

                        g.add((drug_node, RDF.type, BIOLINK.Drug))
                        g.add((drug_node, RDFS.label, Literal(drug_info['name'])))
                        g.add((drug_node, MTR.max_clinical_phase, Literal(row['phase'], datatype=XSD.float)))
                        g.add((drug_node, BIOLINK.targets, gene_node))
                        g.add((drug_node, PROV.wasDerivedFrom, URIRef("https://platform.opentargets.org/")))

            # ==========================================
            # D. PARSE NEW: ASSOCIATED DISEASES
            # ==========================================
            assoc_diseases = target_obj.get('associatedDiseases')
            if assoc_diseases and assoc_diseases.get('rows'):
                disease_rows = assoc_diseases['rows']
                logger.info(f"   -> Found {len(disease_rows)} top associated diseases!")

                for row in disease_rows:
                    disease_info = row.get('disease', {})
                    disease_id = disease_info.get('id')
                    disease_name = disease_info.get('name')

                    if disease_id:
                        # Convert Open Targets EFO/MONDO IDs to valid URIs
                        disease_node = EFO[disease_id]

                        g.add((disease_node, RDF.type, BIOLINK.Disease))
                        g.add((disease_node, RDFS.label, Literal(disease_name)))

                        # Create an OBAN Association Node linking Gene -> Disease
                        assoc_node = MTR[f"OT_Assoc_{gene_symbol}_{disease_id}"]
                        # Ensure the gene has a label (Fixing the previous SHACL error!)
                        g.add((gene_node, RDFS.label, Literal(gene_symbol, datatype=XSD.string)))

                        g.add((assoc_node, RDF.type, BIOLINK.Association))
                        g.add((assoc_node, BIOLINK.has_subject, gene_node))
                        g.add((assoc_node, BIOLINK.has_object, disease_node))
                        g.add((assoc_node, PROV.wasDerivedFrom, URIRef("https://platform.opentargets.org/")))

                        # Extract the 'genetics' datasource score if it exists
                        for ds in row.get('datasourceScores', []):
                            if ds.get('id') == 'genetics':
                                g.add((assoc_node, MTR.ot_genetics_score, Literal(ds.get('score'), datatype=XSD.float)))

        else:
            logger.warning(f"   -> {gene_symbol} not found in Open Targets database.")

        time.sleep(0.1)  # Be polite to the API

    logger.info("--- Open Targets Integration Complete ---")
    return g