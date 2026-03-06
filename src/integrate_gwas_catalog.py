import requests
import time
import urllib.parse
from rdflib import Graph, Literal, RDF, URIRef, Namespace
from rdflib.namespace import RDFS, XSD
from schema_definition import *


BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"

def gwas_api_request(endpoint, params=None):
    full_url = f"{BASE_URL.rstrip('/')}/{endpoint.lstrip('/')}"
    try:
        response = requests.get(full_url, params=params, headers={"Accept": "application/json"})
        if response.status_code == 200:
            return response.json()
        # Handle 404s gracefully without crashing
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

    # Regex to find the first number (handles decimals and signs)
    match = re.search(r"[-+]?\d*\.?\d+", beta_str)

    if match:
        numeric_val = float(match.group())
        # Extract the rest of the string as the unit/direction
        description = beta_str.replace(match.group(), "").strip()

        # CRITICAL FIX: If the text says "decrease", the statistical beta is actually negative
        if "decrease" in description.lower() and numeric_val > 0:
            numeric_val = -numeric_val

        return numeric_val, description

    return None, beta_str # Fallback if no numbers are found

def enrich_graph_with_gwas(g, max_snps=5):
    print("\n--- Starting GWAS Catalog Integration (v2 API) ---")



    snp_nodes = list(g.subjects(RDF.type, MTR.SNP))
    snps_to_process = snp_nodes[:max_snps]

    for i, snp_node in enumerate(snps_to_process):
        # .strip() prevents 404 errors caused by hidden whitespace
        rsid = str(snp_node).split('/')[-1].strip()

        print(f"[{i+1}/{len(snps_to_process)}] Fetching v2 associations for {rsid}...")

        endpoint = "v2/associations"
        params = {
            "rs_id": rsid,
            "sort": "p_value",
            "direction": "asc",
            "size": 10
        }

        associations_data = gwas_api_request(endpoint, params)

        # --- ROBUST ERROR HANDLING ---
        # Check if data exists AND if the nested keys we expect are actually there
        if not associations_data or '_embedded' not in associations_data or 'associations' not in associations_data['_embedded']:
            print("   -> No associated traits found in the GWAS Catalog for this SNP.")
            time.sleep(0.1)
            continue # Safely skip to the next SNP

        associations_list = associations_data['_embedded']['associations']

        # Double check that the list isn't empty
        if len(associations_list) == 0:
            print("   -> No associated traits found in the GWAS Catalog for this SNP.")
            time.sleep(0.1)
            continue

        print(f"   -> Found {len(associations_list)} significant associations!")

        for assoc in associations_list:
            # Extract core fields
            p_value = assoc.get('p_value')
            beta = assoc.get('beta')
            ci_lower = assoc.get('ci_lower')
            ci_upper = assoc.get('ci_upper')
            effect_allele = assoc.get('snp_effect_allele')

            # Extract provenance fields
            pubmed_id = assoc.get('pubmed_id')
            accession_id = assoc.get('accession_id')
            reported_trait = assoc.get('reported_trait')

            efo_traits = assoc.get('efo_traits', [])

            for trait_dict in efo_traits:
                trait_name = trait_dict.get('efo_trait')
                if not trait_name:
                    continue

                # Create Trait Node
                safe_trait_id = urllib.parse.quote(trait_name.strip().replace(" ", "_").replace('"', ''))
                trait_node = EFO[f"Custom_{safe_trait_id}"]

                # Add Trait Core Identity
                g.add((trait_node, RDF.type, BIOLINK.DiseaseOrPhenotypicFeature))
                g.add((trait_node, RDFS.label, Literal(trait_name)))

                # If author reported a different name, save it as an alternative label!
                if reported_trait and reported_trait != trait_name:
                    g.add((trait_node, SKOS.altLabel, Literal(reported_trait)))

                # Create the Association Node
                # We use the accession_id to ensure every unique GWAS study gets its own distinct edge
                unique_assoc_id = accession_id if accession_id else f"Assoc_{rsid}_{safe_trait_id}"
                assoc_node = MTR[f"GWAS_{unique_assoc_id}"]

                g.add((assoc_node, RDF.type, BIOLINK.Association))
                g.add((assoc_node, BIOLINK.has_subject, snp_node))
                g.add((assoc_node, BIOLINK.has_object, trait_node))

                # --- ADDING THE NEW STATISTICAL FIELDS ---
        # --- ADDING THE NEW STATISTICAL FIELDS ---
                if p_value is not None:
                    g.add((assoc_node, MTR.p_value, Literal(p_value, datatype=XSD.float)))

                # --- NEW BETA PARSING LOGIC ---
                if beta is not None:
                    numeric_beta, beta_unit = parse_gwas_beta(beta)

                    if numeric_beta is not None:
                        # Store the mathematically pure float
                        g.add((assoc_node, MTR.beta, Literal(numeric_beta, datatype=XSD.float)))

                    if beta_unit:
                        # Store the unit/context as a separate string so you don't lose the biological meaning!
                        g.add((assoc_node, MTR.beta_unit, Literal(beta_unit)))
                # ------------------------------

                if ci_lower is not None and ci_upper is not None:
                    g.add((assoc_node, MTR.ci_lower, Literal(ci_lower, datatype=XSD.float)))
                    g.add((assoc_node, MTR.ci_upper, Literal(ci_upper, datatype=XSD.float)))
                if effect_allele:
                    g.add((assoc_node, MTR.effect_allele, Literal(effect_allele)))

                # --- ADDING THE NEW PROVENANCE FIELDS ---
                g.add((assoc_node, PROV.wasGeneratedBy, Literal("EBI_GWAS_REST_API_v2")))

                if pubmed_id:
                    # Link the association to the exact PubMed article
                    paper_node = PUBMED[str(pubmed_id)]
                    g.add((assoc_node, BIOLINK.publications, paper_node))

                if accession_id:
                    # Link the association to the specific GWAS Catalog Study
                    study_node = URIRef(f"https://www.ebi.ac.uk/gwas/studies/{accession_id}")
                    g.add((assoc_node, PROV.wasDerivedFrom, study_node))

        time.sleep(0.1)

    print("\n--- GWAS Integration Complete ---")
    return g