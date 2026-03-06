import json
import urllib.parse
from rdflib import Graph, Literal, RDF, URIRef, Namespace
from rdflib.namespace import RDFS, XSD, SKOS
from schema_definition import *

def clean_uri_string(text):
    """Utility to clean strings so they can be valid URIs."""
    return urllib.parse.quote(str(text).strip().replace(" ", "_").replace('"', ''))


def add_safe_literal(graph, subject, predicate, value, datatype=None):
    """Helper to add triples only if the value is not 'NA', None, or empty."""
    if value not in [None, "NA", "", "NaN"]:
        if datatype:
            try:
                # Ensure float/int casting works before adding
                if datatype == XSD.float:
                    val = float(value)
                elif datatype == XSD.integer:
                    val = int(float(value))  # Handle cases like '1.0' for ints
                else:
                    val = value
                graph.add((subject, predicate, Literal(val, datatype=datatype)))
            except (ValueError, TypeError):
                pass  # Skip if datatype conversion fails
        else:
            graph.add((subject, predicate, Literal(value)))


def add_rqtl_to_graph(json_file_path, g):

    with open(json_file_path, 'r') as f:
        data = json.load(f)

    for entry in data:
        ratio_id = entry['ratio_accession']
        ratio_node = MTR[ratio_id]

        # ==========================================
        # 1. DEFINE RATIO NODE & TOP-LEVEL SUMMARIES
        # ==========================================
        g.add((ratio_node, RDF.type, MTR.MetaboliteRatio))
        g.add((ratio_node, RDFS.label, Literal(entry['ratio_name'])))

        add_safe_literal(g, ratio_node, MTR.reaction_distance, entry.get('reaction_distance'), XSD.integer)
        add_safe_literal(g, ratio_node, MTR.max_pgain, entry.get('max_pgain'), XSD.float)

        # 1b. Add LLM Explanations as Annotations
        llm_raw = entry.get('llm_response')
        llm_resp = {}  # Default to empty dict

        # Only try to dig deeper if llm_raw is actually a dictionary
        if isinstance(llm_raw, dict):
            # Note: Using 'reponse' exactly as it is spelled in your JSON schema
            potential_resp = llm_raw.get('reponse', {})
            if isinstance(potential_resp, dict):
                llm_resp = potential_resp
        add_safe_literal(g, ratio_node, MTR.llm_ratio_explanation, llm_resp.get('ratio_explanation'))
        add_safe_literal(g, ratio_node, MTR.llm_ratio_evidence, llm_resp.get('ratio_evidence'))
        add_safe_literal(g, ratio_node, MTR.llm_phenotype_driver, llm_resp.get('phenotype_driver'))
        add_safe_literal(g, ratio_node, MTR.llm_gene_phenotype_relationship,
                         llm_resp.get('gene_ratio_phenotype_relationship'))

        # ==========================================
        # 2. NUMERATOR & DENOMINATOR METABOLITES
        # ==========================================
        metabolite_nodes = {}  # Store for later mQTL linking

        for side, relation in [('numerator_metabolite', MTR.has_numerator),
                               ('denominator_metabolite', MTR.has_denominator)]:
            metab = entry.get(side, {})
            if metab and metab.get('chebi') != "NA":
                chebi_id = str(int(float(metab['chebi'])))
                metab_node = CHEBI[chebi_id]
                metabolite_nodes[side] = metab_node  # Save for mQTL processing

                g.add((metab_node, RDF.type, MTR.Metabolite))
                g.add((metab_node, RDFS.label, Literal(metab.get('name', 'Unknown'))))

                # Link Cross-References
                if metab.get('hmdb') != "NA":
                    g.add((metab_node, SKOS.exactMatch, HMDB[metab['hmdb']]))
                if metab.get('kegg') != "NA":
                    g.add((metab_node, SKOS.exactMatch, MTR[f"KEGG_{metab['kegg']}"]))
                if metab.get('inchikey') != "NA":
                    g.add((metab_node, MTR.inchikey, Literal(metab['inchikey'])))

                # Add LLM specific explanations to the metabolites
                if side == 'numerator_metabolite':
                    add_safe_literal(g, metab_node, MTR.llm_explanation,
                                     llm_resp.get('numerator_metabolite_explanation'))
                else:
                    add_safe_literal(g, metab_node, MTR.llm_explanation,
                                     llm_resp.get('denominator_metabolite_explanation'))

                g.add((ratio_node, relation, metab_node))

        # ==========================================
        # 3. PROCESS SNPs AND ASSOCIATED REGIONS
        # ==========================================
        regions = entry.get('associated_regions', {})
        for rsid, details in regions.items():
            snp_node = DBSNP[rsid]
            g.add((snp_node, RDF.type, MTR.SNP))

            # OBAN Association Node (SNP -> Ratio)
            assoc_node = MTR[f"Assoc_{rsid}_{ratio_id}"]
            g.add((assoc_node, RDF.type, OBAN.association))
            g.add((assoc_node, OBAN.has_subject, snp_node))
            g.add((assoc_node, OBAN.has_object, ratio_node))

            # Extract full Summary Statistics
            stats = details.get('summary_statistics', {})
            add_safe_literal(g, assoc_node, MTR.p_value, stats.get('p_value'), XSD.float)
            add_safe_literal(g, assoc_node, MTR.beta, stats.get('beta'), XSD.float)
            add_safe_literal(g, assoc_node, MTR.se, stats.get('se'), XSD.float)
            add_safe_literal(g, assoc_node, MTR.z_score, stats.get('z'), XSD.float)
            add_safe_literal(g, assoc_node, MTR.maf, stats.get('maf'), XSD.float)
            add_safe_literal(g, assoc_node, MTR.pgain, stats.get('pgain'), XSD.float)
            add_safe_literal(g, assoc_node, MTR.effect_allele, stats.get('effect_allele'))
            add_safe_literal(g, assoc_node, MTR.reference_allele, stats.get('reference_allele'))
            add_safe_literal(g, assoc_node, MTR.numerator_driven, stats.get('numerator_driven'), XSD.boolean)
            add_safe_literal(g, assoc_node, MTR.denominator_driven, stats.get('denominator_driven'), XSD.boolean)

            # Link Closest Genes
            for gene in stats.get('closest_genes', []):
                gene_node = HGNC[clean_uri_string(gene)]
                g.add((snp_node, MTR.closest_gene, gene_node))

            # ==========================================
            # 4. PROCESS MISSENSE VARIANTS IN LD
            # ==========================================
            for missense in details.get('missense_variants_in_ld', []):
                ms_rsid = missense.get('missense_variant')
                if ms_rsid and ms_rsid != "NA":
                    ms_node = DBSNP[ms_rsid]
                    g.add((ms_node, RDF.type, MTR.MissenseVariant))
                    g.add((snp_node, MTR.in_ld_with, ms_node))
                    add_safe_literal(g, ms_node, MTR.ld_score, missense.get('ld'), XSD.float)
                    add_safe_literal(g, ms_node, MTR.consequence, missense.get('Consequence'))

                    if missense.get('gene_name') != "NA":
                        ms_gene_node = HGNC[clean_uri_string(missense['gene_name'])]
                        g.add((ms_node, MTR.affects_gene, ms_gene_node))

            # ==========================================
            # 5. PROCESS CNVs (Copy Number Variants)
            # ==========================================
            for cnv in details.get('cnv_chiara', []):
                cnv_id = cnv.get('ID')
                if cnv_id:
                    cnv_node = MTR[f"CNV_{clean_uri_string(cnv_id)}"]
                    g.add((cnv_node, RDF.type, MTR.CopyNumberVariant))
                    g.add((snp_node, MTR.co_located_with_cnv, cnv_node))  # Link CNV to SNP region

                    add_safe_literal(g, cnv_node, MTR.chromosome, cnv.get('CHR'), XSD.integer)
                    add_safe_literal(g, cnv_node, MTR.freq_cnv, cnv.get('FreqCNV'), XSD.float)
                    add_safe_literal(g, cnv_node, MTR.num_del, cnv.get('NumDel'), XSD.integer)
                    add_safe_literal(g, cnv_node, MTR.num_dup, cnv.get('NumDup'), XSD.integer)

            # ==========================================
            # 6. MASTER LOOP: CAUSAL EVIDENCE (MR & COLOC)
            # ==========================================
            # We loop over all data sources (eQTL, pQTL, GTEx, complex diseases)
            causal_datasets = ['eQTL', 'pQTL', 'GTEx', 'finngen', 'mrlink2_traits', 'yang', 'rosmap', 'western']
            # We loop over the outcomes (rQTLs vs mQTLs)
            target_outcomes = ['ratio', 'numerator_metabolite', 'denominator_metabolite']

            for dataset_name in causal_datasets:
                dataset_data = details.get(dataset_name, {})

                for target_level in target_outcomes:
                    causal_tests = dataset_data.get(target_level, [])

                    for test in causal_tests:
                        exposure_name = test.get('exposure')
                        if not exposure_name or exposure_name == "NA":
                            continue

                        # A. Determine the Outcome Node (Is this acting on the Ratio or an individual Metabolite?)
                        if target_level == 'ratio':
                            outcome_node = ratio_node
                        else:
                            # Use the previously saved metabolite nodes
                            outcome_node = metabolite_nodes.get(target_level)
                            if not outcome_node:
                                continue  # Skip if that metabolite doesn't exist

                        # B. Determine the Exposure Node Type
                        safe_exp = clean_uri_string(exposure_name)
                        if dataset_name in ['eQTL', 'GTEx', 'rosmap']:
                            exposure_node = HGNC[safe_exp]
                            g.add((exposure_node, RDF.type, MTR.Gene))
                        elif dataset_name in ['pQTL', 'western']:
                            exposure_node = MTR[f"Protein_{safe_exp}"]
                            g.add((exposure_node, RDF.type, MTR.Protein))
                        else:
                            # FinnGen, MR-Link, Yang (These are typically complex traits/diseases)
                            exposure_node = MTR[f"Trait_{safe_exp}"]
                            g.add((exposure_node, RDF.type, MTR.ComplexTrait))

                        g.add((exposure_node, RDFS.label, Literal(exposure_name)))

                        # C. Create the Causal Assessment Node
                        assessment_id = f"Causal_{dataset_name}_{safe_exp}_to_{target_level}_{ratio_id}"
                        assessment_node = MTR[assessment_id]

                        g.add((assessment_node, RDF.type, MTR.CausalAssessment))
                        g.add((assessment_node, MTR.dataset_source, Literal(dataset_name)))
                        g.add((assessment_node, MTR.target_level, Literal(target_level)))  # Tags it as mQTL or rQTL
                        g.add((assessment_node, MTR.has_exposure, exposure_node))
                        g.add((assessment_node, MTR.has_outcome, outcome_node))

                        # D. Attach all MR and Colocalization Statistics
                        # IVW (Inverse Variance Weighted) is standard MR
                        add_safe_literal(g, assessment_node, MTR.beta_ivw, test.get('beta_ivw'), XSD.float)
                        add_safe_literal(g, assessment_node, MTR.se_ivw, test.get('se_ivw'), XSD.float)
                        add_safe_literal(g, assessment_node, MTR.p_ivw, test.get('p_ivw'), XSD.float)

                        # PCA and Alpha stats
                        add_safe_literal(g, assessment_node, MTR.beta_pca, test.get('beta_pca'), XSD.float)
                        add_safe_literal(g, assessment_node, MTR.p_pca, test.get('p_pca'), XSD.float)
                        add_safe_literal(g, assessment_node, MTR.alpha, test.get('alpha'), XSD.float)

                        # Colocalization Probabilities (H4 = Both traits share a single causal variant)
                        add_safe_literal(g, assessment_node, MTR.coloc_h3, test.get('PP.H3.abf'), XSD.float)
                        add_safe_literal(g, assessment_node, MTR.coloc_h4, test.get('PP.H4.abf'), XSD.float)
                        add_safe_literal(g, assessment_node, MTR.susie_h4, test.get('susie_max_PP.H4.abf'), XSD.float)

                        # Contextual metadata
                        add_safe_literal(g, assessment_node, MTR.measured_in_tissue, test.get('tissue'))
                        add_safe_literal(g, assessment_node, MTR.snps_overlap, test.get('m_snps_overlap'), XSD.integer)
                        add_safe_literal(g, assessment_node, MTR.variance_explained, test.get('var_explained'),
                                         XSD.float)

    return g