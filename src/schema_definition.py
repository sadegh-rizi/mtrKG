from rdflib import Graph, Literal, RDF, URIRef, Namespace
from rdflib.namespace import RDFS, XSD, OWL, SKOS

# =========================
# Namespaces
# =========================
MTR = Namespace("https://metabolite-ratio-app.unil.ch/rqtl/")
BIOLINK = Namespace("https://w3id.org/biolink/vocab/")
EFO = Namespace("http://www.ebi.ac.uk/efo/")
SIO = Namespace("http://semanticscience.org/resource/")
PROV = Namespace("http://www.w3.org/ns/prov#")
CHEMBL = Namespace("http://identifiers.org/chembl.compound/")
CHEBI = Namespace("http://purl.obolibrary.org/obo/CHEBI_")
HMDB = Namespace("http://identifiers.org/hmdb/")
DBSNP = Namespace("http://identifiers.org/dbsnp/")
REACTOME = Namespace("https://reactome.org/content/detail/")
PUBMED = Namespace("https://pubmed.ncbi.nlm.nih.gov/")
SO = Namespace("http://purl.obolibrary.org/obo/SO_")
SH = Namespace("http://www.w3.org/ns/shacl#")
HGNC = Namespace("http://identifiers.org/hgnc.symbol/")



def bind_namespaces(g: Graph) -> Graph:
    g.bind("mtr", MTR)
    g.bind("biolink", BIOLINK)
    g.bind("efo", EFO)
    g.bind("sio", SIO)
    g.bind("prov", PROV)
    g.bind("chembl", CHEMBL)
    g.bind("chebi", CHEBI)
    g.bind("hmdb", HMDB)
    g.bind("dbsnp", DBSNP)
    g.bind("reactome", REACTOME)
    g.bind("pubmed", PUBMED)
    g.bind("so", SO)
    g.bind("skos", SKOS)
    g.bind("sh", SH)
    return g


def declare_class(g, uri, label, comment=None, superclass=None):
    g.add((uri, RDF.type, OWL.Class))
    g.add((uri, RDFS.label, Literal(label)))
    if comment:
        g.add((uri, RDFS.comment, Literal(comment)))
    if superclass:
        g.add((uri, RDFS.subClassOf, superclass))


def declare_object_property(g, uri, label, domain=None, range_=None, comment=None, superprop=None):
    g.add((uri, RDF.type, OWL.ObjectProperty))
    g.add((uri, RDFS.label, Literal(label)))
    if comment:
        g.add((uri, RDFS.comment, Literal(comment)))
    if domain:
        g.add((uri, RDFS.domain, domain))
    if range_:
        g.add((uri, RDFS.range, range_))
    if superprop:
        g.add((uri, RDFS.subPropertyOf, superprop))


def declare_datatype_property(g, uri, label, domain=None, range_=None, comment=None, superprop=None):
    g.add((uri, RDF.type, OWL.DatatypeProperty))
    g.add((uri, RDFS.label, Literal(label)))
    if comment:
        g.add((uri, RDFS.comment, Literal(comment)))
    if domain:
        g.add((uri, RDFS.domain, domain))
    if range_:
        g.add((uri, RDFS.range, range_))
    if superprop:
        g.add((uri, RDFS.subPropertyOf, superprop))


def build_schema(g: Graph) -> Graph:
    bind_namespaces(g)

    ontology_uri = URIRef("https://metabolite-ratio-app.unil.ch/ontology")
    g.add((ontology_uri, RDF.type, OWL.Ontology))
    g.add((ontology_uri, RDFS.label, Literal("rQTL Metabolic Ratio Ontology")))
    g.add((ontology_uri, RDFS.comment, Literal(
        "Schema for rQTL metabolite ratios, variants, genes, molecular products, and causal evidence."
    )))

    # =========================
    # Classes
    # =========================
    declare_class(g, MTR.MetaboliteRatio, "Metabolite Ratio", "A ratio between two metabolite concentrations.")
    declare_class(g, MTR.Metabolite, "Metabolite", "A small chemical compound.", BIOLINK.ChemicalEntity)
    declare_class(g, MTR.Gene, "Gene", "A genomic locus.", BIOLINK.Gene)
    declare_class(g, MTR.Transcript, "Transcript", "An RNA transcript.", BIOLINK.Transcript)
    declare_class(g, MTR.Protein, "Protein", "A translated gene product.", BIOLINK.Protein)
    declare_class(g, MTR.Phenotype, "Phenotype", "A disease or phenotypic feature.", BIOLINK.DiseaseOrPhenotypicFeature)

    declare_class(g, MTR.SNP, "SNP", "A single nucleotide polymorphism.", BIOLINK.SequenceVariant)
    g.add((MTR.SNP, RDFS.subClassOf, SO["0000694"]))

    declare_class(g, MTR.MissenseVariant, "Missense Variant", "Missense sequence variant.", BIOLINK.SequenceVariant)
    declare_class(g, MTR.CopyNumberVariant, "Copy Number Variant", "Copy number variant.", BIOLINK.SequenceVariant)

    declare_class(
        g,
        MTR.VariantToMetaboliteRatioAssociation,
        "Variant to Metabolite Ratio Association",
        "Association between a variant and a metabolite ratio.",
        BIOLINK.Association,
    )
    declare_class(
        g,
        MTR.CausalAssessment,
        "Causal Assessment",
        "MR and colocalization-based causal assessment.",
        BIOLINK.Association,
    )

    declare_class(g, MTR.DatasetRelease, "Dataset Release", "Generated integrated dataset.", PROV.Entity)
    declare_class(g, MTR.IntegrationRun, "Integration Run", "RDF integration activity.", PROV.Activity)

    # =========================
    # Object properties
    # =========================
    declare_object_property(g, MTR.has_numerator, "has numerator", MTR.MetaboliteRatio, MTR.Metabolite)
    declare_object_property(g, MTR.has_denominator, "has denominator", MTR.MetaboliteRatio, MTR.Metabolite)

    declare_object_property(g, MTR.closest_gene, "closest gene", BIOLINK.SequenceVariant, MTR.Gene)
    declare_object_property(g, MTR.in_ld_with, "in LD with", BIOLINK.SequenceVariant, BIOLINK.SequenceVariant)
    declare_object_property(g, MTR.affects_gene, "affects gene", BIOLINK.SequenceVariant, MTR.Gene)
    declare_object_property(g, MTR.co_located_with_cnv, "co-located with CNV", BIOLINK.SequenceVariant, MTR.CopyNumberVariant)

    declare_object_property(g, MTR.associated_with_ratio, "associated with ratio", BIOLINK.SequenceVariant, MTR.MetaboliteRatio)
    declare_object_property(g, MTR.causal_influence_on, "causal influence on")

    declare_object_property(g, MTR.has_exposure, "has exposure", MTR.CausalAssessment)
    declare_object_property(g, MTR.has_outcome, "has outcome", MTR.CausalAssessment)
    declare_object_property(g, MTR.associated_variant, "associated variant", MTR.CausalAssessment, BIOLINK.SequenceVariant)

    # =========================
    # Metabolite identifiers
    # =========================
    declare_datatype_property(g, MTR.hmdb_id, "HMDB ID", MTR.Metabolite, XSD.string)
    declare_datatype_property(g, MTR.kegg_id, "KEGG ID", MTR.Metabolite, XSD.string)
    declare_datatype_property(g, MTR.metabolon_id, "Metabolon ID", MTR.Metabolite, XSD.string)
    declare_datatype_property(g, MTR.local_accession, "local accession", MTR.Metabolite, XSD.string)
    declare_datatype_property(g, MTR.inchikey, "InChIKey", MTR.Metabolite, XSD.string)

    # =========================
    # Ratio properties
    # =========================
    declare_datatype_property(g, MTR.ratio_formula, "ratio formula", MTR.MetaboliteRatio, XSD.string)
    declare_datatype_property(g, MTR.reaction_distance, "reaction distance", MTR.MetaboliteRatio, XSD.integer)
    declare_datatype_property(g, MTR.max_pgain, "max p-gain", MTR.MetaboliteRatio, XSD.double)

    # LLM annotations
    declare_datatype_property(g, MTR.llm_ratio_explanation, "LLM ratio explanation", MTR.MetaboliteRatio, XSD.string)
    declare_datatype_property(g, MTR.llm_ratio_evidence, "LLM ratio evidence", MTR.MetaboliteRatio, XSD.string)
    declare_datatype_property(g, MTR.llm_phenotype_driver, "LLM phenotype driver", MTR.MetaboliteRatio, XSD.string)
    declare_datatype_property(g, MTR.llm_gene_phenotype_relationship, "LLM gene-phenotype relationship", MTR.MetaboliteRatio, XSD.string)
    declare_datatype_property(g, MTR.llm_explanation, "LLM explanation", MTR.Metabolite, XSD.string)

    # =========================
    # Variant summary statistics
    # =========================
    declare_datatype_property(g, MTR.chromosome, "chromosome", BIOLINK.SequenceVariant, XSD.string)
    declare_datatype_property(g, MTR.position, "position", BIOLINK.SequenceVariant, XSD.integer)
    declare_datatype_property(g, MTR.pos_name, "position name", BIOLINK.SequenceVariant, XSD.string)

    declare_datatype_property(g, MTR.effect_allele, "effect allele", BIOLINK.Association, XSD.string)
    declare_datatype_property(g, MTR.reference_allele, "reference allele", BIOLINK.Association, XSD.string)
    declare_datatype_property(g, MTR.maf, "minor allele frequency", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.beta, "beta", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.standard_error, "standard error", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.z_score, "z-score", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.log_p_value, "log p-value", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.beta_1, "beta numerator", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.se_1, "SE numerator", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.z_1, "z numerator", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.beta_2, "beta denominator", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.se_2, "SE denominator", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.z_2, "z denominator", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.numerator_driven, "numerator driven", BIOLINK.Association, XSD.boolean)
    declare_datatype_property(g, MTR.denominator_driven, "denominator driven", BIOLINK.Association, XSD.boolean)
    declare_datatype_property(g, MTR.already_found, "already found", BIOLINK.Association, XSD.boolean)
    declare_datatype_property(g, MTR.pgain, "p-gain", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.log_pgain, "log p-gain", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.infomap_cluster, "Infomap cluster", BIOLINK.Association, XSD.integer)

    # Missense / LD
    declare_datatype_property(g, MTR.ld_score, "LD score", BIOLINK.SequenceVariant, XSD.double)
    declare_datatype_property(g, MTR.variant_location, "variant location", BIOLINK.SequenceVariant, XSD.string)
    declare_datatype_property(g, MTR.allele_change, "allele change", BIOLINK.SequenceVariant, XSD.string)
    declare_datatype_property(g, MTR.consequence, "consequence", BIOLINK.SequenceVariant, XSD.string)
    declare_datatype_property(g, MTR.codons, "codons", BIOLINK.SequenceVariant, XSD.string)

    # CNV
    declare_datatype_property(g, MTR.num_cnv, "number of CNVs", MTR.CopyNumberVariant, XSD.integer)
    declare_datatype_property(g, MTR.num_dup, "number of duplications", MTR.CopyNumberVariant, XSD.integer)
    declare_datatype_property(g, MTR.num_del, "number of deletions", MTR.CopyNumberVariant, XSD.integer)
    declare_datatype_property(g, MTR.num_neutral, "number of neutral", MTR.CopyNumberVariant, XSD.integer)
    declare_datatype_property(g, MTR.freq_cnv, "CNV frequency", MTR.CopyNumberVariant, XSD.double)
    declare_datatype_property(g, MTR.freq_dup, "duplication frequency", MTR.CopyNumberVariant, XSD.double)
    declare_datatype_property(g, MTR.freq_del, "deletion frequency", MTR.CopyNumberVariant, XSD.double)

    # =========================
    # Causal assessment metadata
    # =========================
    declare_datatype_property(g, MTR.dataset_source, "dataset source", MTR.CausalAssessment, XSD.string)
    declare_datatype_property(g, MTR.target_level, "target level", MTR.CausalAssessment, XSD.string)
    declare_datatype_property(g, MTR.causal_direction, "causal direction", MTR.CausalAssessment, XSD.string)

    declare_datatype_property(g, MTR.exposure_file, "exposure file", MTR.CausalAssessment, XSD.string)
    declare_datatype_property(g, MTR.outcome_file, "outcome file", MTR.CausalAssessment, XSD.string)
    declare_datatype_property(g, MTR.region, "region", MTR.CausalAssessment, XSD.string)
    declare_datatype_property(g, MTR.measured_in_tissue, "measured in tissue", MTR.CausalAssessment, XSD.string)

    declare_datatype_property(g, MTR.variance_explained, "variance explained", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.snps_overlap, "SNPs overlap", MTR.CausalAssessment, XSD.integer)

    declare_datatype_property(g, MTR.alpha, "alpha", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.se_alpha, "SE alpha", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.p_alpha, "p alpha", MTR.CausalAssessment, XSD.double)

    declare_datatype_property(g, MTR.sigma_y, "sigma y", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.se_sigma_y, "SE sigma y", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.p_sigma_y, "p sigma y", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.sigma_x, "sigma x", MTR.CausalAssessment, XSD.double)

    declare_datatype_property(g, MTR.beta_ivw, "beta IVW", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.se_ivw, "SE IVW", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.p_ivw, "p IVW", MTR.CausalAssessment, XSD.double)

    declare_datatype_property(g, MTR.beta_ivw_r, "beta IVW reverse", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.se_ivw_r, "SE IVW reverse", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.p_ivw_r, "p IVW reverse", MTR.CausalAssessment, XSD.double)

    declare_datatype_property(g, MTR.beta_pca, "beta PCA", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.se_pca, "SE PCA", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.p_pca, "p PCA", MTR.CausalAssessment, XSD.double)

    declare_datatype_property(g, MTR.coloc_h1, "coloc H1", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.coloc_h2, "coloc H2", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.coloc_h3, "coloc H3", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.coloc_h4, "coloc H4", MTR.CausalAssessment, XSD.double)

    declare_datatype_property(g, MTR.susie_h1, "SuSiE H1", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.susie_h2, "SuSiE H2", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.susie_h3, "SuSiE H3", MTR.CausalAssessment, XSD.double)
    declare_datatype_property(g, MTR.susie_h4, "SuSiE H4", MTR.CausalAssessment, XSD.double)

    # GWAS-specific datatype properties
    declare_datatype_property(g, MTR.beta_unit, "beta unit", BIOLINK.Association, XSD.string)
    declare_datatype_property(g, MTR.ci_lower, "confidence interval lower bound", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.ci_upper, "confidence interval upper bound", BIOLINK.Association, XSD.double)

    # If not already present
    declare_datatype_property(g, MTR.beta, "beta", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.effect_allele, "effect allele", BIOLINK.Association, XSD.string)

    # Reactome-specific datatype properties
    declare_datatype_property(g, MTR.reactome_source_type, "Reactome source type", BIOLINK.Association, XSD.string)
    declare_datatype_property(g, MTR.reactome_analysis_token, "Reactome analysis token", BIOLINK.Association,
                              XSD.string)
    declare_datatype_property(g, MTR.reactome_species, "Reactome species", BIOLINK.Association, XSD.string)

    declare_datatype_property(g, MTR.reactome_entities_p_value, "Reactome entities p-value", BIOLINK.Association,
                              XSD.double)
    declare_datatype_property(g, MTR.reactome_entities_fdr, "Reactome entities FDR", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.reactome_entities_found, "Reactome entities found", BIOLINK.Association,
                              XSD.integer)
    declare_datatype_property(g, MTR.reactome_entities_total, "Reactome entities total", BIOLINK.Association,
                              XSD.integer)

    declare_datatype_property(g, MTR.reactome_reactions_found, "Reactome reactions found", BIOLINK.Association,
                              XSD.integer)
    declare_datatype_property(g, MTR.reactome_reactions_total, "Reactome reactions total", BIOLINK.Association,
                              XSD.integer)

    declare_class(
        g,
        MTR.MetaboliteLocationAssociation,
        "Metabolite Location Association",
        "Association between a metabolite and a biological location asserted from HMDB.",
        BIOLINK.Association,
    )

    declare_datatype_property(
        g,
        MTR.source_dataset,
        "source dataset",
        BIOLINK.Association,
        XSD.string,
    )

    declare_datatype_property(
        g,
        MTR.location_category,
        "location category",
        BIOLINK.Association,
        XSD.string,
    )

    declare_class(
        g,
        MTR.BiochemicalReaction,
        "Biochemical Reaction",
        "A biochemical reaction from Rhea.",
        BIOLINK.MolecularActivity,
    )

    declare_datatype_property(
        g,
        MTR.rhea_equation,
        "Rhea equation",
        MTR.BiochemicalReaction,
        XSD.string,
    )

    declare_datatype_property(
        g,
        MTR.rhea_ec_number,
        "Rhea EC number",
        MTR.BiochemicalReaction,
        XSD.string,
    )

    declare_datatype_property(
        g,
        MTR.rhea_status,
        "Rhea status",
        MTR.BiochemicalReaction,
        XSD.string,
    )

    declare_datatype_property(
        g, MTR.string_identifier, "STRING identifier", MTR.Gene, XSD.string
    )

    declare_datatype_property(
        g, MTR.string_combined_score, "STRING combined score", BIOLINK.PairwiseGeneToGeneInteraction, XSD.double
    )

    declare_datatype_property(
        g, MTR.string_network_type, "STRING network type", BIOLINK.PairwiseGeneToGeneInteraction, XSD.string
    )

    declare_datatype_property(
        g, MTR.string_experimental_score, "STRING experimental score", BIOLINK.PairwiseGeneToGeneInteraction, XSD.double
    )
    declare_datatype_property(
        g, MTR.string_database_score, "STRING database score", BIOLINK.PairwiseGeneToGeneInteraction, XSD.double
    )
    declare_datatype_property(
        g, MTR.string_textmining_score, "STRING text mining score", BIOLINK.PairwiseGeneToGeneInteraction, XSD.double
    )
    declare_datatype_property(
        g, MTR.string_coexpression_score, "STRING coexpression score", BIOLINK.PairwiseGeneToGeneInteraction, XSD.double
    )
    declare_datatype_property(
        g, MTR.string_neighborhood_score, "STRING neighborhood score", BIOLINK.PairwiseGeneToGeneInteraction, XSD.double
    )
    declare_datatype_property(
        g, MTR.string_fusion_score, "STRING fusion score", BIOLINK.PairwiseGeneToGeneInteraction, XSD.double
    )
    declare_datatype_property(
        g, MTR.string_phylogenetic_score, "STRING phylogenetic profile score", BIOLINK.PairwiseGeneToGeneInteraction,
        XSD.double
    )
    declare_class(
        g,
        MTR.VariantToRegulatoryFeatureOverlapAssociation,
        "Variant to Regulatory Feature Overlap Association",
        "An association recording that a sequence variant overlaps an Ensembl regulatory feature.",
        BIOLINK.Association,
    )

    declare_datatype_property(
        g, MTR.external_feature_id, "external feature ID", BIOLINK.RegulatoryRegion, XSD.string
    )

    declare_datatype_property(
        g, MTR.feature_description, "feature description", BIOLINK.RegulatoryRegion, XSD.string
    )

    declare_datatype_property(
        g, MTR.assembly_name, "assembly name", BIOLINK.NamedThing, XSD.string
    )

    declare_datatype_property(
        g, MTR.query_region, "query region", BIOLINK.Association, XSD.string
    )

    # if not already present elsewhere
    declare_datatype_property(g, MTR.chromosome, "chromosome", BIOLINK.GenomicEntity, XSD.string)
    declare_datatype_property(g, MTR.start_position, "start position", BIOLINK.GenomicEntity, XSD.integer)
    declare_datatype_property(g, MTR.end_position, "end position", BIOLINK.GenomicEntity, XSD.integer)

    declare_datatype_property(
        g, MTR.genomic_location, "genomic location", BIOLINK.SequenceVariant, XSD.string
    )

    # Classes
    declare_class(
        g,
        MTR.TargetSafetyLiabilityAssociation,
        "Target Safety Liability Association",
        "Association between a gene target and a disease/phenotypic feature representing a safety liability.",
        BIOLINK.GeneToPhenotypicFeatureAssociation,
    )

    declare_class(
        g,
        MTR.OpenTargetsDatasourceScore,
        "Open Targets Datasource Score",
        "A per-datasource score attached to an Open Targets association.",
    )

    # Object properties
    declare_object_property(
        g,
        MTR.has_safety_liability,
        "has safety liability",
        MTR.Gene,
        BIOLINK.DiseaseOrPhenotypicFeature,
        superprop=BIOLINK.affects,
    )

    declare_object_property(
        g,
        MTR.about_association,
        "about association",
        MTR.OpenTargetsDatasourceScore,
        BIOLINK.Association,
    )

    # Datatype properties
    declare_datatype_property(g, MTR.open_targets_target_id, "Open Targets target ID", MTR.Gene, XSD.string)
    declare_datatype_property(g, MTR.small_molecule_tractable, "small molecule tractable", MTR.Gene, XSD.boolean)
    declare_datatype_property(g, MTR.antibody_tractable, "antibody tractable", MTR.Gene, XSD.boolean)
    declare_datatype_property(g, MTR.protac_tractable, "PROTAC tractable", MTR.Gene, XSD.boolean)
    declare_datatype_property(g, MTR.other_clinical_tractable, "other clinical tractable", MTR.Gene, XSD.boolean)

    declare_datatype_property(g, MTR.external_event_id, "external event ID", BIOLINK.Association, XSD.string)
    declare_datatype_property(g, MTR.tissue_context, "tissue context", BIOLINK.Association, XSD.string)

    declare_datatype_property(g, MTR.max_clinical_phase, "max clinical phase", BIOLINK.NamedThing, XSD.integer)

    declare_datatype_property(g, MTR.ot_genetics_score, "Open Targets genetics score", BIOLINK.Association, XSD.double)
    declare_datatype_property(g, MTR.datasource_id, "datasource ID", MTR.OpenTargetsDatasourceScore, XSD.string)
    declare_datatype_property(g, MTR.datasource_score, "datasource score", MTR.OpenTargetsDatasourceScore, XSD.double)
    return g