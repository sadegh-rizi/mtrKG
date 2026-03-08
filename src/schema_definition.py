# core/namespaces.py
from rdflib import Namespace
import requests
import time
import urllib.parse
from rdflib import Graph, Literal, RDF, URIRef, Namespace
from rdflib.namespace import RDFS, XSD, OWL, SKOS
import re


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
UP = Namespace("http://purl.uniprot.org/core/")
OBAN = Namespace("http://purl.org/oban/")
SO = Namespace("http://purl.obolibrary.org/obo/SO_")
HGNC = Namespace("http://identifiers.org/hgnc.symbol/")



def build_schema(g):
    """Utility to bind all prefixes to the graph at once."""
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

    ontology_uri = URIRef("https://metabolite-ratio-app.unil.ch/ontology")
    g.add((ontology_uri, RDF.type, OWL.Ontology))
    g.add((ontology_uri, RDFS.label, Literal("rQTL Metabolic Ratio Ontology")))
    g.add((ontology_uri, RDFS.comment, Literal("A schema defining genetically regulated metabolic reactions.")))

    # ==========================================
    # 3. DEFINE CLASSES (The "Nouns")
    # ==========================================
    classes = {
        "MetaboliteRatio": "A ratio between two metabolite concentrations.",
        "Metabolite": "A small chemical compound.",
        "SNP": "Single Nucleotide Polymorphism.",
        "Gene": "A region of DNA.",
        "Phenotype": "An observable clinical trait or disease."
    }

    for cls_name, description in classes.items():
        node = MTR[cls_name]
        g.add((node, RDF.type, OWL.Class))
        g.add((node, RDFS.label, Literal(cls_name)))
        g.add((node, RDFS.comment, Literal(description)))

    # Integrate UniProt Enzyme class (as discussed previously)
    g.add((MTR.Gene, RDFS.subClassOf, UP.Protein))  # A gene maps to a protein/enzyme

    # ==========================================
    # 4. DEFINE OBJECT PROPERTIES (The "Edges" between nodes)
    # Here we define the Domain (starts at) and Range (ends at)
    # ==========================================

    # Ratio -> Metabolite
    g.add((MTR.has_numerator, RDF.type, OWL.ObjectProperty))
    g.add((MTR.has_numerator, RDFS.domain, MTR.MetaboliteRatio))
    g.add((MTR.has_numerator, RDFS.range, MTR.Metabolite))

    g.add((MTR.has_denominator, RDF.type, OWL.ObjectProperty))
    g.add((MTR.has_denominator, RDFS.domain, MTR.MetaboliteRatio))
    g.add((MTR.has_denominator, RDFS.range, MTR.Metabolite))

    # SNP -> Gene
    g.add((MTR.closest_gene, RDF.type, OWL.ObjectProperty))
    g.add((MTR.closest_gene, RDFS.domain, MTR.SNP))
    g.add((MTR.closest_gene, RDFS.range, MTR.Gene))

    g.add((MTR.is_eQTL_for, RDF.type, OWL.ObjectProperty))
    g.add((MTR.is_eQTL_for, RDFS.domain, MTR.SNP))
    g.add((MTR.is_eQTL_for, RDFS.range, MTR.Gene))

    # SNP -> Phenotype
    g.add((MTR.associated_with_trait, RDF.type, OWL.ObjectProperty))
    g.add((MTR.associated_with_trait, RDFS.domain, MTR.SNP))
    g.add((MTR.associated_with_trait, RDFS.range, MTR.Phenotype))

    # ==========================================
    # 5. DEFINE DATATYPE PROPERTIES (The "Edges" to values/numbers)
    # ==========================================

    # Association -> Float
    g.add((MTR.p_value, RDF.type, OWL.DatatypeProperty))
    g.add((MTR.p_value, RDFS.domain, OBAN.association))
    g.add((MTR.p_value, RDFS.range, XSD.float))

    g.add((MTR.beta, RDF.type, OWL.DatatypeProperty))
    g.add((MTR.beta, RDFS.domain, OBAN.association))
    g.add((MTR.beta, RDFS.range, XSD.float))

    g.add((MTR.pGain, RDF.type, OWL.DatatypeProperty))
    g.add((MTR.pGain, RDFS.domain, OBAN.association))
    g.add((MTR.pGain, RDFS.range, XSD.float))

    g.add((MTR.SNP, RDFS.subClassOf, SO["0000694"])) # SO:0000694 is 'SNP'

    # 2. Hybrid Ontology Approach: Map custom edges to Biolink
    g.add((MTR.is_eQTL_for, RDFS.subPropertyOf, BIOLINK.affects_expression_of))
    g.add((MTR.associated_with_trait, RDFS.subPropertyOf, BIOLINK.associated_with))

    # 3. Hybrid Ontology Approach: Map statistics to SIO
    g.add((MTR.p_value, RDFS.subPropertyOf, SIO["000765"])) # SIO:000765 is 'p-value'
    g.add((MTR.beta, RDFS.subPropertyOf, SIO["001078"])) # SIO:001078 is 'beta' (or similar estimate)

    # 4. Add PROV-O metadata to your generated dataset
    g.add((MTR.Dataset, RDF.type, PROV.Entity))
    g.add((MTR.Dataset, PROV.wasGeneratedBy, Literal("Python JSON-to-RDF parser")))

    # --- HYBRID ALIGNMENTS ---

    g.add((OBAN.association, OWL.equivalentClass, BIOLINK.Association))
    g.add((OBAN.has_subject, OWL.equivalentProperty, BIOLINK.has_subject))
    g.add((OBAN.has_object, OWL.equivalentProperty, BIOLINK.has_object))

    # ==========================================
    # STRING DB: PROTEIN-PROTEIN INTERACTIONS
    # ==========================================

    # 1. Define the Interaction Category Class
    g.add((BIOLINK.ProteinProteinInteraction, RDF.type, OWL.Class))
    g.add((BIOLINK.ProteinProteinInteraction, RDFS.subClassOf, BIOLINK.Association))
    g.add((BIOLINK.ProteinProteinInteraction, RDFS.label, Literal("Protein-Protein Interaction")))

    # 2. Define the Direct Interaction Property
    g.add((BIOLINK.interacts_with, RDF.type, OWL.ObjectProperty))
    g.add((BIOLINK.interacts_with, RDFS.label, Literal("interacts with")))
    g.add((BIOLINK.interacts_with, RDFS.domain, MTR.Gene))  # From a Gene/Protein
    g.add((BIOLINK.interacts_with, RDFS.range, MTR.Gene))  # To another Gene/Protein

    # 3. Define the Interaction Score Property (Datatype)
    g.add((MTR.interaction_score, RDF.type, OWL.DatatypeProperty))
    g.add((MTR.interaction_score, RDFS.label, Literal("STRING interaction score")))
    # Assuming you applied the ontology alignment, binding to BIOLINK.Association covers OBAN too
    g.add((MTR.interaction_score, RDFS.domain, BIOLINK.Association))
    g.add((MTR.interaction_score, RDFS.range, XSD.float))

    # ==========================================
    # EPIGENETICS (EWAS): CLASSES & PROPERTIES
    # ==========================================

    # 1. Define the Genomic Entity Class (for CpG Methylation Sites)
    g.add((BIOLINK.GenomicEntity, RDF.type, OWL.Class))
    g.add((BIOLINK.GenomicEntity, RDFS.label, Literal("Genomic Entity")))
    g.add((BIOLINK.GenomicEntity, RDFS.comment,
           Literal("A structural or functional feature of a genome, such as a CpG site.")))

    # 2. Define the Phenotypic Feature Class (for Environmental Traits / Biomarkers)
    # Note: You might already have biolink:Disease, this broadens it to non-disease traits like BMI or smoking
    g.add((BIOLINK.PhenotypicFeature, RDF.type, OWL.Class))
    g.add((BIOLINK.PhenotypicFeature, RDFS.label, Literal("Phenotypic Feature")))

    # 3. Define the 'regulates' Object Property
    g.add((BIOLINK.regulates, RDF.type, OWL.ObjectProperty))
    g.add((BIOLINK.regulates, RDFS.label, Literal("regulates")))

    # Strict typing: Only Genomic Entities (CpGs) can regulate Genes in this context
    g.add((BIOLINK.regulates, RDFS.domain, BIOLINK.GenomicEntity))
    g.add((BIOLINK.regulates, RDFS.range, MTR.Gene))

    # Define the overlaps property for Genomic Entities
    g.add((BIOLINK.overlaps, RDF.type, OWL.ObjectProperty))
    g.add((BIOLINK.overlaps, RDFS.label, Literal("overlaps with")))
    g.add((BIOLINK.overlaps, RDFS.domain, MTR.SNP))
    g.add((BIOLINK.overlaps, RDFS.range, BIOLINK.GenomicEntity))

    # ==========================================
    # PATHWAYS (REACTOME): CLASSES & PROPERTIES
    # ==========================================

    # 1. Define the Pathway Class
    g.add((BIOLINK.Pathway, RDF.type, OWL.Class))
    g.add((BIOLINK.Pathway, RDFS.label, Literal("Biological Pathway")))

    # 2. Define the 'participates_in' Object Property
    g.add((BIOLINK.participates_in, RDF.type, OWL.ObjectProperty))
    g.add((BIOLINK.participates_in, RDFS.label, Literal("participates in")))

    # Strict typing: Only Genes/Proteins can participate in Pathways
    g.add((BIOLINK.participates_in, RDFS.domain, MTR.Gene))
    g.add((BIOLINK.participates_in, RDFS.range, BIOLINK.Pathway))

    # ==========================================
    # PATHWAYS (REACTOME): CLASSES & PROPERTIES
    # ==========================================

    # 1. Define a Parent Class for both Genes and Metabolites
    g.add((BIOLINK.MolecularEntity, RDF.type, OWL.Class))
    g.add((MTR.Gene, RDFS.subClassOf, BIOLINK.MolecularEntity))
    g.add((MTR.Metabolite, RDFS.subClassOf, BIOLINK.MolecularEntity))

    # 2. Define the Pathway Class
    g.add((BIOLINK.Pathway, RDF.type, OWL.Class))
    g.add((BIOLINK.Pathway, RDFS.label, Literal("Biological Pathway")))

    # 3. Define the 'participates_in' Object Property
    g.add((BIOLINK.participates_in, RDF.type, OWL.ObjectProperty))
    g.add((BIOLINK.participates_in, RDFS.label, Literal("participates in")))

    # Update Domain: Now BOTH Genes and Metabolites can participate in pathways!
    g.add((BIOLINK.participates_in, RDFS.domain, BIOLINK.MolecularEntity))
    g.add((BIOLINK.participates_in, RDFS.range, BIOLINK.Pathway))

    # 1. Define the Classes
    g.add((BIOLINK.CellularComponent, RDF.type, OWL.Class))
    g.add((BIOLINK.CellularComponent, RDFS.label, Literal("Cellular Component")))

    g.add((BIOLINK.GrossAnatomicalStructure, RDF.type, OWL.Class))
    g.add((BIOLINK.GrossAnatomicalStructure, RDFS.label, Literal("Tissue or Biofluid")))

    # 2. Define the 'located_in' Object Property
    g.add((BIOLINK.located_in, RDF.type, OWL.ObjectProperty))
    g.add((BIOLINK.located_in, RDFS.label, Literal("located in")))

    # Domain is the Molecular Entity (Metabolite or Protein)
    g.add((BIOLINK.located_in, RDFS.domain, BIOLINK.MolecularEntity))
    # We won't strictly enforce a single range here, as it could be a Cell Component OR Tissue


    return g