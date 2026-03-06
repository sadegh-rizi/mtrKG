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


    return g