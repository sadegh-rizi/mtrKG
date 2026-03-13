import json
import logging
import re
import time
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

from rdflib import Graph, Literal, RDF
from rdflib.namespace import RDFS, XSD, SKOS

from schema_definition import (
    MTR,
    BIOLINK,
    CHEBI,
    HMDB,
    DBSNP,
    PROV,
    build_schema,
)

OUTPUT_DIR = Path("../output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(OUTPUT_DIR / "rqtl_integration.log", mode="w"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

MISSING_VALUES = {None, "", "NA", "NaN", "nan", "null", "None"}
CHEBI_TOKEN = re.compile(r"CHEBI[:_](\d+)", re.IGNORECASE)

CAUSAL_DATASETS = ["eQTL", "pQTL", "GTEx", "finngen", "mrlink2_traits", "yang", "rosmap", "western"]
TARGET_OUTCOMES = ["ratio", "numerator_metabolite", "denominator_metabolite"]

# DIRECTION RULES
LOCAL_TO_TRAIT_DATASETS = {"finngen", "mrlink2_traits"}    # ratio/metabolite --> trait
TRANSCRIPT_TO_LOCAL_DATASETS = {"eQTL", "GTEx"}            # transcript --> ratio/metabolite
PROTEIN_TO_LOCAL_DATASETS = {"pQTL", "western", "rosmap"}  # protein --> ratio/metabolite
OTHER_EXTERNAL_TO_LOCAL = {"yang"}


def clean_uri_string(text: str) -> str:
    return urllib.parse.quote(str(text).strip().replace(" ", "_").replace('"', ""))


def is_missing(value) -> bool:
    return value in MISSING_VALUES


def coerce_value(value, datatype):
    if is_missing(value):
        return None

    if datatype == XSD.boolean:
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            v = value.strip().lower()
            if v in {"true", "1", "yes"}:
                return True
            if v in {"false", "0", "no"}:
                return False
        if isinstance(value, (int, float)):
            return bool(value)
        return None

    if datatype == XSD.integer:
        return int(float(value))

    if datatype in {XSD.float, XSD.double}:
        return float(value)

    return str(value).strip()


def add_safe_literal(graph: Graph, subject, predicate, value, datatype=None) -> None:
    if is_missing(value):
        return
    try:
        if datatype:
            coerced = coerce_value(value, datatype)
            if coerced is None:
                return
            graph.add((subject, predicate, Literal(coerced, datatype=datatype)))
        else:
            graph.add((subject, predicate, Literal(str(value).strip())))
    except (ValueError, TypeError):
        logger.warning("Skipping invalid value for %s %s %r", subject, predicate, value)


def set_safe_literal(graph: Graph, subject, predicate, value, datatype=None) -> None:
    """
    For properties intended to be single-valued.
    """
    if is_missing(value):
        return
    try:
        if datatype:
            coerced = coerce_value(value, datatype)
            if coerced is None:
                return
            graph.set((subject, predicate, Literal(coerced, datatype=datatype)))
        else:
            graph.set((subject, predicate, Literal(str(value).strip())))
    except (ValueError, TypeError):
        logger.warning("Skipping invalid value for %s %s %r", subject, predicate, value)


def normalize_variant_position(raw_position, raw_chromosome=None):
    """
    Normalizes position fields from values like:
      151921111
      "151921111"
      "1:151921111"
      "chr1:151921111"
      "1:151921111-151921111"

    Returns:
      chromosome_str, position_int, genomic_location_str
    """
    chromosome = None if is_missing(raw_chromosome) else str(raw_chromosome).strip()
    position_int = None
    genomic_location = None

    if is_missing(raw_position):
        return chromosome, position_int, genomic_location

    s = str(raw_position).strip()
    genomic_location = s

    if ":" in s:
        chrom_part, pos_part = s.split(":", 1)
        chrom_part = chrom_part.replace("chr", "").strip()

        if not chromosome:
            chromosome = chrom_part

        pos_part = pos_part.split("-")[0].replace(",", "").strip()

        try:
            position_int = int(pos_part)
        except ValueError:
            position_int = None

        return chromosome, position_int, genomic_location

    try:
        position_int = int(float(s))
    except ValueError:
        position_int = None

    if chromosome and position_int is not None:
        genomic_location = f"{chromosome}:{position_int}"

    return chromosome, position_int, genomic_location


def normalize_chebi_id(raw_value):
    if is_missing(raw_value):
        return None
    # Handle native numeric values first (common in JSON: 16347.0).
    if isinstance(raw_value, (int, float)):
        try:
            chebi_num = int(float(raw_value))
            return str(chebi_num) if chebi_num > 0 else None
        except (ValueError, TypeError, OverflowError):
            return None

    raw = str(raw_value).strip()
    if is_missing(raw):
        return None

    # Explicit CURIE/ID formats (CHEBI:16347 or CHEBI_16347).
    token_match = CHEBI_TOKEN.search(raw)
    if token_match:
        chebi_num = int(token_match.group(1))
        return str(chebi_num) if chebi_num > 0 else None

    # Plain numeric strings, including decimal representations like "16347.0".
    try:
        chebi_num = int(float(raw))
        return str(chebi_num) if chebi_num > 0 else None
    except (ValueError, TypeError, OverflowError):
        return None


def mint_ratio_node(ratio_accession: str):
    return MTR[f"ratio/{clean_uri_string(ratio_accession)}"]


def mint_gene_node(symbol: str):
    return MTR[f"gene/{clean_uri_string(symbol)}"]


def mint_transcript_node(symbol: str):
    return MTR[f"transcript/{clean_uri_string(symbol)}"]


def mint_protein_node(symbol: str):
    return MTR[f"protein/{clean_uri_string(symbol)}"]


def mint_trait_node(name: str):
    return MTR[f"trait/{clean_uri_string(name)}"]


def mint_assoc_node(rsid: str, ratio_accession: str):
    return MTR[f"association/variant_ratio/{clean_uri_string(rsid)}__{clean_uri_string(ratio_accession)}"]


def mint_causal_node(dataset_name: str, exposure_name: str, target_level: str, ratio_accession: str, rsid: str):
    return MTR[
        f"causal/{clean_uri_string(dataset_name)}__{clean_uri_string(exposure_name)}__"
        f"{clean_uri_string(target_level)}__{clean_uri_string(ratio_accession)}__{clean_uri_string(rsid)}"
    ]


def mint_metabolite_node(metab: dict):
    chebi_id = normalize_chebi_id(metab.get("chebi"))
    if chebi_id:
        return CHEBI[chebi_id]

    if not is_missing(metab.get("accession")):
        return MTR[f"metabolite/accession/{clean_uri_string(metab['accession'])}"]

    if not is_missing(metab.get("metabolon")):
        return MTR[f"metabolite/metabolon/{clean_uri_string(metab['metabolon'])}"]

    name = metab.get("name", "unknown_metabolite")
    return MTR[f"metabolite/name/{clean_uri_string(name)}"]


def summarize_entry(entry):
    regions = entry.get("associated_regions", {})
    region_count = len(regions) if isinstance(regions, dict) else 0

    missense_count = 0
    cnv_count = 0
    causal_test_count = 0

    if isinstance(regions, dict):
        for _, details in regions.items():
            if not isinstance(details, dict):
                continue

            missense_count += len(details.get("missense_variants_in_ld", []) or [])
            cnv_count += len(details.get("cnv_chiara", []) or [])

            for dataset_name in CAUSAL_DATASETS:
                dataset_data = details.get(dataset_name, {})
                if not isinstance(dataset_data, dict):
                    continue

                for target_level in TARGET_OUTCOMES:
                    tests = dataset_data.get(target_level, [])
                    if isinstance(tests, list):
                        causal_test_count += len(tests)

    return {
        "regions": region_count,
        "missense": missense_count,
        "cnv": cnv_count,
        "causal_tests": causal_test_count,
    }


def ensure_typed_labeled_cached(graph: Graph, cache: set, node, class_uri, label=None):
    """
    Avoid re-adding the same type/label over and over.
    """
    cache_key = (node, class_uri)
    if cache_key not in cache:
        graph.add((node, RDF.type, class_uri))
        cache.add(cache_key)

    if label and not is_missing(label):
        graph.set((node, RDFS.label, Literal(str(label).strip(), datatype=XSD.string)))


def add_dataset_provenance(graph: Graph, json_file_path: str):
    dataset_node = MTR["dataset/rqtl_integrated"]
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_node = MTR[f"activity/integration_run/{run_id}"]
    source_node = MTR[f"source/{clean_uri_string(Path(json_file_path).name)}"]

    graph.add((dataset_node, RDF.type, MTR.DatasetRelease))
    graph.set((dataset_node, RDFS.label, Literal("Integrated rQTL dataset", datatype=XSD.string)))

    graph.add((run_node, RDF.type, MTR.IntegrationRun))
    graph.set((run_node, RDFS.label, Literal(f"Integration run {run_id}", datatype=XSD.string)))

    graph.add((source_node, RDF.type, PROV.Entity))
    graph.set((source_node, RDFS.label, Literal(str(json_file_path), datatype=XSD.string)))

    graph.add((dataset_node, PROV.wasGeneratedBy, run_node))
    graph.add((run_node, PROV.used, source_node))
    graph.add((run_node, PROV.generated, dataset_node))


def extract_llm_response(entry: dict) -> dict:
    llm_raw = entry.get("llm_response")
    if not isinstance(llm_raw, dict):
        return {}

    candidate = llm_raw.get("reponse", llm_raw.get("response", {}))
    return candidate if isinstance(candidate, dict) else {}


def add_metabolite_metadata(g: Graph, metab_node, metab: dict):
    set_safe_literal(g, metab_node, MTR.hmdb_id, metab.get("hmdb"), XSD.string)
    set_safe_literal(g, metab_node, MTR.kegg_id, metab.get("kegg"), XSD.string)
    set_safe_literal(g, metab_node, MTR.local_accession, metab.get("accession"), XSD.string)
    set_safe_literal(g, metab_node, MTR.metabolon_id, metab.get("metabolon"), XSD.string)
    set_safe_literal(g, metab_node, MTR.inchikey, metab.get("inchikey"), XSD.string)

    if not is_missing(metab.get("hmdb")):
        g.add((metab_node, SKOS.exactMatch, HMDB[str(metab["hmdb"]).strip()]))


def resolve_local_outcome_node(target_level: str, ratio_node, metabolite_nodes: dict):
    if target_level == "ratio":
        return ratio_node
    return metabolite_nodes.get(target_level)


def create_exposure_outcome_nodes(g: Graph, dataset_name: str, target_level: str, exposure_name: str,
                                  ratio_node, metabolite_nodes: dict, node_cache: set):
    local_node = resolve_local_outcome_node(target_level, ratio_node, metabolite_nodes)
    if local_node is None:
        return None, None, None

    if dataset_name in LOCAL_TO_TRAIT_DATASETS:
        trait_node = mint_trait_node(exposure_name)
        ensure_typed_labeled_cached(g, node_cache, trait_node, MTR.Phenotype, exposure_name)
        return local_node, trait_node, "local_to_trait"

    if dataset_name in TRANSCRIPT_TO_LOCAL_DATASETS:
        gene_node = mint_gene_node(exposure_name)
        transcript_node = mint_transcript_node(exposure_name)

        ensure_typed_labeled_cached(g, node_cache, gene_node, MTR.Gene, exposure_name)
        ensure_typed_labeled_cached(g, node_cache, transcript_node, MTR.Transcript, f"{exposure_name} RNA")
        g.add((gene_node, BIOLINK.transcribed_to, transcript_node))

        return transcript_node, local_node, "external_to_local"

    if dataset_name in PROTEIN_TO_LOCAL_DATASETS:
        gene_node = mint_gene_node(exposure_name)
        protein_node = mint_protein_node(exposure_name)

        ensure_typed_labeled_cached(g, node_cache, gene_node, MTR.Gene, exposure_name)
        ensure_typed_labeled_cached(g, node_cache, protein_node, MTR.Protein, exposure_name)
        g.add((gene_node, BIOLINK.has_gene_product, protein_node))

        return protein_node, local_node, "external_to_local"

    trait_node = mint_trait_node(exposure_name)
    ensure_typed_labeled_cached(g, node_cache, trait_node, MTR.Phenotype, exposure_name)
    return trait_node, local_node, "external_to_local"


def add_rqtl_to_graph(
    json_file_path: str,
    g: Graph = None,
    checkpoint_every: int = 250,
    log_graph_every: int = 50,
    region_log_every: int = 25,
) -> Graph:
    if g is None:
        g = Graph()

    build_schema(g)
    add_dataset_provenance(g, json_file_path)

    logger.info("--- Starting rQTL JSON Integration (%s) ---", json_file_path)

    try:
        with open(json_file_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        logger.info("Loaded JSON with %d entries.", len(data))
    except Exception as exc:
        logger.exception("Failed to load JSON file: %s", exc)
        return g

    # caches to avoid redundant graph churn
    node_cache = set()

    for i, entry in enumerate(data):
        ratio_accession = entry.get("ratio_accession")
        if is_missing(ratio_accession):
            logger.warning("Skipping entry %d without ratio_accession", i)
            continue

        ratio_node = mint_ratio_node(str(ratio_accession))
        ratio_name = entry.get("ratio_name", ratio_accession)

        entry_start = time.perf_counter()
        summary = summarize_entry(entry)

        logger.info(
            "[%d/%d] Processing %s | regions=%d missense=%d cnv=%d causal_tests=%d",
            i + 1,
            len(data),
            ratio_accession,
            summary["regions"],
            summary["missense"],
            summary["cnv"],
            summary["causal_tests"],
        )

        if summary["regions"] > 100 or summary["causal_tests"] > 5000:
            logger.warning(
                "Large entry detected: %s | regions=%d causal_tests=%d",
                ratio_accession,
                summary["regions"],
                summary["causal_tests"],
            )

        # =========================
        # 1. Ratio node
        # =========================
        ensure_typed_labeled_cached(g, node_cache, ratio_node, MTR.MetaboliteRatio, ratio_name)
        set_safe_literal(g, ratio_node, MTR.ratio_formula, ratio_name, XSD.string)
        set_safe_literal(g, ratio_node, MTR.reaction_distance, entry.get("reaction_distance"), XSD.integer)
        set_safe_literal(g, ratio_node, MTR.max_pgain, entry.get("max_pgain"), XSD.double)

        llm_resp = extract_llm_response(entry)
        set_safe_literal(g, ratio_node, MTR.llm_ratio_explanation, llm_resp.get("ratio_explanation"), XSD.string)
        set_safe_literal(g, ratio_node, MTR.llm_ratio_evidence, llm_resp.get("ratio_evidence"), XSD.string)
        set_safe_literal(g, ratio_node, MTR.llm_phenotype_driver, llm_resp.get("phenotype_driver"), XSD.string)
        set_safe_literal(
            g,
            ratio_node,
            MTR.llm_gene_phenotype_relationship,
            llm_resp.get("gene_ratio_phenotype_relationship"),
            XSD.string
        )

        # =========================
        # 2. Metabolites
        # =========================
        metabolite_nodes = {}

        for side, relation in [
            ("numerator_metabolite", MTR.has_numerator),
            ("denominator_metabolite", MTR.has_denominator),
        ]:
            metab = entry.get(side, {})
            if not isinstance(metab, dict):
                continue

            metab_node = mint_metabolite_node(metab)
            metabolite_nodes[side] = metab_node

            ensure_typed_labeled_cached(g, node_cache, metab_node, MTR.Metabolite, metab.get("name", side))
            add_metabolite_metadata(g, metab_node, metab)

            if side == "numerator_metabolite":
                set_safe_literal(
                    g, metab_node, MTR.llm_explanation,
                    llm_resp.get("numerator_metabolite_explanation"), XSD.string
                )
            else:
                set_safe_literal(
                    g, metab_node, MTR.llm_explanation,
                    llm_resp.get("denominator_metabolite_explanation"), XSD.string
                )

            g.add((ratio_node, relation, metab_node))

        # =========================
        # 3. Associated regions
        # =========================
        regions = entry.get("associated_regions", {})
        if not isinstance(regions, dict):
            continue

        for r_idx, (rsid, details) in enumerate(regions.items(), start=1):
            if is_missing(rsid):
                continue

            if r_idx % region_log_every == 0:
                logger.info(
                    "   ratio=%s | region %d/%d | rsid=%s",
                    ratio_accession,
                    r_idx,
                    len(regions),
                    rsid,
                )

            snp_node = DBSNP[str(rsid).strip()]
            ensure_typed_labeled_cached(g, node_cache, snp_node, MTR.SNP, str(rsid).strip())

            stats = details.get("summary_statistics", {}) if isinstance(details, dict) else {}

            set_safe_literal(g, snp_node, MTR.pos_name, stats.get("pos_name"), XSD.string)

            norm_chr, norm_pos, norm_location = normalize_variant_position(
                raw_position=stats.get("position"),
                raw_chromosome=stats.get("chromosome"),
            )

            set_safe_literal(g, snp_node, MTR.chromosome, norm_chr, XSD.string)
            set_safe_literal(g, snp_node, MTR.position, norm_pos, XSD.integer)
            set_safe_literal(g, snp_node, MTR.genomic_location, norm_location, XSD.string)

            assoc_node = mint_assoc_node(str(rsid), str(ratio_accession))
            ensure_typed_labeled_cached(
                g,
                node_cache,
                assoc_node,
                MTR.VariantToMetaboliteRatioAssociation,
                f"{rsid} associated with {ratio_name}"
            )

            g.add((assoc_node, BIOLINK.subject, snp_node))
            g.add((assoc_node, BIOLINK.object, ratio_node))
            g.add((assoc_node, BIOLINK.predicate, MTR.associated_with_ratio))

            set_safe_literal(g, assoc_node, MTR.effect_allele, stats.get("effect_allele"), XSD.string)
            set_safe_literal(g, assoc_node, MTR.reference_allele, stats.get("reference_allele"), XSD.string)
            set_safe_literal(g, assoc_node, MTR.maf, stats.get("maf"), XSD.double)
            set_safe_literal(g, assoc_node, MTR.beta, stats.get("beta"), XSD.double)
            set_safe_literal(g, assoc_node, MTR.standard_error, stats.get("se"), XSD.double)
            set_safe_literal(g, assoc_node, MTR.z_score, stats.get("z"), XSD.double)
            set_safe_literal(g, assoc_node, MTR.log_p_value, stats.get("log_pval"), XSD.double)

            set_safe_literal(g, assoc_node, MTR.beta_1, stats.get("beta_1"), XSD.double)
            set_safe_literal(g, assoc_node, MTR.se_1, stats.get("se_1"), XSD.double)
            set_safe_literal(g, assoc_node, MTR.z_1, stats.get("z_1"), XSD.double)
            set_safe_literal(g, assoc_node, MTR.beta_2, stats.get("beta_2"), XSD.double)
            set_safe_literal(g, assoc_node, MTR.se_2, stats.get("se_2"), XSD.double)
            set_safe_literal(g, assoc_node, MTR.z_2, stats.get("z_2"), XSD.double)

            set_safe_literal(g, assoc_node, MTR.numerator_driven, stats.get("numerator_driven"), XSD.boolean)
            set_safe_literal(g, assoc_node, MTR.denominator_driven, stats.get("denominator_driven"), XSD.boolean)
            set_safe_literal(g, assoc_node, MTR.already_found, stats.get("already_found"), XSD.boolean)
            set_safe_literal(g, assoc_node, MTR.pgain, stats.get("pgain"), XSD.double)
            set_safe_literal(g, assoc_node, MTR.log_pgain, stats.get("log_pgain"), XSD.double)
            set_safe_literal(g, assoc_node, MTR.infomap_cluster, stats.get("infomap_cluster"), XSD.integer)

            for gene_symbol in stats.get("closest_genes", []):
                if is_missing(gene_symbol):
                    continue
                gene_node = mint_gene_node(str(gene_symbol))
                ensure_typed_labeled_cached(g, node_cache, gene_node, MTR.Gene, str(gene_symbol))
                g.add((snp_node, MTR.closest_gene, gene_node))

            # =========================
            # 4. Missense variants in LD
            # =========================
            for missense in details.get("missense_variants_in_ld", []):
                if not isinstance(missense, dict):
                    continue

                ms_rsid = missense.get("missense_variant")
                if is_missing(ms_rsid):
                    continue

                ms_node = DBSNP[str(ms_rsid).strip()]
                ensure_typed_labeled_cached(g, node_cache, ms_node, MTR.MissenseVariant, str(ms_rsid).strip())
                g.add((snp_node, MTR.in_ld_with, ms_node))

                set_safe_literal(g, ms_node, MTR.ld_score, missense.get("ld"), XSD.double)
                set_safe_literal(g, ms_node, MTR.variant_location, missense.get("Location"), XSD.string)
                set_safe_literal(g, ms_node, MTR.allele_change, missense.get("Allele"), XSD.string)
                set_safe_literal(g, ms_node, MTR.consequence, missense.get("Consequence"), XSD.string)
                set_safe_literal(g, ms_node, MTR.codons, missense.get("Codons"), XSD.string)

                gene_name = missense.get("gene_name", missense.get("Gene"))
                if not is_missing(gene_name):
                    ms_gene_node = mint_gene_node(str(gene_name))
                    ensure_typed_labeled_cached(g, node_cache, ms_gene_node, MTR.Gene, str(gene_name))
                    g.add((ms_node, MTR.affects_gene, ms_gene_node))

            # =========================
            # 5. CNVs
            # =========================
            for cnv in details.get("cnv_chiara", []):
                if not isinstance(cnv, dict):
                    continue

                cnv_id = cnv.get("ID")
                if is_missing(cnv_id):
                    continue

                cnv_node = MTR[f"cnv/{clean_uri_string(cnv_id)}"]
                ensure_typed_labeled_cached(g, node_cache, cnv_node, MTR.CopyNumberVariant, str(cnv_id))
                g.add((snp_node, MTR.co_located_with_cnv, cnv_node))

                set_safe_literal(g, cnv_node, MTR.chromosome, cnv.get("CHR"), XSD.string)
                set_safe_literal(g, cnv_node, MTR.position, cnv.get("POS"), XSD.integer)
                set_safe_literal(g, cnv_node, MTR.num_cnv, cnv.get("NumCNV"), XSD.integer)
                set_safe_literal(g, cnv_node, MTR.num_dup, cnv.get("NumDup"), XSD.integer)
                set_safe_literal(g, cnv_node, MTR.num_del, cnv.get("NumDel"), XSD.integer)
                set_safe_literal(g, cnv_node, MTR.num_neutral, cnv.get("NumNeutral"), XSD.integer)
                set_safe_literal(g, cnv_node, MTR.freq_cnv, cnv.get("FreqCNV"), XSD.double)
                set_safe_literal(g, cnv_node, MTR.freq_dup, cnv.get("FreqDup"), XSD.double)
                set_safe_literal(g, cnv_node, MTR.freq_del, cnv.get("FreqDel"), XSD.double)

            # =========================
            # 6. MR / coloc
            # =========================
            for dataset_name in CAUSAL_DATASETS:
                dataset_data = details.get(dataset_name, {})
                if not isinstance(dataset_data, dict):
                    continue

                dataset_total = 0
                for target_level in TARGET_OUTCOMES:
                    tests = dataset_data.get(target_level, [])
                    if isinstance(tests, list):
                        dataset_total += len(tests)

                if dataset_total > 0:
                    logger.info(
                        "      rsid=%s | dataset=%s | tests=%d",
                        rsid,
                        dataset_name,
                        dataset_total,
                    )

                for target_level in TARGET_OUTCOMES:
                    causal_tests = dataset_data.get(target_level, [])
                    if not isinstance(causal_tests, list):
                        continue

                    for test in causal_tests:
                        if not isinstance(test, dict):
                            continue

                        exposure_name = test.get("exposure")
                        if is_missing(exposure_name):
                            continue

                        exposure_node, outcome_node, direction = create_exposure_outcome_nodes(
                            g=g,
                            dataset_name=dataset_name,
                            target_level=target_level,
                            exposure_name=str(exposure_name),
                            ratio_node=ratio_node,
                            metabolite_nodes=metabolite_nodes,
                            node_cache=node_cache,
                        )

                        if exposure_node is None or outcome_node is None:
                            continue

                        assessment_node = mint_causal_node(
                            dataset_name=dataset_name,
                            exposure_name=str(exposure_name),
                            target_level=target_level,
                            ratio_accession=str(ratio_accession),
                            rsid=str(rsid),
                        )

                        ensure_typed_labeled_cached(
                            g,
                            node_cache,
                            assessment_node,
                            MTR.CausalAssessment,
                            f"{dataset_name} causal assessment: {exposure_name} ({direction})"
                        )

                        g.add((assessment_node, MTR.has_exposure, exposure_node))
                        g.add((assessment_node, MTR.has_outcome, outcome_node))
                        g.add((assessment_node, MTR.associated_variant, snp_node))

                        g.add((assessment_node, BIOLINK.subject, exposure_node))
                        g.add((assessment_node, BIOLINK.object, outcome_node))
                        g.add((assessment_node, BIOLINK.predicate, MTR.causal_influence_on))

                        # single-valued metadata
                        set_safe_literal(g, assessment_node, MTR.dataset_source, dataset_name, XSD.string)
                        set_safe_literal(g, assessment_node, MTR.target_level, target_level, XSD.string)
                        set_safe_literal(g, assessment_node, MTR.causal_direction, direction, XSD.string)

                        set_safe_literal(g, assessment_node, MTR.exposure_file, test.get("exposure_file"), XSD.string)
                        set_safe_literal(g, assessment_node, MTR.outcome_file, test.get("outcome_file"), XSD.string)
                        set_safe_literal(g, assessment_node, MTR.region, test.get("region"), XSD.string)
                        set_safe_literal(g, assessment_node, MTR.measured_in_tissue, test.get("tissue"), XSD.string)

                        set_safe_literal(g, assessment_node, MTR.variance_explained, test.get("var_explained"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.snps_overlap, test.get("m_snps_overlap"), XSD.integer)

                        set_safe_literal(g, assessment_node, MTR.alpha, test.get("alpha"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.se_alpha, test.get("se(alpha)"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.p_alpha, test.get("p(alpha)"), XSD.double)

                        set_safe_literal(g, assessment_node, MTR.sigma_y, test.get("sigma_y"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.se_sigma_y, test.get("se(sigma_y)"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.p_sigma_y, test.get("p(sigma_y)"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.sigma_x, test.get("sigma_x"), XSD.double)

                        set_safe_literal(g, assessment_node, MTR.beta_ivw, test.get("beta_ivw"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.se_ivw, test.get("se_ivw"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.p_ivw, test.get("p_ivw"), XSD.double)

                        set_safe_literal(g, assessment_node, MTR.beta_ivw_r, test.get("beta_ivw_r"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.se_ivw_r, test.get("se_ivw_r"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.p_ivw_r, test.get("p_ivw_r"), XSD.double)

                        set_safe_literal(g, assessment_node, MTR.beta_pca, test.get("beta_pca"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.se_pca, test.get("se_pca"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.p_pca, test.get("p_pca"), XSD.double)

                        set_safe_literal(g, assessment_node, MTR.coloc_h1, test.get("PP.H1.abf"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.coloc_h2, test.get("PP.H2.abf"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.coloc_h3, test.get("PP.H3.abf"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.coloc_h4, test.get("PP.H4.abf"), XSD.double)

                        set_safe_literal(g, assessment_node, MTR.susie_h1, test.get("susie_max_PP.H1.abf"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.susie_h2, test.get("susie_max_PP.H2.abf"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.susie_h3, test.get("susie_max_PP.H3.abf"), XSD.double)
                        set_safe_literal(g, assessment_node, MTR.susie_h4, test.get("susie_max_PP.H4.abf"), XSD.double)

        elapsed = time.perf_counter() - entry_start
        logger.info("Finished %s in %.2fs", ratio_accession, elapsed)

        if (i + 1) % log_graph_every == 0:
            logger.info("Graph currently contains %d triples", len(g))

        if checkpoint_every and (i + 1) % checkpoint_every == 0:
            checkpoint_path = OUTPUT_DIR / f"rqtl_checkpoint_{i+1}.ttl"
            g.serialize(checkpoint_path, format="turtle")
            logger.info("Checkpoint saved to %s", checkpoint_path)

    logger.info("--- rQTL Integration Complete ---")
    logger.info("Final graph contains %d triples", len(g))
    return g


if __name__ == "__main__":
    graph = add_rqtl_to_graph("../data/json_files/GCST90200330_GCST90200020.json")
    output_path = OUTPUT_DIR / "rqtl.ttl"
    graph.serialize(output_path, format="turtle")
    logger.info("Serialized graph to %s", output_path)
