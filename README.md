# mtrKG

mtrKG is a **Metabolite Ratio Knowledge Graph** project that integrates local rQTL outputs with public resources (GWAS Catalog, Open Targets, STRING, Ensembl Regulatory Build, Reactome, HMDB, Rhea), then supports SPARQL analysis, SHACL validation, and drug repurposing experiments.

## What this repository contains

1. A schema and integration pipeline that builds RDF/Turtle knowledge graphs.
2. Analysis queries (`.rq`) and notebook workflows for hypothesis-driven exploration.
3. Validation shapes and reports (PySHACL).
4. Optional ML-based link prediction for drug repurposing (PyKEEN/TransE).

## End-to-end workflow (high level)

1. Ingest local rQTL JSON files into an RDF graph.
2. Enrich the graph with external biomedical sources.
3. Serialize the graph to Turtle (`output/mtrKG_01.ttl` / `output/mtrKG.ttl`).
4. Load the graph into GraphDB (repository usually named `mtrKG`).
5. Run SPARQL analysis from `analysis/*.rq`.
6. Save query outputs as CSV files for downstream reporting.
7. Validate graph quality with SHACL constraints.

## Repository map

| Path | Purpose | Key contents |
|---|---|---|
| `src/` | Core code and notebooks | Integrators, schema, SPARQL engines, build/analysis/validation notebooks |
| `analysis/` | SPARQL and rule assets | `002_query.rq` ... `014_query.rq`, `q10_query_construct.rq`, `rqtl_ruleset.pie`, generated `*.csv` |
| `data/` | Source datasets and DB artifacts | `json_files/` (rQTL JSONs), GWAS/HMDB exports, SQLite assets |
| `output/` | Generated outputs | KG Turtle files, integration logs/reports, drug predictions, validation outputs |
| `validation/` | SHACL shapes and reports | `rqtl_shapes.ttl`, `schema_validation.ttl`, validation reports |
| `doc/` | Documentation assets | Drawings, pathway visuals, exported HTML/PDF diagrams |
| `notebook/` | Auxiliary notebook area | Additional/legacy notebook content |
| `sandbox/` | Experimental work area | Scratch notebooks, prototypes, intermediate artifacts |

## Main notebooks

| Notebook | Purpose |
|---|---|
| `src/create_mtrKG.ipynb` | Main KG build and enrichment pipeline |
| `src/analyse_graph.ipynb` | Runs SPARQL analyses over GraphDB and exports query outputs |
| `src/validation.ipynb` | Runs PySHACL validation and writes reports |
| `src/utility.ipynb` | Utility analysis/visualization snippets |
| `src/notebook.ipynb`, `src/notebook_01.ipynb` | Broader exploratory/legacy workflows |

## Main Python modules

| Module | Responsibility |
|---|---|
| `src/schema_definition.py` | Declares namespaces, classes, and predicates; builds ontology/schema scaffold |
| `src/integrate_rQTLs.py` | Ingests local rQTL JSON and creates core ratio/variant/causal structures |
| `src/integrate_gwas_catalog.py` | Enriches SNPs with GWAS Catalog associations |
| `src/integrate_open_targets.py` | Adds target tractability, diseases, liabilities, and known drugs |
| `src/integrate_string.py` | Adds gene-gene interaction edges from STRING |
| `src/integrate_encode.py` | Adds SNP overlap with Ensembl regulatory/motif features |
| `src/integrate_reactome.py` | Adds gene/metabolite pathway context from Reactome |
| `src/integrate_HMDB.py` | Adds HMDB metabolite location knowledge |
| `src/integrate_rhea.py` | Adds reaction participation via Rhea SPARQL |
| `src/integrate_ewas.py` | Optional EWAS enrichment module |
| `src/graphdb_engine.py` | Executes SPARQL against GraphDB and returns DataFrames |
| `src/execute_sparql.py` | Executes SPARQL directly on local RDFLib graph |
| `src/predict_drugs.py` | Trains TransE embeddings and generates drug repurposing candidates |

## Data assets

- `data/json_files/`: local rQTL JSON inputs (5,095 files in this workspace snapshot).
- `data/hmdb_metabolites.xml`, `data/hmdb_proteins.xml`: HMDB bulk XML resources.
- `data/gwas-catalog-download-associations-alt-full.tsv`: GWAS association table.
- `data/Human-GEM.xml`: metabolic model resource.
- `data/instance/metabolite-ratio-app.sqlite`: local app database asset.
- `data/populate_db_with_jsons.py`: helper script for filling the SQLite app schema.

## Setup

There is currently no pinned `requirements.txt` in the repository root. A practical environment for the current code typically needs:

- Python
- `rdflib`, `pandas`, `requests`, `SPARQLWrapper`
- `pyshacl` (validation)
- `numpy`, `scipy`, `torch`, `pykeen` (drug prediction)
- `tqdm` (notebook progress bars)
- `networkx`, `matplotlib` (utility visualization notebook)

Example installation:

```bash
pip install rdflib pandas requests SPARQLWrapper pyshacl numpy scipy torch pykeen tqdm networkx matplotlib
```

## Building the KG

Recommended entry point: `src/create_mtrKG.ipynb`.

The notebook pipeline applies integrations in this sequence:

1. rQTL local JSON ingestion
2. GWAS Catalog
3. Open Targets
4. STRING
5. Ensembl Regulatory Build
6. Reactome (genes, then metabolites)
7. HMDB
8. Rhea
9. Serialize graph to Turtle
10. Optional drug repurposing

Typical output files are written under `output/` and `output/integration/`:

- `mtrKG.ttl` / `mtrKG_01.ttl`
- `*_integration.log`
- `*_mapping_report.csv`
- `drug_repurposing_predictions.csv`

## SPARQL analysis workflow

Primary analysis notebook: `src/analyse_graph.ipynb`.

Queries are stored in `analysis/*.rq` and executed via `query_graphdb()` in `src/graphdb_engine.py`.

Current behavior:

1. Reads each `.rq` file.
2. Prints the exact SPARQL query text being executed.
3. Executes the query on GraphDB (`repo_name="mtrKG"` by default).
4. Saves SELECT results as CSV with the same base name:
   - `analysis/002_query.rq` -> `analysis/002_query.csv`
   - `analysis/013_query_count_triples.rq` -> `analysis/013_query_count_triples.csv`

GraphDB prerequisites:

- GraphDB is running locally.
- Repository exists (usually `mtrKG`).
- The generated Turtle graph is loaded into that repository.

## Validation workflow

Validation assets:

- `validation/rqtl_shapes.ttl`
- `validation/schema_validation.ttl`

Validation execution is shown in `src/validation.ipynb` using `pyshacl.validate(...)`.

Typical outputs:

- `output/validation_report.ttl`
- `output/validation_text.txt`

## Rule-based reasoning asset

`analysis/rqtl_ruleset.pie` contains rule templates for deriving additional relations, including:

- putative SNP-to-gene implication
- ratio-to-pathway inference
- putative disease implication chains
- candidate therapeutic target links

This file is intended for rule-engine workflows (for example, GraphDB rule sets).

## Drug repurposing workflow

`src/predict_drugs.py`:

1. Loads structural triples from Turtle.
2. Trains a TransE embedding model with PyKEEN.
3. Scores disease-drug proximity in embedding space.
4. Writes predictions to `output/drug_repurposing_predictions.csv`.

## Notes and practical considerations

- Several integration scripts call external APIs; internet access and rate limits matter.
- Some outputs and logs are large; avoid committing regenerated artifacts unintentionally.
- Relative paths in notebooks are usually written with execution from the `src/` context.
- If you run outside notebooks, ensure your working directory keeps file paths consistent.
