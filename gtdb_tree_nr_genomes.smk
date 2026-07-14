# =========================
# GTDB-Tk: classify -> de novo
# =========================

configfile: "config/gtdb_config.yml"

# conda env config vars set GTDBTK_DATA_PATH=/data/san/data0/databases/gtdbtk/release226/release226

from pathlib import Path

# ---- load config ----
BASE_DIR   = Path(config["base_dir"])
GENOME_DIR = BASE_DIR / config["genomes_dir"]

THREADS        = config.get("threads", 8)
OUTGROUP_TAXON = config["outgroup_taxon"]
CONDA_ENV      = config["conda_env"]

# ---- derived paths ----
CLASSIFY_DIR = BASE_DIR / "classify"
DENOVO_DIR   = BASE_DIR / "denovo_out"

SUMMARY      = CLASSIFY_DIR / "gtdbtk.bac120.summary.tsv"
CLASSIFY_LOG = CLASSIFY_DIR / "gtdbtk.log"

CUSTOM_TAX   = BASE_DIR / "custom_taxonomy.txt"
DENOVO_LOG   = DENOVO_DIR / "gtdbtk.log"

# =========================
# final targets
# =========================
rule all:
    input:
        str(SUMMARY),
        str(CUSTOM_TAX),
        str(DENOVO_LOG)

# =========================
# GTDB-Tk classify
# =========================
rule gtdbtk_classify:
    input:
        genomes_dir=lambda wc: directory(str(GENOME_DIR))
    output:
        summary=str(SUMMARY),
        log=str(CLASSIFY_LOG)
    params:
        out=str(CLASSIFY_DIR)
    threads:
        THREADS
    conda:
        CONDA_ENV
    shell:
        r"""
		export GTDBTK_DATA_PATH={config[gtdbtk_data_path]}
        mkdir -p {params.out}
        gtdbtk classify_wf \
          --skip_ani_screen \
          --genome_dir {input.genomes_dir} \
          --out_dir {params.out} \
          --cpus {threads}
        """

# =========================
# build custom taxonomy
# =========================
rule make_custom_taxonomy:
    input:
        summary=str(SUMMARY)
    output:
        tax=str(CUSTOM_TAX)
    shell:
        r"""
        cut -f1,2 {input.summary} > {output.tax}
        """

# =========================
# GTDB-Tk de novo tree
# =========================
rule gtdbtk_denovo:
    input:
        genomes_dir=lambda wc: directory(str(GENOME_DIR)),
        custom_tax=str(CUSTOM_TAX)
    output:
        log=str(DENOVO_LOG)
    params:
        out=str(DENOVO_DIR),
        outgroup=OUTGROUP_TAXON
    threads:
        THREADS
    conda:
        CONDA_ENV
    shell:
        r"""
		export GTDBTK_DATA_PATH={config[gtdbtk_data_path]}
        mkdir -p {params.out}
        gtdbtk de_novo_wf \
          --genome_dir {input.genomes_dir} \
          --out_dir {params.out} \
          --bacteria \
          --outgroup_taxon "{params.outgroup}" \
          --skip_gtdb_refs \
          --custom_taxonomy_file {input.custom_tax} \
          --cpus {threads}
        """
