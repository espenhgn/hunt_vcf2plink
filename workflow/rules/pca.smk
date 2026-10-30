PCA_REF = repo_path(config["PCA_REF_PREFIX"])


rule pca_overlap:
    """Reference PCA variants present in the data. Fails below PCA_MIN_OVERLAP,
    which catches build or variant ID mismatches."""
    input:
        pvar=f"{MERGED}.pvar",
        eigenvec=f"{PCA_REF}.eigenvec.allele",
    output:
        ids=f"{PCA}.pca_overlap.ids",
        summary=f"{QC}/pca_overlap.tsv",
    params:
        min_overlap=config["PCA_MIN_OVERLAP"],
    log:
        f"{LOGS}/pca_overlap.log",
    shell:
        "python {SCRIPTS}/pca_overlap.py --pvar {input.pvar} --eigenvec {input.eigenvec}"
        " --min-overlap {params.min_overlap}"
        " --out {output.ids} --summary {output.summary} > {log} 2>&1"


rule pca_project:
    """Project samples onto the reference PCA basis (as in src/pca/inhouse_pca_insample.sh)."""
    input:
        pgen=fileset(MERGED, PGEN),
        ids=rules.pca_overlap.output.ids,
        acount=f"{PCA_REF}.acount",
        eigenvec=f"{PCA_REF}.eigenvec.allele",
    output:
        f"{PCA}.sscore",
    params:
        plink_mem=plink_mem("pca_project"),
        in_prefix=MERGED,
        out_prefix=PCA,
        last_col=config["N_PC"] + 5,
    threads: res("pca_project", "threads")
    resources:
        mem_mb=res("pca_project", "mem_mb"),
    log:
        f"{LOGS}/pca_project.log",
    shell:
        "{PLINK2} --pfile {params.in_prefix}"
        " --extract {input.ids}"
        " --read-freq {input.acount}"
        " --score {input.eigenvec} 2 5 'header-read' 'no-mean-imputation' 'variance-standardize'"
        " --score-col-nums 6-{params.last_col}"
        " --threads {threads} --memory {params.plink_mem}"
        " --out {params.out_prefix} > {log} 2>&1"


rule pcs_table:
    """IID and PC1..PCn as a plain TSV."""
    input:
        rules.pca_project.output,
    output:
        f"{PCA}.pcs.tsv",
    params:
        n_pc=config["N_PC"],
    log:
        f"{LOGS}/pcs_table.log",
    shell:
        "python {SCRIPTS}/pcs_table.py {input} {params.n_pc} > {output} 2> {log}"
