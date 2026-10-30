rule report:
    """Variant and sample counts per step, one row per chromosome."""
    input:
        sex=rules.make_sex_file.output.summary,
        liftover=[f"{QC}/liftover/chr{c}.tsv" for c in CHROMS],
        align=[f"{QC}/align/chr{c}.tsv" for c in CHROMS],
        b37_logs=[f"{LOGS}/vcf2pgen_b37/chr{c}.log" for c in CHROMS],
        merged=f"{MERGED}.pvar",
        pca=rules.pca_overlap.output.summary,
    output:
        f"{QC}/report.tsv",
    log:
        f"{LOGS}/report.log",
    shell:
        "python {SCRIPTS}/report.py --sex {input.sex} --liftover {input.liftover}"
        " --align {input.align} --b37-logs {input.b37_logs}"
        " --merged-pvar {input.merged} --pca {input.pca} > {output} 2> {log}"
