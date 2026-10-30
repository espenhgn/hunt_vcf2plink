rule liftover_b38:
    """GRCh37 -> GRCh38 with triple-liftOver (as in src/plink2_liftover.sh)."""
    input:
        fileset(B37, PGEN),
    output:
        pgen=temp(fileset(LIFTED, PGEN)),
        summary=f"{QC}/liftover/chr{{chr}}.tsv",
    params:
        plink_mem=plink_mem("liftover_b38"),
        in_prefix=B37,
        out_prefix=LIFTED,
        triple_liftover=repo_path(config["TRIPLE_LIFTOVER"]),
        work=f"{WORK}/liftover_b38/chr{{chr}}",
    threads: res("liftover_b38", "threads")
    resources:
        mem_mb=res("liftover_b38", "mem_mb"),
    log:
        f"{LOGS}/liftover_b38/chr{{chr}}.log",
    shell:
        "bash {SCRIPTS}/liftover.sh {PLINK2} {params.triple_liftover}"
        " {params.in_prefix} {params.out_prefix} {wildcards.chr}"
        " {threads} {params.plink_mem} {params.work} {output.summary} > {log} 2>&1"


rule finalize_b38:
    """GRCh38 file set with all lifted variants. Variant IDs are set to
    chr:pos:ref:alt in GRCh38 coordinates (alleles as in the input, not aligned to
    the reference sequence). The PAR is split off chrX again, in b38 coordinates;
    IDs are assigned before the split, so PAR variant IDs start with 'X:'."""
    input:
        fileset(LIFTED, PGEN),
    output:
        fileset(B38, PGEN),
    params:
        plink_mem=plink_mem("finalize_b38"),
        in_prefix=LIFTED,
        out_prefix=B38,
        split_par=lambda w: "--split-par b38" if w.chr == "X" else "",
    threads: res("finalize_b38", "threads")
    resources:
        mem_mb=res("finalize_b38", "mem_mb"),
    log:
        f"{LOGS}/finalize_b38/chr{{chr}}.log",
    shell:
        "{PLINK2} --pfile {params.in_prefix}"
        " --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing"
        " {params.split_par}"
        " --make-pgen"
        " --threads {threads} --memory {params.plink_mem}"
        " --out {params.out_prefix} > {log} 2>&1"
