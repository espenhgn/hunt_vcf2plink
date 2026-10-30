rule align_b38:
    """SNPs aligned to the GRCh38 reference sequence, with chr:pos:ref:alt IDs
    (as in src/align_ids_toGRCh38.sh)."""
    input:
        pgen=fileset(B38, PGEN),
        fasta=config["REF_FASTA_GRCH38"],
    output:
        pgen=fileset(ALIGNED, PGEN),
        summary=f"{QC}/align/chr{{chr}}.tsv",
    params:
        plink_mem=plink_mem("align_b38"),
        in_prefix=B38,
        out_prefix=ALIGNED,
        work=f"{WORK}/align_b38/chr{{chr}}",
    threads: res("align_b38", "threads")
    resources:
        mem_mb=res("align_b38", "mem_mb"),
    log:
        f"{LOGS}/align_b38/chr{{chr}}.log",
    shell:
        "bash {SCRIPTS}/align.sh {PLINK2} {input.fasta}"
        " {params.in_prefix} {params.out_prefix} {wildcards.chr}"
        " {threads} {params.plink_mem} {params.work} {output.summary} > {log} 2>&1"
