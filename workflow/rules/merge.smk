rule merge_aligned:
    """Merge per-chromosome aligned file sets (as in src/merge_chr_GRCh38.sh, but
    writing pgen, which keeps dosages; .bed is made from it by pgen2bed)."""
    input:
        [fileset(ALIGNED, PGEN, chr=c) for c in CHROMS],
    output:
        pgen=fileset(MERGED, PGEN),
        mergelist=f"{MERGED}.mergelist.txt",
    params:
        plink_mem=plink_mem("merge_aligned"),
        prefixes=[ALIGNED.format(chr=c) for c in CHROMS],
        out_prefix=MERGED,
    threads: res("merge_aligned", "threads")
    resources:
        mem_mb=res("merge_aligned", "mem_mb"),
    log:
        f"{LOGS}/merge_aligned.log",
    shell:
        "printf '%s\\n' {params.prefixes} > {output.mergelist} && "
        "{PLINK2} --pmerge-list {output.mergelist} pfile"
        " --delete-pmerge-result"
        " --make-pgen"
        " --threads {threads} --memory {params.plink_mem}"
        " --out {params.out_prefix} > {log} 2>&1"
