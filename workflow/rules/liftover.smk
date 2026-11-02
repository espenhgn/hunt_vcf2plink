LIFTOVER_DIR = f"{OUT}/resources/liftover"


localrules:
    setup_liftover,


rule setup_liftover:
    """UCSC liftOver (official build for this platform) and the hg19ToHg38 chain file."""
    output:
        liftover=f"{LIFTOVER_DIR}/liftOver",
        chain=f"{LIFTOVER_DIR}/hg19ToHg38.over.chain.gz",
    params:
        liftover=config.get("LIFTOVER_BINARY") or "",
        liftover_md5=config.get("LIFTOVER_MD5") or "",
        url_base=config["LIFTOVER_URL_BASE"],
        platforms=config["LIFTOVER_PLATFORMS"],
        chain=config["CHAIN_FILE"],
        chain_md5=config.get("CHAIN_MD5") or "",
    log:
        f"{LOGS}/setup_liftover.log",
    script:
        "../scripts/setup_liftover.py"


rule liftover_b38:
    """GRCh37 -> GRCh38 with UCSC liftOver, strand flipping SNVs in inverted regions
    (triple-liftOver method; as in src/plink2_liftover.sh)."""
    input:
        pgen=fileset(B37, PGEN),
        liftover=rules.setup_liftover.output.liftover,
        chain=rules.setup_liftover.output.chain,
    output:
        pgen=temp(fileset(LIFTED, PGEN)),
        summary=f"{QC}/liftover/chr{{chr}}.tsv",
    params:
        plink_mem=plink_mem("liftover_b38"),
        in_prefix=B37,
        out_prefix=LIFTED,
        work=f"{WORK}/liftover_b38/chr{{chr}}",
    threads: res("liftover_b38", "threads")
    resources:
        mem_mb=res("liftover_b38", "mem_mb"),
    log:
        f"{LOGS}/liftover_b38/chr{{chr}}.log",
    shell:
        "bash {SCRIPTS}/liftover.sh {PLINK2} {input.liftover} {input.chain}"
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
