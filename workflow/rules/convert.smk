def vcf_file(chr):
    return str(Path(config["VCF_DIR"]) / config["VCF_FILE"].format(chr=chr))


rule make_sex_file:
    """Sex (1=male, 2=female) per sample from the covariates file, for --update-sex."""
    input:
        covariates=config["COVARIATES_FILE"],
        vcf=vcf_file(CHROMS[0]),
    output:
        sex=f"{QC}/sex.txt",
        summary=f"{QC}/sex.summary.tsv",
    params:
        id_col=config["COVARIATES_ID_COL"],
        sex_col=config["COVARIATES_SEX_COL"],
        min_match=config["MIN_SEX_MATCH"],
    log:
        f"{LOGS}/make_sex_file.log",
    shell:
        "python {SCRIPTS}/make_sex_file.py"
        " --covariates {input.covariates} --vcf {input.vcf}"
        " --id-col {params.id_col} --sex-col {params.sex_col}"
        " --min-match {params.min_match}"
        " --out {output.sex} --summary {output.summary} > {log} 2>&1"


rule vcf2pgen_b37:
    """VCF -> pgen with import QC (as in src/vcf2pfile.sh). PAR is split off chrX
    (b37 coordinates) so that male PAR genotypes are diploid during QC."""
    input:
        vcf=lambda w: vcf_file(w.chr),
        sex=rules.make_sex_file.output.sex,
    output:
        fileset(B37, PGEN) if config.get("KEEP_GRCH37") else temp(fileset(B37, PGEN)),
    params:
        plink_mem=plink_mem("vcf2pgen_b37"),
        prefix=B37,
        qc=config["QC"],
    threads: res("vcf2pgen_b37", "threads")
    resources:
        mem_mb=res("vcf2pgen_b37", "mem_mb"),
    log:
        f"{LOGS}/vcf2pgen_b37/chr{{chr}}.log",
    shell:
        "{PLINK2} --vcf {input.vcf} 'dosage={params.qc[dosage_field]}'"
        " --hard-call-threshold {params.qc[hard_call_threshold]}"
        " --exclude-if-info 'INFO<{params.qc[min_info]}'"
        " --maf {params.qc[maf]}"
        " --geno {params.qc[geno]}"
        " --hwe {params.qc[hwe]}"
        " --split-par b37"
        " --update-sex {input.sex}"
        " --make-pgen"
        " --threads {threads} --memory {params.plink_mem}"
        " --out {params.prefix} > {log} 2>&1"


rule pgen2bed:
    """PLINK1 binary file set (hard calls) from a pgen file set. Numeric chromosome
    codes (23=X, 25=PAR) for compatibility with PLINK1 tools."""
    input:
        multiext("{prefix}", *PGEN),
    output:
        multiext("{prefix}", *BED),
    params:
        plink_mem=plink_mem("pgen2bed"),
        prefix="{prefix}",
    threads: res("pgen2bed", "threads")
    resources:
        mem_mb=res("pgen2bed", "mem_mb"),
    log:
        "{prefix}.bed.log",
    shell:
        "{PLINK2} --pfile {params.prefix} --make-bed --output-chr 26"
        " --threads {threads} --memory {params.plink_mem}"
        " --out {params.prefix} > {log} 2>&1"
