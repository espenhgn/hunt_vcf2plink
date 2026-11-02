# In-sample PCA (as in src/imputation.stats.sh): LD-pruned common autosomal variants,
# KING relatedness, PCA in unrelated samples, and projection of all samples onto
# that basis. LD pruning is per chromosome (pruning windows never span chromosomes).

PCA_CFG = config["PCA"]
PCA_SOURCE = {"GRCh38": B38, "GRCh38.aligned": ALIGNED}[PCA_CFG["input"]]
PCA_CHROMS = [c for c in CHROMS if c.isdigit() and 1 <= int(c) <= 22]
PCA_WORK = f"{OUT}/tmp/pca"

if not PCA_CHROMS:
    raise ValueError("PCA needs at least one autosome in CHROMS")


def pca_keep(wildcards, input):
    """LD pruning for the PCA basis uses unrelated samples only, if PCA.unrelated."""
    return f"--keep {input.keep}" if "keep" in input.keys() else ""


def pca_prune_input(wildcards):
    inputs = {"pgen": fileset(PCA_SOURCE, PGEN, chr=wildcards.chr)}
    if wildcards.stage == "pca" and PCA_CFG["unrelated"]:
        inputs["keep"] = rules.king.output.unrelated
    return inputs


rule pca_ld_prune:
    """LD-pruned common variants of one chromosome (stage "king": all samples;
    stage "pca": unrelated samples), in a random subset of PCA.ld_samples samples."""
    input:
        unpack(pca_prune_input),
    output:
        pgen=temp(multiext(f"{PCA_WORK}/{{stage}}/chr{{chr}}", *PGEN)),
        prune_in=temp(f"{PCA_WORK}/{{stage}}/chr{{chr}}.prune.in"),
    wildcard_constraints:
        stage="king|pca",
        chr=r"[0-9]+",
    params:
        plink_mem=plink_mem("pca_ld_prune"),
        in_prefix=lambda wc: PCA_SOURCE.format(chr=wc.chr),
        out_prefix=f"{PCA_WORK}/{{stage}}/chr{{chr}}",
        keep=pca_keep,
        maf=PCA_CFG["maf"],
        geno=PCA_CFG["geno"],
        window=PCA_CFG["ld_window"],
        r2=PCA_CFG["ld_r2"],
        n_samples=PCA_CFG["ld_samples"],
        seed=PCA_CFG["seed"],
    threads: res("pca_ld_prune", "threads")
    resources:
        mem_mb=res("pca_ld_prune", "mem_mb"),
    log:
        f"{LOGS}/pca_ld_prune/{{stage}}.chr{{chr}}.log",
    shell:
        "({PLINK2} --pfile {params.in_prefix} {params.keep}"
        " --thin-indiv-count {params.n_samples} --seed {params.seed}"
        " --maf {params.maf} --geno {params.geno} --rm-dup exclude-all"
        " --indep-pairwise {params.window} {params.r2}"
        " --threads {threads} --memory {params.plink_mem}"
        " --out {params.out_prefix} && "
        "{PLINK2} --pfile {params.in_prefix} --extract {output.prune_in}"
        " --make-pgen"
        " --threads {threads} --memory {params.plink_mem}"
        " --out {params.out_prefix}) > {log} 2>&1"


rule pca_merge:
    """Merge the per-chromosome LD-pruned file sets of one stage."""
    input:
        lambda wc: [f"{PCA_WORK}/{wc.stage}/chr{c}{e}" for c in PCA_CHROMS for e in PGEN],
    output:
        pgen=temp(multiext(f"{PCA_WORK}/{{stage}}/merged", *PGEN)),
        mergelist=temp(f"{PCA_WORK}/{{stage}}/mergelist.txt"),
    wildcard_constraints:
        stage="king|pca",
    params:
        plink_mem=plink_mem("pca_merge"),
        prefixes=lambda wc: [f"{PCA_WORK}/{wc.stage}/chr{c}" for c in PCA_CHROMS],
        # --pmerge-list needs at least two file sets
        source=lambda wc, output: (
            f"--pmerge-list {output.mergelist} pfile --delete-pmerge-result"
            if len(PCA_CHROMS) > 1 else f"--pfile {PCA_WORK}/{wc.stage}/chr{PCA_CHROMS[0]}"),
        out_prefix=f"{PCA_WORK}/{{stage}}/merged",
    threads: res("pca_merge", "threads")
    resources:
        mem_mb=res("pca_merge", "mem_mb"),
    log:
        f"{LOGS}/pca_merge/{{stage}}.log",
    shell:
        "printf '%s\\n' {params.prefixes} > {output.mergelist} && "
        "{PLINK2} {params.source}"
        " --make-pgen"
        " --threads {threads} --memory {params.plink_mem}"
        " --out {params.out_prefix} > {log} 2>&1"


rule king:
    """KING kinship of sample pairs above PCA.king_cutoff (default: 3rd degree), and
    a maximal set of samples without such pairs (plink2 --king-cutoff-table)."""
    input:
        pgen=multiext(f"{PCA_WORK}/king/merged", *PGEN),
    output:
        kin0=f"{PCA}.kin0",
        unrelated=f"{PCA}.king.cutoff.in.id",
        related=f"{PCA}.king.cutoff.out.id",
        summary=f"{QC}/king.tsv",
    params:
        plink_mem=plink_mem("king"),
        in_prefix=f"{PCA_WORK}/king/merged",
        out_prefix=PCA,
        cutoff=PCA_CFG["king_cutoff"],
    threads: res("king", "threads")
    resources:
        mem_mb=res("king", "mem_mb"),
    log:
        f"{LOGS}/king.log",
    shell:
        "({PLINK2} --pfile {params.in_prefix}"
        " --make-king-table --king-table-filter {params.cutoff}"
        " --threads {threads} --memory {params.plink_mem}"
        " --out {params.out_prefix} && "
        "{PLINK2} --pfile {params.in_prefix}"
        " --king-cutoff-table {output.kin0} {params.cutoff}"
        " --threads {threads} --memory {params.plink_mem}"
        " --out {params.out_prefix}) > {log} 2>&1 && "
        "awk -v OFS='\\t' '!/^#/ {{n[FILENAME]++}}"
        " END {{print \"king_related_pairs\", n[ARGV[1]] + 0;"
        " print \"king_unrelated_samples\", n[ARGV[2]] + 0;"
        " print \"king_excluded_samples\", n[ARGV[3]] + 0}}'"
        " {output.kin0} {output.unrelated} {output.related} > {output.summary}"


def pca_basis_input(wildcards):
    inputs = {"pgen": multiext(f"{PCA_WORK}/pca/merged", *PGEN)}
    if PCA_CFG["unrelated"]:
        inputs["keep"] = rules.king.output.unrelated
    return inputs


rule pca_basis:
    """PCA (variant weights and allele counts) in unrelated samples, or all samples
    if not PCA.unrelated."""
    input:
        unpack(pca_basis_input),
    output:
        eigenvec=f"{PCA}.pca_basis.eigenvec",
        eigenval=f"{PCA}.pca_basis.eigenval",
        eigenvec_allele=f"{PCA}.pca_basis.eigenvec.allele",
        acount=f"{PCA}.pca_basis.acount",
        summary=f"{QC}/pca_basis.tsv",
    params:
        plink_mem=plink_mem("pca_basis"),
        in_prefix=f"{PCA_WORK}/pca/merged",
        out_prefix=f"{PCA}.pca_basis",
        keep=pca_keep,
        approx="approx" if PCA_CFG["approx"] else "",
        n_pc=config["N_PC"],
    threads: res("pca_basis", "threads")
    resources:
        mem_mb=res("pca_basis", "mem_mb"),
    log:
        f"{LOGS}/pca_basis.log",
    shell:
        "{PLINK2} --pfile {params.in_prefix} {params.keep}"
        " --freq counts"
        " --pca allele-wts {params.approx} {params.n_pc} vcols=chrom,ref,alt"
        " --threads {threads} --memory {params.plink_mem}"
        " --out {params.out_prefix} > {log} 2>&1 && "
        "awk -v OFS='\\t' '!/^#/ {{n[FILENAME]++}}"
        " END {{print \"pca_variants\", n[ARGV[1]] + 0; print \"pca_basis_samples\", n[ARGV[2]] + 0}}'"
        " {params.in_prefix}.pvar {output.eigenvec} > {output.summary}"


rule pca_project:
    """Project all samples onto the PCA basis."""
    input:
        pgen=multiext(f"{PCA_WORK}/pca/merged", *PGEN),
        acount=rules.pca_basis.output.acount,
        eigenvec_allele=rules.pca_basis.output.eigenvec_allele,
    output:
        f"{PCA}.sscore",
    params:
        plink_mem=plink_mem("pca_project"),
        in_prefix=f"{PCA_WORK}/pca/merged",
        out_prefix=PCA,
        last_col=config["N_PC"] + 5,
    threads: res("pca_project", "threads")
    resources:
        mem_mb=res("pca_project", "mem_mb"),
    log:
        f"{LOGS}/pca_project.log",
    shell:
        "{PLINK2} --pfile {params.in_prefix}"
        " --read-freq {input.acount}"
        " --score {input.eigenvec_allele} 2 5 'header-read' 'no-mean-imputation' 'variance-standardize'"
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
