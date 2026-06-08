#!/usr/bin/env snakemake -s

# Snakefile for ntRoot pipeline
import os
import shutil

onsuccess:
    shutil.rmtree(".snakemake", ignore_errors=True)

ruleorder: exome_mask_fasta > exome_masked_provided

# Read parameters from config or set default values
reference=config["reference"]
draft_base = os.path.basename(os.path.realpath(reference))
reads_prefix_full=config["reads"] if "reads" in config else ""
reads_prefix=os.path.basename(reads_prefix_full)
k=config["kmer"] if "kmer" in config else None

genomes = config["genomes"] if "genomes" in config else ""
genome_prefix = ".".join([os.path.basename(os.path.realpath(genome)).removesuffix(".fa").removesuffix(".fasta").removesuffix(".fna") for genome in genomes])

# Common parameters
t = config["threads"] if "threads" in config else 4
b = config["b_prefix"] + "_" if "b_prefix" in config and config["b_prefix"] != "" else ""


# ntEdit parameters
z = config["z_param"] if "z_param" in config else 100
v = config["verbose"] if "verbose" in config else 0
j = config["j_param"] if "j_param" in config else 3
Y = config["Y_param"] if "Y_param" in config else 0.55
l = config["l_vcf"] if "l_vcf" in config else ""
cutoff = config["cutoff"] if "cutoff" in config else 0

# Ancestry inference parameters
tile_size = config["tile_size"] if "tile_size" in config else 5000000

# Third-party VCF parameters
input_vcf = config["input_vcf"] if "input_vcf" in config else None
input_vcf_basename = os.path.basename(os.path.realpath(input_vcf)) if input_vcf else "None"
strip_info = config["strip_info"] if "strip_info" in config else None

# Exome parameters
exome = config["exome"] if "exome" in config else False
masked = config["masked"] if "masked" in config else False
exome_bed = config["exome_bed"] if "exome_bed" in config else ""
exome_bed_base = os.path.basename(os.path.realpath(exome_bed)) if exome_bed else "None"

# time command
mac_time_command = "command time -l -o"
linux_time_command = "command time -v -o"
time_command = mac_time_command if os.uname().sysname == "Darwin" else linux_time_command

rule all:
    input: f"{reads_prefix}_ntedit_k{k}_variants.vcf_ancestry-predictions_tile{tile_size}.tsv"

rule ntroot_genome:
    input: f"{genome_prefix}_ntedit_k{k}_variants.vcf_ancestry-predictions_tile{tile_size}.tsv"

rule ntroot_reads:
    input: f"{reads_prefix}_ntedit_k{k}_variants.vcf_ancestry-predictions_tile{tile_size}.tsv"

rule ntroot_reads_exome:
    input: f"{reads_prefix}_ntedit_k{k}_exome_variants.vcf_ancestry-predictions_tile{tile_size}.tsv"

rule ntroot_genome_lai:
    input: f"{genome_prefix}_ntedit_k{k}_variants.vcf_ancestry-predictions-tile-resolution_tile{tile_size}.tsv"

rule ntroot_reads_lai:
    input: f"{reads_prefix}_ntedit_k{k}_variants.vcf_ancestry-predictions-tile-resolution_tile{tile_size}.tsv"

rule ntroot_reads_exome_lai:
    input: f"{reads_prefix}_ntedit_k{k}_exome_variants.vcf_ancestry-predictions-tile-resolution_tile{tile_size}.tsv"

rule ntroot_input_vcf:
    input: f"{input_vcf_basename}.cross-ref.vcf_ancestry-predictions_tile{tile_size}.tsv"

rule ntroot_input_vcf_lai:
    input: f"{input_vcf_basename}.cross-ref.vcf_ancestry-predictions-tile-resolution_tile{tile_size}.tsv"

rule ntedit_reads:
    input: 
        reference = reference
    output:
        out_vcf = f"{reads_prefix}_ntedit_k{k}_variants.vcf",
        out_fa = temp(f"{reads_prefix}_ntedit_k{k}_edited.fa"),
        out_changes = temp(f"{reads_prefix}_ntedit_k{k}_changes.tsv"),
        out_bf = temp(f"{reads_prefix}_k{k}.bf")
    params:
        benchmark = f"{time_command} ntedit_snv_k{k}.time",
        params = f"-k {k} -t {t} -z {z} -j {j} -Y {Y}",
        cutoff = f"--cutoff {cutoff}" if cutoff > 0 else "--solid",
        vcf_input = f"-l {l}" if l else ""
    shell:
        "{params.benchmark} run-ntedit snv --reference {input.reference} --reads {reads_prefix_full} {params.params} "
        "{params.vcf_input} {params.cutoff}"

rule ntedit_exome_reads:
    input: 
        reference = f"masked_{draft_base}"
    output:
        out_vcf = f"{reads_prefix}_ntedit_k{k}_exome_variants.vcf",
        out_fa = temp(f"{reads_prefix}_ntedit_k{k}_edited.fa"),
        out_changes = temp(f"{reads_prefix}_ntedit_k{k}_changes.tsv"),
        out_bf = temp(f"{reads_prefix}_k{k}.bf")
    params:
        benchmark = f"{time_command} ntedit_snv_k{k}.time",
        params = f"-k {k} -t {t} -z {z} -j {j} -Y {Y}",
        cutoff = f"--cutoff {cutoff}" if cutoff > 0 else "--solid",
        vcf_input = f"-l {l}" if l else "",
        out_vcf_tmp = f"{reads_prefix}_ntedit_k{k}_variants.vcf"
    shell:
        """
        {params.benchmark} run-ntedit snv --reference {input.reference} --reads {reads_prefix_full} {params.params} {params.vcf_input} {params.cutoff}
        mv {params.out_vcf_tmp} {output.out_vcf}
        """

rule samtools_faidx:
    input: reference = reference
    output: out_fai = f"{draft_base}.fai"
    params:
        benchmark = f"{time_command} samtools_faidx_{draft_base}.time"
    run:
        if input.reference.endswith(".gz"):
            shell("{params.benchmark} gunzip -c {input.reference} |samtools faidx -o {output.out_fai} -")
        else:
            shell("{params.benchmark} samtools faidx -o {output.out_fai} {input.reference}")

rule exome_masked_provided:
    input: 
        reference = reference
    output:
        masked_reference = f"masked_{draft_base}"
    shell:
        "ln -sf {input.reference} {output.masked_reference}"

rule exome_sort_bed:
    input:
        bed = exome_bed,
        ref_fai = rules.samtools_faidx.output
    output:
        sorted_bed = temp(f"{exome_bed_base}.sorted.bed")
    shell:
        "bedtools sort -i {input.bed} -faidx {input.ref_fai} > {output.sorted_bed}"

rule exome_complement_bed:
    input:
        sorted_bed = f"{exome_bed_base}.sorted.bed",
        ref_fai = rules.samtools_faidx.output
    output:
        complement_bed = temp(f"{exome_bed_base}.sorted.bed.complement.bed")
    shell:
        "bedtools complement -i {input.sorted_bed} -g {input.ref_fai} > {output.complement_bed}"

rule exome_mask_fasta:
    input:
        reference = reference,
        complement_bed = rules.exome_complement_bed.output
    output:
        masked_reference = f"masked_{draft_base}"
    run:
        if input.reference.endswith(".gz"):
            shell("bedtools maskfasta -fi <(gunzip -c {input.reference}) -bed {input.complement_bed} -fo tmp_{output.masked_reference}")
            shell("gzip -c tmp_{output.masked_reference} > {output.masked_reference}")
            shell("rm tmp_{output.masked_reference}")
        else:
            shell("bedtools maskfasta -fi {input.reference} -bed {input.complement_bed} -fo {output.masked_reference}")

rule ntedit_genome:
    input: 
        reference = reference,
        genomes = genomes
    output:
        out_vcf = f"{genome_prefix}_ntedit_k{k}_variants.vcf",
        out_fa = temp(f"{genome_prefix}_ntedit_k{k}_edited.fa"),
        out_changes = temp(f"{genome_prefix}_ntedit_k{k}_changes.tsv"),
        out_bf = temp(f"{genome_prefix}_k{k}.bf")
    params:
        benchmark = f"{time_command} ntedit_snv_k{k}.time",
        params = f"-k {k} -t {t} -z {z} -j {j} -Y {Y}",
        vcf_input = f"-l {l}" if l else ""
    shell:
        "{params.benchmark} run-ntedit snv --reference {reference} --genome {input.genomes} {params.params} "
        " {params.vcf_input}"


rule ancestry_prediction:
    input: 
        vcf = "{vcf}"
    output: 
        predictions = "{vcf}_ancestry-predictions_tile{tile_size}.tsv"
    params:
        benchmark = f"{time_command} ancestry_prediction_tile{tile_size}.time" if input_vcf_basename else f"{time_command} ancestry_prediction_k{k}_tile{tile_size}.time",
        tile_size = tile_size,
        verbosity = v
    shell:
        "{params.benchmark} ntRootAncestryPredictor.pl -f {input.vcf} -t {params.tile_size} -v {params.verbosity}"


rule ancestry_prediction_lai:
    input: 
        vcf = "{vcf}",
        ref_fai = f"{draft_base}.fai"
    output: 
        lai_output = "{vcf}_ancestry-predictions-tile-resolution_tile{tile_size}.tsv"
    params:
        benchmark = f"{time_command} ancestry_prediction_tile{tile_size}.time" if input_vcf_basename else f"{time_command} ancestry_prediction_k{k}_tile{tile_size}.time",
        tile_size = tile_size,
        verbosity = v
    shell:
        "{params.benchmark} ntRootAncestryPredictor.pl -f {input.vcf} -t {params.tile_size} -v {params.verbosity} -r 1 -i {input.ref_fai}"

rule sort_vcf_input:
    input: vcf = f"{input_vcf}"
    output: vcf_sorted = temp(f"{input_vcf_basename}_sorted.vcf")
    params:
        benchmark = f"{time_command} sort_vcf_{input_vcf_basename}.time",
        cat_cmd = "gunzip -c" if f"{input_vcf}".endswith(".gz") else "cat"
    shell:
        """{params.benchmark} sh -c '(echo "##fileformat=VCFv4.2" ; {params.cat_cmd} {input.vcf} |grep -v "^#" |sort -k1,1 -k2,2n) > {output.vcf_sorted}'"""

rule sort_vcf_l:
    input: vcf = l
    output: temp(f"{os.path.basename(os.path.realpath(l))}_sorted.tmp.vcf")
    params:
        benchmark = f"{time_command} sort_vcf_l.time",
        cat_cmd = "gunzip -c" if f"{l}".endswith(".gz") else "cat"
    shell:
        """{params.benchmark} sh -c "(echo '##fileformat=VCFv4.2' ; {params.cat_cmd} {input.vcf} | awk '\$5 !~ /^</' | grep -v '^#' |sort -k1,1 -k2,2n) > {output}" """

rule bedtools_intersect:
    input:         
        sorted_vcf = f"{input_vcf_basename}_sorted.vcf",
        sorted_ref_vars = f"{os.path.basename(os.path.realpath(l))}_sorted.tmp.vcf"
    output:
        bedtools = temp(f"{input_vcf_basename}.bedtools-intersect.bed")
    params:
        benchmark = f"{time_command} bedtools_intersect_{input_vcf_basename}.time"
    shell:
        "{params.benchmark} bedtools intersect -loj -sorted -a {input.sorted_vcf} -b {input.sorted_ref_vars} > {output.bedtools}"

rule cross_reference_vcf:
    input: 
        vcf = f"{input_vcf}",
        ref_vars = l,
        bedtools = f"{input_vcf_basename}.bedtools-intersect.bed"
    output: f"{input_vcf_basename}.cross-ref.vcf"
    params:
        benchmark = f"{time_command} cross_reference_vcf_{input_vcf_basename}.time",
        prefix=f"{input_vcf_basename}.cross-ref",
        strip = "--strip" if strip_info else ""
    shell: 
        "{params.benchmark} ntroot_cross_reference_vcf.py -b {input.bedtools} --vcf {input.vcf} --vcf_l {input.ref_vars} -p {params.prefix} {params.strip}"