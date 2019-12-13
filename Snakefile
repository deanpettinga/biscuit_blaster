import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.4.4")

#-------------------------------------------------------------------------------
rule all:
    input:
        expand("analysis/{sample}.sorted.markdup.bam", sample=config["sample"]),
        expand("analysis/{sample}.sorted.markdup.bam.bai", sample=config["sample"]),
        expand("analysis/{sample}.sorted.markdup.bam.flagstat", sample=config["sample"]),

#-------------------------------------------------------------------------------

rule index:
    params:
        config["reference"]
    output:
        "reference/ref.fa",
        "reference/ref.fa.bis.amb",
        "reference/ref.fa.bis.ann",
        "reference/ref.fa.bis.pac",
        "reference/ref.fa.dau.bwt",
        "reference/ref.fa.dau.sa",
        "reference/ref.fa.par.bwt",
        "reference/ref.fa.par.sa",
    log:
        "logs/wget_default_index.log",
        "logs/biscuit_index.log",
    singularity:
        "docker://dpettinga/biscuit_blaster"
    wrapper:
        "file:wrappers/biscuit_index"

#-------------------------------------------------------------------------------
rule biscuit_blaster:
    input:
        fq1 = expand("raw_data/{sample}", sample=config["fq1"]),
        ref = "reference/ref.fa",
        bis_amb = "reference/ref.fa.bis.amb",
        bis_ann = "reference/ref.fa.bis.ann",
        bis_pac = "reference/ref.fa.bis.pac",
        dau_bwt = "reference/ref.fa.dau.bwt",
        dau_sa = "reference/ref.fa.dau.sa",
        par_bwt = "reference/ref.fa.par.bwt",
        par_sa = "reference/ref.fa.par.sa",
    output:
        result = "analysis/{sample}.sorted.markdup.bam",
        disc = "analysis/{sample}.disc.sam",
        split = "analysis/{sample}.split.sam",
        interleave = "analysis/{sample}.interleave.fastq",
    params:
        lib_type = config["lib_type"]
        #fq2 = expand("raw_data/{fq2}", fq2=config["fq2"]),
    threads:
        config["threads"]
    log:
        biscuit= "logs/biscuit.{sample}.log",
        samblaster= "logs/samblaster.{sample}.log",
        view= "logs/samtools_view.{sample}.log",
        sort= "logs/samtools_sort.{sample}.log",
    singularity:
        "docker://dpettinga/biscuit_blaster"
    wrapper:
        "file:wrappers/biscuit_blaster"

#-------------------------------------------------------------------------------

rule post_process:
    input:
        result= "analysis/{sample}.sorted.markdup.bam",
        disc= "analysis/{sample}.disc.sam",
        split= "analysis/{sample}.split.sam",
        interleave= "analysis/{sample}.interleave.fastq",
    params:
        disc= "analysis/{sample}.disc.bam"
    output:
        result= "analysis/{sample}.sorted.markdup.bam.bai",
        disc="analysis/{sample}.disc.bam",
        disc_bai= "analysis/{sample}.disc.bam.bai",
        split= "analysis/{sample}.split.bam",
        split_bai= "analysis/{sample}.split.bam.bai",
        interleave= "analysis/{sample}.interleave.fastq.gz",
    threads:
        config["threads"]
    singularity:
        "docker://dpettinga/biscuit_blaster"
    shell:
        """
        samtools sort -o {output.disc} -O BAM -@ {threads} {input.disc}
        samtools index {output.disc}
        samtools sort -o {output.split} -O BAM -@ {threads} {input.split}
        samtools index {output.split}
        pigz -p {threads} {input.interleave}
        samtools index {input.result}
        """

rule flagstat_result:
    input:
        "analysis/{sample}.sorted.markdup.bam",
    output:
        "analysis/{sample}.sorted.markdup.bam.flagstat",
    singularity:
        "docker://dpettinga/biscuit_blaster"
    shell:
        "samtools flagstat {input} > {output}"
