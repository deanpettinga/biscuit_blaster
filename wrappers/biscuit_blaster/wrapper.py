__author__ = "Dean Pettinga"
__copyright__ = "Copyright 2019, Dean Pettinga"
__email__ = "dean.pettinga@vai.org"
__license__ = "MIT"

import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# get fq1
fq1 = snakemake.input.get("fq1")
assert fq1 is not None, "input-> fq1 is a required input parameter"


# get fq2 if it exists
if len(snakemake.params.get("fq2","")) != 0:
    fq2 = snakemake.params.get("fq2","")
else:
    fq2 = str("")

input_reads = " ".join([str(fq1), str(fq2)])

shell(
    "biscuit align -M "
    "-t {snakemake.threads} "
    "-b {snakemake.params.lib_type} "
    "{snakemake.input.ref} "
    "{input_reads} "
    "2> {snakemake.log.biscuit} | "
    "samblaster -M --excludeDups --addMateTags "
    "-d {snakemake.output.disc} "
    "-s {snakemake.output.split} "
    "-u {snakemake.output.interleave} "
    "2> {snakemake.log.samblaster} | "
    "samtools view -hbu -F 4 -q 20 2> {snakemake.log.view} | "
    "samtools sort -@ {snakemake.threads} -m 5G -o {snakemake.output.result} -O BAM - 2> {snakemake.log.sort}")
