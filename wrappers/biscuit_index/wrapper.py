__author__ = "Dean Pettinga"
__copyright__ = "Copyright 2019, Dean Pettinga"
__email__ = "dean.pettinga@vai.org"
__license__ = "MIT"

import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

assert (len(snakemake.params) == 1), "Only 1 param allowed:'DEFAULT' or 'path/to/reference.fa'. Please set reference in src/config.yaml"
assert snakemake.params[0] is not None, "Only 1 param allowed:'DEFAULT' or 'path/to/reference.fa'. Please set reference in src/config.yaml"

if snakemake.params[0] == "DEFAULT_HUMAN":
    if not os.path.isfile("reference/GRCh38.primary_assembly.genome.fa"):
        shell("wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz 2> {snakemake.log[0]}| gunzip -c > reference/GRCh38.primary_assembly.genome.fa ")
    else:
        shell("ln -sr reference/GRCh38.primary_assembly.genome.fa reference/ref.fa")
        shell("biscuit index reference/ref.fa 2> {snakemake.log[1]}")

elif snakemake.params[0] == "DEFAULT_MOUSE":
    if not os.path.isfile("reference/GRCm38.primary_assembly.genome.fa"):
        shell("wget -O -ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz | gunzip -c > reference/GRCm38.primary_assembly.genome.fa ")
    else:
        shell("ln -sr reference/GRCm38.primary_assembly.genome.fa reference/ref.fa")
        shell("biscuit index reference/ref.fa 2> {snakemake.log}")
else:
    shell("ln -s {snakemake.params[0]} reference/ref.fa")
    shell("ln -s {snakemake.params[0]}.bis.amb reference/ref.fa.bis.amb")
    shell("ln -s {snakemake.params[0]}.bis.ann reference/ref.fa.bis.ann")
    shell("ln -s {snakemake.params[0]}.bis.pac reference/ref.fa.bis.pac")
    shell("ln -s {snakemake.params[0]}.dau.bwt reference/ref.fa.dau.bwt")
    shell("ln -s {snakemake.params[0]}.dau.sa reference/ref.fa.dau.sa")
    shell("ln -s {snakemake.params[0]}.par.bwt reference/ref.fa.par.bwt")
    shell("ln -s {snakemake.params[0]}.par.sa reference/ref.fa.par.sa")


# if snakemake.params[0] == "DEFAULT":
#     shell_1 = "wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz | gunzip -c > reference/GRCh38.primary_assembly.genome.fa "
#     shell_2 = "ln -s reference/GRCh38.primary_assembly.genome.fa reference/ref.fa"
#     shell_3 = "biscuit index -p reference/ref.fa 2> {log}"
#     shell_4 = ""
#     shell_5 = ""
#     shell_6 = ""
#     shell_7 = ""
#     shell_8 = ""
#     shell_9 = ""
# else:
#     shell_1 = "ln -s {snakemake.input[0]} reference/ref.fa"
#     shell_2 = "ln -s {snakemake.input[0]}.bis.amb reference/ref.fa.bis.amb"
#     shell_3 = "ln -s {snakemake.input[0]}.bis.ann reference/ref.fa.bis.ann"
#     shell_4 = "ln -s {snakemake.input[0]}.bis.pac reference/ref.fa.bis.pac"
#     shell_5 = "ln -s {snakemake.input[0]}.dau.bwt reference/ref.fa.dau.bwt"
#     shell_6 = "ln -s {snakemake.input[0]}.dau.sa reference/ref.fa.dau.sa"
#     shell_7 = "ln -s {snakemake.input[0]}.par.bwt reference/ref.fa.par.bwt"
#     shell_8 = "ln -s {snakemake.input[0]}.par.sa reference/ref.fa.par.sa"
#     shell_9 = "biscuit index -p reference/ref.fa 2> {log}"
#
# shell("{shell_1}"
#     "{shell_2}"
#     "{shell_3}"
#     "{shell_4}"
#     "{shell_5}"
#     "{shell_6}"
#     "{shell_7}"
#     "{shell_8}"
#     "{shell_9}")
