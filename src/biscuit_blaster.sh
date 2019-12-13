#! /bin/bash

LIB_TYPE=${1?Error: No library type given}
SAMPLE_NAME=${2?Error: No sample name given}
REFERENCE=${3?Error: No reference.fa given}
THREADS=${4?Error: Number of threads not given}
FQ1=${5?Error: No fq1 given}
FQ2=${6:}

if [${FQ2} == ""]
then
  snakemake \
  -s Snakefile \
  --config \
  lib_type=${LIB_TYPE} \
  sample=${SAMPLE_NAME} \
  fq1=${FQ1} \
  reference=${REFERENCE} \
  threads=${THREADS} \
  --use-singularity \
  --singularity-args "-B /secondary,/primary"

elif [${FQ2} != ""]
then
  snakemake \
  -s Snakefile \
  --config \
  lib_type=${LIB_TYPE} \
  sample=${SAMPLE_NAME} \
  fq1=${FQ1} \
  fq2=${FQ2} \
  reference=${REFERENCE} \
  threads=${THREADS} \
  --use-singularity \
  --singularity-args "-B /secondary,/primary"

fi
