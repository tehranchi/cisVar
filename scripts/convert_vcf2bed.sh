#!/bin/bash
# Convert vcf to bed, only including SNVs

module load bedops
in=$1
out=$(echo "${1}" | sed 's/vcf.comm.gz/comm.bed.gz/')
# Convert to bed and drop indels
pigz -dc ${in} | vcf2bed --snvs | pigz > ${out}
