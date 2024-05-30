setwd("")

## source the functions
source('haplotypeBglV5genogrm.R')

#make custom made haploblocks
haploblocks <- makehaploblocks()

## after phasing your data with beagle extract the haplotypes for each block
haplos.allele <- makehaplotypes(phasedbgl='/home/afeesa/paper3/data/haploview_dat/GENOimphased_ssa12.vcf.gz.recode.vcf',mapinfohap='/home/afeesa/paper3/data/haploview_dat/GENOimp_filtssa12.map',hapblocks=haploblocks)

