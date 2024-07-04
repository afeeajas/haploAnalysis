setwd("")

## source the functions
source('haplotypeBglV5genogrm.R')

#make custom made haploblocks
haploblocks <- makehaploblocks(mapinfohap='GENOimp.map',nbsize=10000)

#or make haploblocks based on Haploview output
haploblocks <-getblockhaploview(haploviewfile="haploviewOutput")

## with phased beagle data, extract the haplotypes for each block
haplos.allele <- makehaplotypes(phasedbgl='GENOphased.vcf',mapinfohap='GENOimp.map',hapblocks=haploblocks)

