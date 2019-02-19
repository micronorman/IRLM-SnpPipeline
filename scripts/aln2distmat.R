#!/usr/bin/env Rscript

library(ape)

args<-commandArgs(trailingOnly = TRUE)

# Load alignment
aln_seq<-read.dna(args[1], format = "fasta")

# Calculate distance matrix
distmat<-dist.dna(aln_seq, model = "N")

write.csv(as.matrix(distmat),"")