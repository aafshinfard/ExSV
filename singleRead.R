40096820

library("BSgenome.Ecoli.NCBI.20080805");
library("BSgenome.Hsapiens.UCSC.hg38");
library("BSgenome.Mmusculus.UCSC.mm10");
library("RSVSim")

Ref_Human = BSgenome.Hsapiens.UCSC.hg38$chr19;
#Ref_Human = BSgenome.Hsapiens.UCSC.hg38$chr19;
#Chr19_main = substr(as.character(Ref_Human),1,length(Ref_Human))

genome2 = readDNAStringSet( "/home/ameer/ExactSV/Chr19.fa" )
targ1 = readDNAStringSet( "/home/ameer/Desktop/SV_RefTarget/hg19MJR/Targ_Chr19_1.fasta" )
targ2 = readDNAStringSet( "/home/ameer/Desktop/SV_RefTarget/hg19MJR/Targ_Chr19_2.fasta" )
