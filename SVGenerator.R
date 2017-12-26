library("BSgenome.Ecoli.NCBI.20080805");
library("BSgenome.Hsapiens.UCSC.hg38");
library("BSgenome.Mmusculus.UCSC.mm10");
library("RSVSim")

Ref_Human = BSgenome.Hsapiens.UCSC.hg38$chr19;
#Ref_Human = BSgenome.Hsapiens.UCSC.hg38$chr19;
#Chr19_main = substr(as.character(Ref_Human),1,length(Ref_Human))

genome2 = readDNAStringSet( "/home/ameer/ExactSV/Chr19.fa" )

class(Ref_Human)
class(genome2$chr19)
Ref_Human = genome2$chr19

Ref_Genome = DNAStringSet(Ref_Human)
names(Ref_Genome) = c("Hchr19")
Ref_Genome
writeXStringSet(Ref_Genome, "~/Desktop/SV_RefTarget/Data/Ref_chr19.fasta",format="fasta");

############################################################
############################## HOMO SV

sv_main = simulateSV(output=NA, genome=Ref_Genome, dels=2, invs = 2, ins=2, size = c(100,400), bpSeqSize=9, seed=100, verbose=FALSE) # seed = const -> same simulations
metadata(sv_main)


############################################################
############################## HETERO SV

sv1 = simulateSV(output=NA, genome=sv_main, dels=1, invs = 1, ins=1, size = c(150,350), bpSeqSize=9, seed=102, verbose=FALSE) # seed = const -> same simulations
metadata(sv1)
#sv2 = simulateSV(output=NA, genome=sv_main, dels=3, invs = 3, ins=3, size = c(80,160,350), bpSeqSize=9, seed=104, verbose=FALSE) # seed = const -> same simulations
sv2 = sv_main
metadata(sv2)
writeXStringSet(sv1, "~/Desktop/SV_RefTarget/Data/Targ_1.fasta",format="fasta");
writeXStringSet(sv2, "~/Desktop/SV_RefTarget/Data/Targ_2.fasta",format="fasta");



rsv1 = sv1$Hchr19
rsv2 = sv2$Hchr19
strsv1 = substr(as.character(rsv1),1,length(rsv1))
strsv2 = substr(as.character(rsv2),1,length(rsv2))

########## generating SNPs in Target Genomes
#Targ_Genome = Ref_Genome;
G1 = sum(length(rsv1))
G2 = sum(length(rsv2))
SNP_rate = .001

p = SNP_rate/(SNP_rate + 1) 
total1 = floor(G1*SNP_rate)
total2 = floor(G2*SNP_rate)

SNP_intervals1 = rgeom(2*total1,p)+1
SNP_intervals2 = rgeom(2*total2,p)+1
SNP_cum1 = cumsum(SNP_intervals1)
SNP_cum2 = cumsum(SNP_intervals2)
SNP_locs1 = SNP_cum1[SNP_cum1 < G1]
SNP_locs2 = SNP_cum2[SNP_cum2 < G2]

### Creating the target Genome
nucleotide2numbers = function (Nucleotide){
  switch(Nucleotide,A = 0, C=1,G=2, T= 3,N=4)
}

number2nucleotide = function (number){
  switch(number+1,"A", "C","G", "T","N")
}

mutations = function ( Ref_SNP, substitution_matrix = matrix()){
  seq = unlist(strsplit(as.character(Ref_SNP),split=""))
  print(dim(substitution_matrix))
  if (sum(dim(substitution_matrix)) == 2){
    t = length(Ref_SNP)
    r_mutations = sample(x=c(1,2,3),size=t,replace=T)
    x = sapply(X=seq,FUN=nucleotide2numbers)
    y = (x+r_mutations) %% 4
    y[x==4]=4
    SNPs = sapply(X=y,FUN=number2nucleotide) 
    return(SNPs)
    
  }
}
chr1_snps = rsv1[SNP_locs1]
chr2_snps = rsv2[SNP_locs1]
SNPs1 = mutations(chr1_snps)
SNPs2 = mutations(chr2_snps)
Tar_chr1_SNP = DNAString(paste(SNPs1,collapse=""))
Tar_chr2_SNP = DNAString(paste(SNPs2,collapse=""))
Tar_chr1 = rsv1
Tar_chr2 = rsv2
Tar_chr1[SNP_locs1] = Tar_chr1_SNP
Tar_chr2[SNP_locs1] = Tar_chr2_SNP

Targ_Chr19_1 = DNAStringSet(x = Tar_chr1)
Targ_Chr19_2 = DNAStringSet(x = Tar_chr2)
names(Targ_Chr19_1) = c("Chr19_1")
names(Targ_Chr19_2) = c("Chr19_2")
#sum(Tar_chr1 == chr1)

class(Targ_Chr19_1)
writeXStringSet(Targ_Chr19_1, "~/Rcodes2016/Files/Targ_Chr19_1.fasta",format="fasta");
writeXStringSet(Targ_Chr19_2, "~/Rcodes2016/Files/Targ_Chr19_2.fasta",format="fasta");
metadata(sv_main)
metadata(sv1)
metadata(sv2)

######################### MORE SV Sim

simulateSV()
Targ_Genome2 = simulateSV(output=NA, genome=Ref_Genome, dels=2,invs = 2,ins=2,size = c(70,400),bpSeqSize=9,seed=100, verbose=FALSE) # seed = const -> same simulations


# novel insertion

#### RSVSim structural variation simulation 

#Normal SV simulation:
Targ_Genome2 = simulateSV(output=NA, genome=Targ_Genome, dels=2,invs = 2,ins=2,size = c(70,400),bpSeqSize=9,seed=100, verbose=FALSE) # seed = const -> same simulations
metadata(Targ_Genome2)

Tar_chr1 = Targ_Genome2$chr19SV
novel_insertion = substr(as.character(Ref_Ecoli),200001,200500)
Tar_chr1 = c(Tar_chr1[1:90000],DNAString(novel_insertion),Tar_chr1[90001:length(Tar_chr1)])
Targ_Genome3 = DNAStringSet(x = Tar_chr1)
metadata(Targ_Genome3)=metadata(Targ_Genome2)
a = metadata(Targ_Genome3)$inversions
a = a[1,]
a[1,1] = "novel_insertion";
a = metadata(Targ_Genome3)$inversions
a = a[1,]
a[1,1] = "novel_insertion";
a[1,3] = 90001;
a[1,4] = 90500;
a[1,5] = 500;
a=a[1,1:5];
names(Targ_Genome3)="chr19SV";
metadata(Targ_Genome3)$novel_insertions = a;
metadata(Targ_Genome3)
writeXStringSet(Targ_Genome3, "~/Rcodes2016/Files/Target_SV_chr19.fasta",format="fasta");

# Complex SV simulation:
regioninv = GRanges(c( IRanges(15174256,15174756)),seqnames="chr19SV")
Targ_Genome2 = simulateSV(output=NA, genome=Targ_Genome,invs = 1,regionsInvs =regioninv ,size = c(500),
                          bpSeqSize=9, verbose=FALSE) # seed = const -> same simulations
regioninv = GRanges(c( IRanges(16100553,16100852)),seqnames="chr19SV")
Targ_Genome2 = simulateSV(output=NA, genome=Targ_Genome2,invs = 1,regionsInvs =regioninv ,size = c(300),
                          bpSeqSize=9, verbose=FALSE) # seed = const -> same simulations
regioninv = GRanges(c( IRanges(19634731,19634800)),seqnames="chr19SV")
Targ_Genome2 = simulateSV(output=NA, genome=Targ_Genome2,invs = 1,regionsInvs =regioninv ,size = c(70),
                          bpSeqSize=9, verbose=FALSE) # seed = const -> same simulations
Targ_Genome2 -> ttg

regiondel = GRanges(c(IRanges(40215321,40215390) ),seqnames="chr19SV")
Targ_Genome3 = simulateSV(output=NA, genome=Targ_Genome2,dels = 1,regionsDels =regiondel ,size = c(70),
                          bpSeqSize=9, verbose=FALSE) # seed = const -> same simulations
regiondel = GRanges(c(IRanges(16100353,16100653) ),seqnames="chr19SV")
Targ_Genome3 = simulateSV(output=NA, genome=Targ_Genome3,dels = 1,regionsDels =regiondel ,size = c(300),
                          bpSeqSize=9, verbose=FALSE) # seed = const -> same simulations

regiondel = GRanges(c(IRanges(31068081,31068581) ),seqnames="chr19SV")
Targ_Genome3 = simulateSV(output=NA, genome=Targ_Genome3,dels = 1,regionsDels =regiondel ,size = c(500),
                          bpSeqSize=9, verbose=FALSE) # seed = const -> same simulations





regionins = GRanges(IRanges(15174456,15174956),seqnames="chr19SV",chrB="chr19SV",startB=3068381)
names(regionins)="KnownIns"
Targ_Genome4 = simulateSV(output=NA, genome=Targ_Genome3,regionsIns = regionins ,sizeIns = 500,
                          bpSeqSize=9, random=FALSE) # seed = const -> same simulations





#regionins = GRanges(c(IRanges(16100553,16100852)),seqnames="chr19SV",chrB="chr19SV",startB=c(2850010))
#Targ_Genome2 = simulateSV(output=NA, genome=Targ_Genome2,invs = 1,regionsInvs =IRanges(15174256,15174756) ,size = c(300),
                          #bpSeqSize=9, verbose=FALSE) # seed = const -> same simulations
regionins = GRanges(IRanges(25634631,25634700),seqnames="chr19SV",chrB="chr19SV",startB=33358025)
names(regionins)="KnownIns"
Targ_Genome4 = simulateSV(output=NA, genome=Targ_Genome4,regionsIns = regionins ,sizeIns = 70,
                          bpSeqSize=9, random=FALSE) # seed = const -> same simulations

#regionins = GRanges(c(IRanges(25634400,25634899)),seqnames="chr19SV",chrB="chr19SV",startB=c(22358025))




#regionins = GRanges(c(IRanges(15174456,15174956)),seqnames="chr19SV",chrB="chr19SV",startB=c(30168381))

Targ_Genome2 = simulateSV(output=NA, genome=Targ_Genome,invs = 1,regionsInvs =IRanges(15174256,15174756) ,size = c(500),
                          bpSeqSize=9, verbose=FALSE) # seed = const -> same simulations


Targ_Genome2 = simulateSV(output=NA, genome=Targ_Genome, dels=3,regionsDels =regiondel ,invs = 3,regionsInvs =regioninv ,ins=3,regionsIns = regionins,size = c(500,300,70),
                          bpSeqSize=9, verbose=FALSE) # seed = const -> same simulations
Targ_Genome2
metadata(Targ_Genome4)

#Targ_Genome2 = simulateSV(output=NA, genome=Targ_Genome, dels=2,invs = 2,ins=2,size = c(70,400),
#                          bpSeqSize=9,seed=101, verbose=FALSE) # seed = const -> same simulations
#Targ_Genome2
#metadata(Targ_Genome2)
#Targ_Genome3 = simulateSV(output=NA, genome=Targ_Genome2, dels=1,invs = 1,ins=1,size = c(150,500),
#                          bpSeqSize=9,seed=100, verbose=FALSE) # seed = const -> same simulations
#Targ_Genome3
#metadata(Targ_Genome3)
Targ_Genome2 = Targ_Genome4

Tar_chr1 = Targ_Genome2$chr19SV
novel_insertion = substr(as.character(Ref_Ecoli),200001,200500)
Tar_chr1 = c(Tar_chr1[1:90000],DNAString(novel_insertion),Tar_chr1[90001:length(Tar_chr1)])
Targ_Genome3 = DNAStringSet(x = Tar_chr1)
names(Targ_Genome) = c("chr1")
names(Targ_Genome) = c("chr19")
metadata(Targ_Genome3)=metadata(Targ_Genome2)
a = metadata(Targ_Genome3)$inversions
a = a[1,]
a[1,1] = "novel_insertion";
a[1,3] = 90001;
a[1,4] = 90500;
a[1,5] = 500;
a=a[1,1:5];
metadata(Targ_Genome3)$novel_insertions = a;
writeXStringSet(Targ_Genome3, "~/Rcodes2016/Files/Target_compSV_chr19.fasta",format="fasta");




######
Tar_chr1 = Targ_Genome2$chr1
Tar_chr2 = Targ_Genome2$chr2

G1 = length(Tar_chr1)
G2 = length(Tar_chr1)



#### now we want to generate reads from this target genome

Depth = 30
L = 150

N1 = floor(G1*Depth/L)
N2 = floor(G2*Depth/L)

lambda1 = N1/G1
lambda2 = N2/G2

Read_intervals1 = rgeom(2*N1,lambda1)+1
Read_intervals2 = rgeom(2*N2,lambda2)+1
Read_cum1 = cumsum(Read_intervals1)
Read_cum2 = cumsum(Read_intervals2)
Read_locs1 = Read_cum1[Read_cum1 < G1-L+1]
Read_locs2 = Read_cum2[Read_cum2 < G2-L+1]


v1 = Views(Tar_chr1,start=Read_locs1,width=L)
v2 = Views(Tar_chr2,start=Read_locs2,width=L)

Reads_chr1 = toString(v1)
Reads_chr2 = toString(v2)

Reads_chr1 = strsplit(Reads_chr1,split=", ")
Reads_chr2 = strsplit(Reads_chr2,split=", ")
Reads_chr1 = unlist(Reads_chr1)
Reads_chr2 = unlist(Reads_chr2)

readFile = file("SVreads.fq","w")
QScore = paste(rep("I",L),collapse="")
for (i in 1:length(Reads_chr1)){
  writeLines(paste("@",Read_locs1[i],sep=""),con=readFile)
  writeLines(Reads_chr1[i],con=readFile)
  writeLines("+",con=readFile)
  writeLines(QScore,con=readFile)
}
for (i in 1:length(Reads_chr2)){
  writeLines(paste("@",Read_locs2[i]+length(Tar_chr1),sep=""),con=readFile)
  writeLines(Reads_chr2[i],con=readFile)
  writeLines("+",con=readFile)
  writeLines(QScore,con=readFile)
}

close(con=readFile)


########## 

# Read Generation:



# others :
# ########## Insertion (cut only)
# sim2 = simulateSV(output=NA, genome=genome, ins=3, sizeIns=5, bpSeqSize=6,
#                   seed=246, verbose=FALSE)
# sim2
# metadata(sim2)
# 
# ########## Insertion (Copy and cut)
# sim3 = simulateSV(output=NA, genome=genome, ins=3, sizeIns=5, percCopiedIns=0.66,
#                   bpSeqSize=6, seed=246, verbose=FALSE)
# sim3
# metadata(sim3)
# 
# ########## Inversions (with different size) -> invert in place
# sim4 = simulateSV(output=NA, genome=genome, invs=3, sizeInvs=c(2,4,6),
#                   bpSeqSize=6, seed=456, verbose=FALSE)
# sim4
# metadata(sim4)
# 
# ########## Tandem Duplication
# sim5 = simulateSV(output=NA, genome=genome, dups=1, sizeDups=6, maxDups=3,
#                   bpSeqSize=6, seed=3456, verbose=FALSE)
# sim5
# metadata(sim5)
# 
# ########## Translocation (from 5' or 3' ends - balanced or unbalanced -> percBalancedTrans (balancing ratio))
# sim6 = simulateSV(output=NA, genome=genome,trans=1, bpSeqSize=6, seed=123, verbose=FALSE)
# sim5
# metadata(sim5)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 






