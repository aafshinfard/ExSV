#####
library("RSVSim")
genome = DNAStringSet(
  c("AAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTT",
    "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCC"))
names(genome) = c("chr1","chr2")
genome

simulateSV() #hg19 by default

########## Deletion
sim = simulateSV(output=NA, genome=genome, dels=3, sizeDels=10,
                 bpSeqSize=6, verbose=FALSE) # seed = const -> same simulations
sim
metadata(sim)

########## Insertion (cut only)
sim2 = simulateSV(output=NA, genome=genome, ins=3, sizeIns=5, bpSeqSize=6,
                 seed=246, verbose=FALSE)
sim2
metadata(sim2)

########## Insertion (Copy and cut)
sim3 = simulateSV(output=NA, genome=genome, ins=3, sizeIns=5, percCopiedIns=0.66,
                 bpSeqSize=6, seed=246, verbose=FALSE)
sim3
metadata(sim3)

########## Inversions (with different size) -> invert in place
sim4 = simulateSV(output=NA, genome=genome, invs=3, sizeInvs=c(2,4,6),
                 bpSeqSize=6, seed=456, verbose=FALSE)
sim4
metadata(sim4)

########## Tandem Duplication
sim5 = simulateSV(output=NA, genome=genome, dups=1, sizeDups=6, maxDups=3,
                 bpSeqSize=6, seed=3456, verbose=FALSE)
sim5
metadata(sim5)

########## Translocation (from 5' or 3' ends - balanced or unbalanced -> percBalancedTrans (balancing ratio))
sim6 = simulateSV(output=NA, genome=genome,trans=1, bpSeqSize=6, seed=123, verbose=FALSE)
sim5
metadata(sim5)