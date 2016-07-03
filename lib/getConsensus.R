#Consensus Generator : Uses BioStrings and Decipher package.
#Takes the fasta file generated in extension files from SMALT generated sam/bam file 
#Generates consensus sequence file from a set of sequences present in a fasta file.
#Author: Singh KS
#Date: 3rd July 2016

# Function to check whether package is installed
is.installed <- function(x){
  is.element(x, installed.packages()[,1])
}

checkInstallation <- function(x) {
  # check if package x is installed; If not install it
  if (!is.installed(x)) {
    cat('Package ',x,' is not found!', '\n')
    cat('Installing ',x, '\n')
    install.packages(x)
  } else
    {
      cat('Package ',x,' is found!', '\n')  
    }
}

#stop("Package not found")
#get the consensus of sequences present in fasta file
getConsensus <- function(file, out) {
  #
  checkInstallation("Biostrings")
  checkInstallation("DECIPHER")
  #
  readFasta <- readDNAStringSet(file)
  msa <- AlignSeqs(dna)
  consensus <- ConsensusSequence(msa)
  #
  writeXStringSet(consensus, file=out)

}
library("Biostrings")
library("Biostrings") 
getConsensus("left.fasta", "left.consensus")
getConsensus("right.fasta", "right.consensus")