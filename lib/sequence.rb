#!/usr/bin/ruby

require 'utility'
require 'avail'

module Sequence

	class DataFormatError < IOError
    	def self.message
      		puts "Data format error -- check input file"
      		exit
    	end
  	end


	def Sequence.checkSequence(seq)
		
		count = 0
		
		if seq.include?("N")
			
			count = seq.count("N")
			puts "WARNING:: #{count} N characters are found in your partial. It will have serious implication on the quality of mapping. Please remove it and try again!"
			exit
		
		else
			puts "Partial is OK!"			
		end	

	end

	def Sequence.countFastaSeqs(file)
		raise DataFormatError.message unless Utility.fileType(file) == "fasta"
		count = 0
		File.open(file.to_s, "r") do |file|
			file.each_line do |line|
				if line.start_with?(">")
					count += 1
				else
					next
				end	
			end	
		end	

		return count

	end	

	def Sequence.countHits(blastOutFile)
		
		Utility.navigate("DB")

		if Utility.checkFileExist(blastOutFile)

			cmd = "awk '{print $2}' #{blastOutFile} | wc -l"
			Utility.executeCmd(cmd)
		
		else
			
			puts "ERROR::Blast output file is not found!"
			exit
		
		end	

	end	

	def Sequence.validateBam(seqHeader)
		cmd="samtools quickcheck -v #{seqHeader}.bam && echo 'BAM looks OK' || echo 'WARNING::BAM file failed check!!'"
		Utility.executeCmd(cmd)
	end

	def Sequence.createIndex(seqHeader)
		cmd="smalt index -k 5 -s 1 #{seqHeader} #{seqHeader}.fasta"
		Utility.executeCmd(cmd)
	end	

	def Sequence.runMapper(seqHeader, fastq=[])
		cmd="smalt map -f samsoft:x -l pe -o #{seqHeader}.sam -p -x #{seqHeader} #{fastq[0]} #{fastq[1]}"
		Utility.executeCmd(cmd)
	end	
	
	def Sequence.convertSamToBam(seqHeader)
		cmd="samtools view -b -S -o #{seqHeader}.bam #{seqHeader}.sam"
		Utility.executeCmd(cmd)
	end

	def Sequence.sortBam(seqHeader)
		cmd="samtools sort -o #{seqHeader}.sorted.bam #{seqHeader}.bam"
		Utility.executeCmd(cmd)
	end

	def Sequence.indexBam(seqHeader)
		cmd="samtools index #{seqHeader}.sorted.bam"
		Utility.executeCmd(cmd)
	end

	def Sequence.createConsensus(seqHeader)
		cmd="samtools mpileup -vf #{seqHeader}.fasta #{seqHeader}.sorted.bam | bcftools call -m -O z - > #{seqHeader}.vcf.gz"
		Utility.executeCmd(cmd)

		cmd="bcftools index #{seqHeader}.vcf.gz"
		Utility.executeCmd(cmd)

		cmd="bcftools consensus -i -f #{seqHeader}.fasta #{seqHeader}.vcf.gz > #{seqHeader}_1st_consensus.fasta"
		Utility.executeCmd(cmd)
	end	

	def Sequence.getEndsOverlapingReads(seqHeader, pass=nil)

		fileObject=""
		filename=""
			if pass.eql?(1)
				filename = "#{seqHeader}_1st_consensus.fasta"
			elsif pass.eql?(2)
				filename = "#{seqHeader}_2nd_consensus.fasta"
			elsif pass.eql?(3)
				filename = "#{seqHeader}_3rd_consensus.fasta"	
			else
				filename = "#{seqHeader}.fasta"
			end			

		seqProp = Utility.getSequenceLength(filename)

		len = seqProp[seqHeader].to_i
			1.times do |i|
				fileObject = "#{seqHeader}\t1\t20\n"
				Avail.writeFile("left.bed", fileObject, headerPrefix = nil)
			end	
			1.times do |i|
				innerValue = len-20
				fileObject = "#{seqHeader}\t#{innerValue}\t#{len}\n"
				Avail.writeFile("right.bed", fileObject, headerPrefix = nil)
			end	

		cmd="samtools view -L left.bed TRINITY_DN10_c0_g1_i1_1st_consensus_sorted.bam | awk '/^@/{ next; } { print ">"$1; print $10 }' - > left.fasta"
		Utility.executeCmd(cmd)	

		cmd="samtools view -L right.bed TRINITY_DN10_c0_g1_i1_1st_consensus_sorted.bam | awk '/^@/{ next; } { print ">"$1; print $10 }' - > right.fasta"
		Utility.executeCmd(cmd)

	end	

	def Sequence.getConsensus

		File.open("getConsensus.R", "w") do |file|
	
			file << "is.installed <- function(x){\n"
		  	file << "is.element(x, installed.packages()[,1])\n"
			file << "}\n"
			file << "checkInstallation <- function(x) {\n"
		  	file << "if (!is.installed(x)) {\n"
		    file << "cat(\'Package \',x,\' is not found!\', \'\n\')\n"
		    file << "cat(\'Installing \',x, \'\n\')\n"
		    file << "install.packages(x)\n"
		  	file << "} else\n"
		    file << "{\n"
		    file << "cat(\'Package \',x,\' is found!\', \'\n\')\n"  
		    file << "}\n"
			file << "}\n"
			file << "getConsensus <- function(file, out) {\n"
			file << "checkInstallation(\"Biostrings\")\n"
			file << "checkInstallation(\"DECIPHER\")\n"
			file << "readFasta <- readDNAStringSet(file)\n"
			file << "msa <- AlignSeqs(dna)\n"
			file << "consensus <- ConsensusSequence(msa)\n"
			file << "writeXStringSet(consensus, file=out)\n"
			file << "}\n"
			file << "library(\"Biostrings\")"
			file << "library(\"Biostrings\")"
			file << "getConsensus(\"left.fasta\", \"left.consensus\")\n"
			file << "getConsensus(\"right.fasta\", \"right.consensus\")\n"
		end

		Utility.runRscript("getConsensus.R")

	end	


end	