#!/usr/bin/ruby

require 'utility'

module Sequence

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

	def Sequence.validateSam(samFile)
		
		
	end
	
	def Sequence.convertSamToBam(samFile)


	end

	def Sequence.sortBam(bamFile)


	end

	def Sequence.indexBam(bamFile)


	end

	def Sequence.countPair(bamFile)


	end

	def Sequence.seqComparison(fasta1, fasta2)


	end

	def Sequence.checkSeqLength(fastaFile)


	end	

	def Sequence.findORF(fastaFile)


	end	

end	