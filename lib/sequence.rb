#!/usr/bin/ruby

require 'utility'
require 'avail'

module Sequence

	@step=Time.new

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
end	
