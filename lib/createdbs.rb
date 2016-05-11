
require 'utility'

class Createdbs

	def formatReadFiles(file)
		#convert from zipped to fastq
		#fastq to fasta
		begin
			if file.is_a?(Hash)
				file.each do |value|
					convertQ2A(Utility.unzip(value))
				end	
			else
				Puts "Input file is not paired!"
			end	
		rescue Exception => e
			e.backtrace
			puts "Oopss! Faulty file conversion. Check the raw read files!"
		end	
	end	

	def useBlast



	end	

	def useNhmmer


	end

	def queryReads


	end

	def fetchPairs


	end	




end	
