#!/usr/bin/ruby

require 'utility'
require 'fileutils'
require 'avail'

class Querydb
	
	attr_accessor :filenames, :prefix, :step	

	def initialize(prefix)
    	 @filenames = Array.new
    	 @prefix = prefix
    	 @step=Time.new
 	end

	def useBlast(dir, partials)
		
		begin
			#############
			#Blast specific parameters
			evalue="0.01".to_f
			outfmt="6".to_i
			maxTarget="50000".to_i
			count=16
			############

			Avail.moveFile(Utility::APP_ROOT, dir, partials) if Utility.const_defined?(:APP_ROOT)
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Partials file moved to DB directory"

			currentPath = Utility.navigate(dir)

			@filenames = Utility.splitFasta(partials)
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Fasta file is fragmented in individual sequence files!"
			
			@filenames.each_with_index do |file, index|
				currentPath = Utility.navigate(dir)

				if File.exist?("#{file}.fasta")

					cmd = "blastn -query #{file}.fasta -db #{@prefix} -out #{file}.out -evalue #{evalue} -outfmt #{outfmt} -max_target_seqs #{maxTarget} -num_threads #{count}"
					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Executing blast search for #{file}"
					Avail.executeCmd(cmd)

					cmd = "awk '{print $2}' #{file}.out | sed '1d' > #{file}.list"
					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Processing hits list..."
					Avail.executeCmd(cmd)

					if Avail.check_blast_out("#{file}.list") == true
						puts "#################################"	
						puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  ERROR::#{file} BLAST search returned no hits!!"
						puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Aborting #{file}"
						puts "#################################"
						temps = ["#{file}.fasta","#{file}.out","#{file}.list"]
						Avail.expunge(temps)
						next
					else 
						count = `wc -l "#{file}.list"`.strip.split(' ')[0].to_i
						puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  #{count} reads were found related with #{file}"		
					end	

					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  moving temp files....!!"
					Utility.createAndMoveFiles(file, "fasta", "Data")
					Utility.createAndMoveFiles(file, "out", "Data")
					Utility.createAndMoveFiles(file, "list", "Data")	

					currentPath = Utility.navigate("Data")	
					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  #{currentPath}"
				
					#cmd = "seqtk subseq #{@prefix}_R1_trimmed.fastq ./#{file}/#{file}.list > ./#{file}/#{file}_filtered_R1.fastq"
		
					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Fetching R1 mates..."	
					cmd = "seqtk subseq #{@prefix}_R1_trimmed.fastq ./#{file}/#{file}.list > ./#{file}/#{file}_filtered_R1.fastq"
					Avail.executeCmd(cmd)

					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Generated #{file}_filtered_R1.fasta" if File.exist? "#{file}/#{file}_filtered_R1.fasta"

					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Fetching R2 mates..."
					cmd = "seqtk subseq #{@prefix}_R2_trimmed.fastq ./#{file}/#{file}.list > ./#{file}/#{file}_filtered_R2.fasta"	
					Avail.executeCmd(cmd)

					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Generated #{file}_filtered_R2.fasta" if File.exist? "#{file}/#{file}_filtered_R2.fasta"
				
				else
					
					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  ERROR::#{file} related fasta file does not exist!"
					exit
				end	
			end	
		rescue => err
			$stderr.puts "Exception: #{err}\n\n"
			err.backtrace.each { |l| $stderr.puts l + "\n" }
			err	
		end		
	end	

	def useNhmmer(dir, partials)

		begin
			##################
			#nhmmer specific parameters			
			evalue="0.01".to_f
			count="16".to_i
			#################

			Avail.moveFile(Utility::APP_ROOT, dir, partials) if Utility.const_defined?(:APP_ROOT)
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Partials file moved to DB directory"

			currentPath = Utility.navigate(dir)
			fileNames = Utility.splitFasta(partials)
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Fasta file is fragmented in individual sequence files!"

			#fileNames has all the names of fasta files. This moethod iterates over this array and perform nhmmer.
			fileNames.each_with_index do |file, index|

				if File.exist?("#{file}.fasta")
					
					cmd = "nhmmer --tblout #{file}.out --noali --incE #{evalue} -E #{evalue} --max --dna --cpu #{count} #{file}.fasta #{@prefix}.fasta"
					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Executing nhmmer search for #{file}"
					Avail.executeCmd(cmd)

					cmd = "awk '{print $2}' #{file}.out | sed '1d' > #{file}.list"
					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Processing hits list..."
					Avail.executeCmd(cmd)

					#count = 0
					#File.open(filename) {|f| count = f.read.count("\n")}	
					if Avail.check_blast_out("#{file}.list") == true
						puts "#################################"	
						puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  ERROR::#{file} NHMMER search returned no hits!!"
						puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Aborting #{file}"
						puts "#################################"
						temps = ["#{file}.fasta","#{file}.out","#{file}.list"]
						Avail.expunge(temps)
						next
					else 
						count = `wc -l "#{file}.list"`.strip.split(' ')[0].to_i
						puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  #{count} reads were found related with #{file}"		
					end
	                #currentPath = Utility.navigate("Data")
	                #puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  #{currentPath}"
						
					cmd = "seqtk subseq #{@prefix}.fasta #{file}.list > #{file}.filtered_R1.fa"
					puts = "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Fetching R1 mates..."
					Avail.executeCmd(cmd)

					cmd = "seqtk subseq #{@prefix}.fasta #{file}.list > #{file}.filtered_R2.fa"
					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Fetching R2 mates..."	
					Avail.executeCmd(cmd)

					##Need to move all the files in temp/Data directory. Files ending with *.list *.fasta and store it under file folder.
					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  moving temp files....!!"
	                Utility.createAndMoveFiles(file, "fasta", "Data")
	                Utility.createAndMoveFiles(file, "out", "Data")
	                Utility.createAndMoveFiles(file, "list", "Data")
	                Utility.createAndMoveFiles(file, "filtered_R1.fa", "Data")
	                Utility.createAndMoveFiles(file, "filtered_R2.fa", "Data")

	                currentPath = Utility.navigate("Data")

				else
					puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  ERROR::#{file} related fasta file does not exist!"
					exit
				end	
			end
		rescue => err
			$stderr.puts "Exception: #{err}\n\n"
			err.backtrace.each { |l| $stderr.puts l + "\n" }
			err	
		end	
	end	
end	
