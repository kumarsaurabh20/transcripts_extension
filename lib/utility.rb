#!/usr/bin/ruby

require 'prerequisite'
require 'pathname'
require 'fileutils'
require 'avail'

module Utility

	APP_ROOT = File.expand_path("../../", __FILE__) 
	@step=Time.new
		

	class BadFileError < StandardError
    end

	# Check and set path to all the required tools and software
	#Use prerequisite shell script
	# find out if some programmes are not installed.
	def Utility.checkAndSetTools
		check = Prerequisite.new
		check.checkTools
		exit
	end		

	def Utility.createDbFasta(prefix, array)
		
		Avail.makeDir("DB")		
		
		if array.is_a?(Array)
			
			cmd = "cat #{array[0]} #{array[1]} > #{prefix}.fasta"
			
			Avail.executeCmd(cmd) 
		else 
			
			raise "Expected Array Value to create a final fasta file!"	
		end
	end

	# check the write permission of $workDir before building of the work directory
	def Utility.checkPermissions(file)
		
		if File.exist?(file) && File.executable?(file)

			return true
		else
			return false
		end	
	end

	def Utility.fileType(file)	
		
		ex = File.extname(file)
		
		if ex.eql?(".fasta") || ex.eql?(".fa") || ex.eql?(".fsa")
			
			return "fasta"
		
		elsif ex.eql?(".gz")
			
			return "zipped"
		
		elsif ex.eql?(".fastq")
			
			return "fastq"
		else
			
			return "other"
		end
	end

	def Utility.navigate(folder)
		
		temp = ""
		setpath = Avail.setDir
		FileUtils.cd(setpath)

		if File.directory?(folder)
			temp = File.join(setpath, folder)
			FileUtils.cd(temp)
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Current working directory is moved to #{folder}"
		else
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  ERROR::Data folder is not found!"
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  If Data folder is not available, create a folder and name it Data and dump all your raw data files!"
			exit
		end

		return temp
	end

	def Utility.mergeZippedFiles(file1, file2, prefix)
		r1 = File.expand_path(file1)
		puts r1.to_s
		r2 = File.expand_path(file2)
		puts r2.to_s
		name = File.basename(file1, "*.fastq.gz")
		outfile = "#{prefix}_merged.fastq"
		cmd = "cat #{r1} #{r2} > #{outfile}"
		puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Merging zipped fastq files!"
		Avail.executeCmd(cmd)
		puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  This is the ourput of mergeZipped method #{outfile}"
		return outfile
	end	

	def Utility.unzip(file)		
		#Avail.navigate(dir)
		cmd = "gzip -d #{file}"
		Avail.executeCmd(cmd)
	end	

	# check whether fasta file exist and how many sequences it contains
	def Utility.countSeq(file)
		cmd = "grep -c '^>' #{file}"
		count = `#{cmd}`		
		return count
	end

	def Utility.convertQ2A(file, prefix)
		name = File.basename(file, ".*")
		outfile = "#{name}.fasta"
		if name.empty?
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Error::FastQ File not found!"
		else
			if defined? prefix
				cmd = "seqtk seq -A #{file} > #{prefix}.fasta"
				Avail.executeCmd(cmd)			
				return "#{prefix}.fasta"
			else
				cmd = "seqtk seq -A #{file} > #{outfile}"
				Avail.executeCmd(cmd)	
				return "#{outfile}.fasta"
			end
		end	
	end	

	def Utility.directory_exists?(dir)
		
		Avail.setDir
		
		if File.directory?(dir)
			Avail.navigate(dir)
		else
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  #{dir} does not exist. Check again!"
		end	
	end

	def Utility.checkFileExist(file)

		if File.exist?(file) and File.executable?(file)
			return true
		else
			return false
		end		

	end	


	def Utility.createAndMoveFiles(file, ext, dir)
		
		#fileBasename = File.basename(File.expand_path(file), ".*")
		dirname = File.join("#{APP_ROOT}", dir, file)

		if Dir.exist? dirname
			temp = File.expand_path("#{file}.#{ext}", "#{APP_ROOT}/DB")
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Transferring file #{temp}"
			FileUtils.mv temp, dirname, :force => true
			
		else
			Avail.createDir(dir, file)
			temp = File.expand_path("#{file}.#{ext}", "#{APP_ROOT}/DB")
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Transferring file #{temp}"
			FileUtils.mv temp, dirname, :force => true
		end	
		
	end

	def Utility.writeFile(outputFile, seqObject, headerPrefix)
		
		File.open(outputFile.to_s, 'a') do |file|
		
			seqObject.each_with_index do |element, index|
			
				file.write ">#{headerPrefix}#{index}\n#{element}\n"
			
			end	
		
		end

	end

	def Utility.readSingleFasta(inputFile)
		
		fasta = Hash.new
		puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Reading fasta file..."
		sequence = ""
		header = ""

		File.open(inputFile.to_s, "r") do |f|
	  		
	  		f.each_line do |line|

	  			if line.match('>')

	  				header = line.chomp

	  			else

	  				sequence.concat(line.chomp)

	  			end	
	  		end
		end

		return header, sequence

	end

	def Utility.splitFasta(input)
		
		temp = Array.new

		File.open(input.to_s, "r") do |file|

			header = Array.new	

			file.each_line do |line|				

				if line.include? ">"
					
					temp << line.split(" ")[0].tr(">", "")

					header = line.split(" ")

					Avail.fileCreate(header, line)

				else	
					
					Avail.fileCreate(header, line)
				
				end	

			end	

		end	

		return temp

	end		

	def Utility.getSequenceLength(file)

		data = Hash.new
		header=""
		len=nil

		 File.open(file.to_s, "r") do |file|
		 	
		 	file.each_line do |line|
		 		if line.include?(">")
		 			header = line
		 			next
		 		else
		 			len = line.length
		 		end	
		 		data[header] = len
		 	end	
		 end
		 return data	
	end	

	def Utility.intro 
		puts "							" 
		puts "###############################################################################################"
		puts "							"  
		puts "Welcome to transcript extension workflow"
		puts "This interactive ruby program will aid you in the extension of incompletely/partially"
		puts "assembled transcripts (generated by de-novo transcriptome assemblers). This can be "
		puts "achieved by first creating Database A (for sequence searching) and Database B (for " 
		puts "mate fetching) Once you retrieve a targeted subset of paired-end reads, use Geneious"
		puts "software to further map the targeted subset of reads against partial transcript. After"
		puts "mapping, take the consensus sequence and run this script again to fetch consensus sequence"
		puts "specific reads. Repeat these steps until you get a desired extension of your partial transcript"
		puts " "
		puts "For any queries, please contact me on this email or raise an issue on GitHub."
		puts ""	  
		puts "Kumar Saurabh Singh"
		puts "k.saurabh-singh@exeter.ac.uk"
		puts "							"
		puts "#################################################################################################"
		puts "							"
	end

end	
