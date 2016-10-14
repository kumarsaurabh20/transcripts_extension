#!/usr/bin/ruby	

require 'utility'
require 'optparse'
require 'optparse/time'
require 'ostruct'
require 'createAndQueryDb'

class Extension

	VERSION = 1.0
	
	#begin  
	#  	raise 'A test exception.'  
	#	rescue Exception => e  
	#	  	puts e.message  
	#	  	puts e.backtrace.inspect  
	#end  


	class FileError < StandardError  
	end 

	#FileError = Class.new(Exception)

	def options(args) 
		options = OpenStruct.new
		options.transcripts = nil
		options.raw_reads = {}
		options.algorithm = nil
		options.outfile = nil
		options.prefix = nil
		options.vital = true
		options.threads = 1

		parser = OptionParser.new do |opts|
			opts.banner = "Usage: itransmap [options]"
			opts.separator ""
      		opts.separator "Required options:"
      		opts.separator "-----------------"

			opts.on('-s', '--transcripts <string>', 'incompletely assembled or partial transcripts. Accepts a fasta file with single partial transcript or a consensus sequence.') do |transcripts|
					check = Utility.fileType(transcripts)
					if check.eql?("fasta")
						options.transcripts = transcripts;
					else
						puts "Error::Check the partial file type!"
						exit
					end		
			end

			opts.on('-r','--raw <pe_1,pe_2>','Paired end raw RNAseq data. The final R1 and R2 file are merged files (replicates and groups).') do |pe|
					
					reads = pe.split(',')

					########################
					#Check for file paths if on file names are provided then expand the names
					########################

					check = Utility.fileType(reads[0].to_s)
					present = Utility.checkPermissions(reads[0].to_s)
					present2 = Utility.checkPermissions(reads[1].to_s)

					if check.eql?("zipped") and present == true and present2 == true
						options.raw_reads[:pe_1] = reads[0]
						options.raw_reads[:pe_2] = reads[1]
					else
						puts "Error::The raw read files are either not gz compressed or the files are not available!"
						exit
					end		
					#Utility.checkPermissions					
			end

			opts.on('-a', '--algorithm <blast|hmmer>', 'Search algorithm (blast for faster searches and hmmer for a sensitive search).') do |algorithm|
					options.algorithm = algorithm	
			end

			#opts.on('-o', '--outfile <String>', 'Final output fasta file with extended transcripts.') do |outfile|
			#		options.outfile = outfile;	
			#end

			opts.on('-p', '--prefix <String>', 'Database Prefix.') do |prefix|
					options.outfile = prefix	
			end

			#opts.on('-t', '--num_threads <int>', 'Number of threads to run the program (Default: number of partial sequences equal to number of threads).') do |thread|
			#		options.threads = thread;	
			#end

			opts.separator ""
      		opts.separator "Common options:"
      		opts.separator "---------------"

			opts.on_tail('-h', '--help', 'Displays help') do |help|
		        puts opts
		        puts " "
		        exit
		    end

			opts.on_tail('-k', '--precheck', 'Use this option to check the availability of required tools in the system') do |precheck|
				Utility.checkAndSetTools
				exit
			end

			opts.on_tail('-d', '--database <database_prefix>', 'Use this option to create blast and hmmer database') do |prefix|
				create = CreateAndQueryDb.new(prefix)
				create.formatReadFiles("partials.fasta", "blast")
				exit
			end

			#opts.on_tail('-m', '--mapping <mapping_prefix>', 'Use this option to map the reads against partial transcripts') do |prefix|
			#	create = CreateAndQueryDb.new(prefix)
			#	create.formatReadFiles("partials.fasta", "blast")
			#	exit
			#end
			
		    opts.on_tail('-v', '--version', 'Show version') do |version|
		        puts "Extension workflow: version #{VERSION}"
		        exit
		    end

		    opts.on_tail('-i', '--info', 'Display tool information') do |info|
		        Utility.intro
		        exit
		    end

		    opts.on_tail('-c', '--copyright', 'Display copyright information') do |copyright|
		        Utility.copyright
		        exit
		    end

		    opts.on(' ', ' ', '') do |nothing|
				puts opts
				puts " "
		        exit
			end

		# options specifications go here
			#f !args
			#		opts.on_tail('-h', '--help', 'Displays help') do
			#        puts opts
			#        exit
		    #	end
		    #end	
		end	
		parser.parse!
		puts options.to_s
	end

	# if no working directory is set, use current directory	
end

