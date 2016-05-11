#!/usr/bin/ruby

require 'utility'
require 'optparse'
require 'optparse/time'
require 'ostruct'

class Itransmap

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
		options.vital = true
		options.threads = 1

		parser = OptionParser.new do |opts|
			opts.banner = "Usage: itransmap [options]"
			opts.separator ""
      		opts.separator "Required options:"
      		opts.separator "-----------------"

			opts.on('-s', '--transcripts <string>', 'incompletely assembled or partial transcripts. Accepts a multifasta file.') do |transcripts|
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
					check = Utility.fileType(reads[0].to_s)
					present = Utility.checkPermissions(reads[0].to_s)
					present2 = Utility.checkPermissions(reads[1].to_s)

					if check.eql?("zipped") && present == true && present2 == true
						options.raw_reads[:pe_1] = reads[0];
						options.raw_reads[:pe_2] = reads[1];
					else
						puts "Error::Either the file is not gz compressed or the files are not available!"
						exit
					end		
					#Utility.checkPermissions					
			end
			opts.on('-a', '--algorithm <blast|hmmer>', 'Search algorithm (blast for faster searches and hmmer for a sensitive search).') do |algorithm|
					options.algorithm = algorithm;	
			end
			opts.on('-o', '--outfile <String>', 'Final output fasta file with extended transcripts.') do |outfile|
					options.outfile = outfile;	
			end
			opts.on('-t', '--num_threads <int>', 'Number of threads to run the program (Default: number of partial sequences equal to number of threads).') do |thread|
					options.threads = thread;	
			end

			opts.separator ""
      		opts.separator "Common options:"
      		opts.separator "---------------"

			opts.on('-p', '--precheck', 'Use this option to check the availability of required tools in the system') do |precheck|
				Utility.checkAndSetTools
			end
			opts.on_tail('-h', '--help', 'Displays help') do
		        puts opts
		        exit
		    end
		    opts.on_tail('-v', '--version', 'Show version') do
		        puts "iTransMap: version #{VERSION}"
		        exit
		    end
		    opts.on_tail('-i', '--info', 'Display tool information') do
		        Utility.intro
		        exit
		    end
		    opts.on_tail('-c', '--copyright', 'Display copyright information') do
		        Utility.copyright
		        exit
		    end
		    opts.on(' ', ' ', '') do
				puts opts
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

