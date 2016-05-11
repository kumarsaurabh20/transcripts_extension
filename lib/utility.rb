#!/usr/bin/ruby

require 'prerequisite'

module Utility

	BadRunError=Class.new(Exception)

	# Check and set path to all the required tools and software
	#Use prerequisite shell script
	# find out if some programmes are not installed.
	def Utility.checkAndSetTools
		check = Prerequisite.new
		check.checkTools
		exit
	end	

	# check the write permission of $workDir before building of the work directory
	def Utility.checkPermissions(file)
		#path = File.absolute_Path(file)
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

	def Utility.unzip(file)
		cmd = "gzip -d #{file}"
		%x[ #{cmd} ]
		if $?.exitstatus == 0
			return true
		else
			raise BadRunError, "Error: Oops. Can not run the shell command!"
		end	
	end	

	# check whether fasta file exist and how many sequences it contains
	def Utility.countSeq(file)
		cmd = "grep -c '^>' #{file}"
		count = `#{cmd}`
		return count
	end

	def convertQ2A(file)
		name = File.basename(file, ".*")
		if name.empty?
			puts "Error::FastQ File not found!"
		else
			cmd = "seqtk seq -A #{file} > #{name}.fasta"
			`#{cmd}`
			if $?.exitstatus == 0
				puts "FastA generated!"
			else
				raise BadRunError, "Error: Oops. Can not run the shell command!"
			end	 	  			
		end	
	end	

	def expunge(file)

	end	



	# check whether the necessary perl scripts exist and cand be found


	def Utility.intro 
		puts "							" 
		puts "##############################################################################"
		puts "							"  
		puts "Welcome to iTransMap"
		puts "An interactive ruby program for the extension of incompletely/partially"
		puts "assembled transcripts (generated by de-novo transcriptome assemblers)"
		puts "using whole transcriptome (RNA-Seq) data!"
		puts " 																		"	  
		puts "Kumar Saurabh Singh"
		puts "k.saurabh-singh@exeter.ac.uk"
		puts "							"
		puts "###############################################################################"
		puts "							"
	end

	def Utility.copyright
		puts "Copyright 2016 Singh KS"
		puts "Licensed under the Apache License, version 2.0 (the \"License\"); You may not use this file except in compliance with the License. You may obtain a copy of the license at http://www.apache.org/licenses/LICENSE-2.0 Unless required by applicable law or agreed to in writing, software distributed under the license is distributed on an \"AS IS\" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the license"
	end  


end	