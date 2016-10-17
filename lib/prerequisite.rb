#!/usr/bin/ruby

require 'date'
require 'pp'

class Prerequisite
#check for 
#seqtk, sickle, ncbi-blast, nhmmer
#Date.strptime("2012-09-21 19:45:48","%Y-%m-%d %H:%M:%S")
VERSION="1.0.0"
@@softwarePath = Hash.new

	def checkTools
		software = %w[seqtk sickle blastn makeblastdb nhmmer] 
		software.each do |tool|
				toolExists(tool) 
		end	
		writeConfig
		puts " -------------- "
		@@softwarePath.each do |key, value|	
			#print "#{key}\t\t\t#{value}\n"

			printf("%-10s %-10s\n", key, value)
		end	
		puts " -------------- "
		puts "System check completed. Please check the config file!"
	end

	def toolExists(params)

		cmd = "which #{params}"
		path = %x[ #{cmd} ]

		if $?.exitstatus == 0
	  		@@softwarePath[params] = path.strip	
		else
			if path.empty?
				@@softwarePath[params] = "WARNING::Not Found. Please check the installation!"
			else
				displayErr(params)
				exit
			end				
		end					
	end

	def writeConfig
		current = Time.new
		#puts "Current Time : " + current.inspect
		#puts @@softwarePath.inspect
		File.open("itransmap.config", 'a') do |file|
			file.write "-------------------------------\n"
			file.write "#{current.inspect}\n"
			file.write "-------------------------------\n"
			file.write "iTransMap\tVersion-#{VERSION}\n"
			file.write "-------------------------------\n"
			file.write "Tool\tPath\n"
			file.write "-------------------------------\n"
			@@softwarePath.each do |key, value|	
				file.printf("%-10s %-10s\n", key, value)
			end	
		end	
	end	

	def displayErr(params) 
		puts "ERROR::Something terrible happend. Contact on k.saurabh-singh@exeter.ac.uk"
	end

end	

#pc = Prerequisite.new
#pc.checkTools
