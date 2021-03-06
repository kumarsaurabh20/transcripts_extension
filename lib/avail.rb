#!/usr/bin/ruby

require 'fileutils'

module Avail

	BadRunError=Class.new(Exception)
	ArgumentError=Class.new(StandardError)
	NoMethodError=Class.new(NameError)
    SequenceFormatError=Class.new(Exception)

  # Error raised when FASTA file is malformed
	class DataFormatError < IOError
    	def message
      		"Data format error -- check input file"
    	end
  	end

  	@step = Time.now


	def Avail.writeFile(outputFile, fileObject, headerPrefix = nil)		
		if fileObject.is_a? String || fileObject.integer?
				File.open(outputFile.to_s, 'a') do |file|
					fileObject.each do |line|
						file.write "#{fileObject}\n"
					end	
				end
		elsif fileObject.is_a? Array 
			File.open(outputFile.to_s, 'a') do |file|
				fileObject.each_with_index do |element, index|
					file.write "#{headerPrefix}#{index}\n#{element}\n"
				end	
			end
		elsif fileObject.is_a? Hash	
				File.open(outputFile.to_s, 'a') do |file|
					fileObject.each do |key, value|
						file.write "#{key}\n#{value}\n"
					end	
				end
		else
			raise "Something terrible has happened while writing file. Try again!!"
		end	
	end

	def Avail.check_blast_out(file)
		return File.zero?(file)
	end

	def Avail.fileCreate(header, line)
		File.open("#{header[0].tr(">", "")}.fasta", "a") do |f|
			f.puts "#{line}"
		end	
	end	

	def Avail.executeCmd(cmd)

		`#{cmd}` if cmd.is_a?(String)
			
		if $?.exitstatus == 0
			return true
		else
			raise BadRunError, "Error: Oops..Bad command execution!"
		end		

	end

	def Avail.moveFile(src, dest, file)
		sourcePath = "#{src}/#{file}"
		destPath = "#{File.expand_path("../../", __FILE__)}/#{dest}"
		FileUtils.mv sourcePath, destPath, :force => true, :verbose => true
		FileUtils.chmod 0755, "#{destPath}/#{file}"
	end	

	def Avail.createDir(basedir=nil, dir)
		begin
			temp = File.expand_path("../../", __FILE__)

			if basedir != nil	
				path = File.join(temp, basedir)
				FileUtils.cd(path)
				puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Changed to #{path} directory"
			else
				FileUtils.cd(temp)
				puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Changed to #{path} directory"
			end	
			
			FileUtils.mkdir dir if !Dir.exist?(dir)
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Created #{dir}"

		rescue Exception => e
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Warning::#{dir} already exist!"
		end		
	end	

	def Avail.setDir
		checkPath = File.expand_path("../../init.rb", __FILE__)
		return File.dirname(checkPath) 
	end

	def Avail.checkDirAtRootAndMake(dir)
		if File.directory?(dir)
			FileUtils.rm_rf(dir)
			FileUtils.mkdir dir
			return true
		else	
			FileUtils.mkdir dir
			return true
		end	
	end	

	def Avail.expunge(files=[])
		FileUtils.rm files	
	end	
end	
