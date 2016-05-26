#!/usr/bin/ruby

require 'fileutils'

module Avail

	BadRunError=Class.new(Exception)

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

	def Avail.createDir(dir)

		path = File.expand_path("../../", __FILE__)
		FileUtils.cd(path)
		
		begin

			Dir.exist?(dir)
			FileUtils.mkdir dir
			puts "#{dir} directory created!"		

		rescue Exception => e

			e.backtrace
			puts "Waoo, #{dir} already exist!"
		
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

	#http://code.tutsplus.com/tutorials/ruby-for-newbies-working-with-directories-and-files--net-18810
	def Avail.makeDir(dir)
		
		setpath = self.setDir
		
		FileUtils.cd(setpath)
		
		self.checkDirAtRootAndMake(dir)
		
		return true	
	end

	#Multiple arguments
	#http://stackoverflow.com/questions/831077/how-do-i-pass-multiple-arguments-to-a-ruby-method-as-an-array
	def Avail.expunge(file)

	end

end	


