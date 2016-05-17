
require 'utility'

class Createdbs

	def formatReadFiles(file, prefix)
		begin
			if file.is_a?(Hash)
				file.each do |value|
					name = File.basename(value, ".*")
					convertQ2A(Utility.unzip(value, "DB", "../"))
					puts "Decompressiong #{value}..."
					temp = File.join(name, ".fastq")
					convertQ2A(temp)
					faFiles = Dir.glob("*.fasta")
					if faFiles.size == 2 && faFiles.is_a?(Array)
						createDbFasta(prefix, array)
					else
						puts "Cannot Create merge fasta file."
					end	
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
		
		temp = makeDir("blastDB", "./../")
		Dir.chdir(temp)



	end	

	def useNhmmer


	end

	def queryReads


	end

	def fetchPairs


	end	




end	

#convert from zipped to fastq
		#fastq to fasta

		#irb(main):002:0> accepted_formats = [".txt", ".pdf"]
		#irb(main):003:0> File.extname("example.pdf") # get the extension

		#Dir.chdir("/var/spool/mail")
		#Dir.entries("testdir")
		#Dir.foreach("testdir") {|x| puts "Got #{x}" }
		#Dir.getwd           #=> "/tmp"
		#Dir.pwd             #=> "/tmp"	

		#Dir["config.?"]                     #=> ["config.h"]
		#Dir.glob("config.?")                #=> ["config.h"]
		#Dir.glob("*.[a-z][a-z]")            #=> ["main.rb"]
		#Dir.glob("*.[^r]*")                 #=> ["config.h"]
		#Dir.glob("*.{rb,h}")                #=> ["main.rb", "config.h"]
		#Dir.glob("*")                       #=> ["config.h", "main.rb"]
		#Dir.glob("*", File::FNM_DOTMATCH)   #=> [".", "..", "config.h", "main.rb"]

		#rbfiles = File.join("**", "*.rb")
		#Dir.glob(rbfiles)                   #=> ["main.rb",
                                    #    "lib/song.rb",
                                    #    "lib/song/karaoke.rb"]

        #libdirs = File.join("**", "lib")
		#Dir.glob(libdirs)                   #=> ["lib"]

		#librbfiles = File.join("**", "lib", "**", "*.rb")
		#Dir.glob(librbfiles)                #=> ["lib/song.rb",
                                    #    "lib/song/karaoke.rb"]

		#librbfiles = File.join("**", "lib", "*.rb")
		#Dir.glob(librbfiles)                #=> ["lib/song.rb"]                            

		#Dir.mkdir(File.join(Dir.home, ".foo"), 0700) #=> 0
