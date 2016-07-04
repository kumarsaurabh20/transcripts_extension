
require 'utility'
require 'fileutils'
require 'avail'
require 'sequence'
require 'querydb'

class CreateAndQueryDb

	BadRunError=Class.new(Exception)
	ArgumentError=Class.new(StandardError)
	NoMethodError=Class.new(NameError)
	


	attr_accessor :prefix, :step	

	def initialize(prefix)
    	 @prefix = prefix
    	 @step=Time.new
 	end

	def formatReadFiles(partials, searchAlgo)

		outfasta = ""
		
		begin

			currentPath = Utility.navigate("Sample_data")

			files = Dir.glob("*.fastq.gz")
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  #{files.inspect}" 
			
			if files.is_a?(Array) and files.size == 2

				puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  #{files.inspect} is an array!"

				outfastq = Utility.mergeZippedFiles(files[0].to_s, files[1].to_s, @prefix)
				puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  #{outfastq.to_s}"
				puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Converting from FASTQ to FASTA..."
				outfasta = Utility.convertQ2A(outfastq, @prefix)
				puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Files are reformatted to FASTA!"
				puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  #{outfasta.to_s}"
				Avail.createDir(nil, "DB")
				Avail.moveFile(currentPath, "DB", outfasta)
			else

				puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  ERROR::Data files could not be found!"
			
			end	

		rescue Exception => e

			e.backtrace
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  ERROR::Something bad happened!!"
			exit
		
		end

		Utility.navigate("DB")

		if Utility.checkFileExist(outfasta)

			name = File.basename(outfasta, ".*")
			createDbA(outfasta)
					
		else
			
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  ERROR::Cannot move merged fasta file to DB folder!"
			exit
		
		end

		createDbB

		query = Querydb.new(@prefix)

		if searchAlgo.eql?("blast")

			query.useBlast("DB", partials)

		elsif searchAlgo.eql?("nhmmer")	
			
			query.useNhmmer("DB", partials)

		else

			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Please enter the search algorithm!"
			raise ArgumentError

		end	

	end	

	def createDbA(file)

		name = File.basename(file, ".*")
		
		if Utility.checkFileExist(file)

			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Creating BLAST DB of raw sequence reads...."
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Check blast log file for detail BLAST DB description!"

			cmd = "makeblastdb -in #{@prefix}.fasta -input_type fasta -dbtype nucl -title #{@prefix} -out #{@prefix} -logfile #{@prefix}.log"
			Avail.executeCmd(cmd)
			#Avail.executeCmd("cat #{prefix}.log")	
		
		else
			
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  ERROR::Database directory does not exist!"
			exit
		
		end		

		puts " "
		puts "-------------------------------"
		puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Database A successfully created"
		puts "-------------------------------"
		puts " "

	end	

	def createDbB

		currentPath = Utility.navigate("Sample_data")

		r1 = "#{@prefix}_R1_trimmed.fastq"
		r2 = "#{@prefix}_R2_trimmed.fastq"

		if currentPath.include? "/Sample_data"
			
			rawReads = Dir.glob("*.gz")
			#sickle pe -f R1.fastq.gz -r R2.fastq.gz -t sanger -o R1_sickle_20.fastq -p R2_sickle_20.fastq -s lost-mate.fastq -q 20 -l 20
			#rawreads.each do |file|
			#	Utility.unzip(file)
			#end	
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  #{rawReads.inspect}"
			cmd = "sickle pe -f #{rawReads[0].to_s} -r #{rawReads[1].to_s} -t sanger -o #{r1} -p #{r2} -s singles.fastq"
			
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Creating database B from raw sequence reads...."
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Processing reads by trimming low quality bases (Threashold for quality is 20 and for length is 20bp)..."
			Avail.executeCmd(cmd)
			#FileUtils.rm "singles.fastq"

			puts " "
			puts "-------------------------------"
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Database B successfully created"
			puts "-------------------------------"
			puts " "
		
		else
			
			puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  ERROR::Possibly wandering in wrong directory!"
			exit
		
		end		

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
