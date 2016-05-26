
require 'utility'
require 'fileutils'
require 'avail'

class CreateAndQueryDb

	BadRunError=Class.new(Exception)

	def formatReadFiles(prefix)

		outfasta = ""
		
		begin

			currentPath = Utility.navigate("Sample_data")

			files = Dir.glob("*.fastq.gz")
			puts files.inspect
			
			if files.is_a?(Array) and files.size == 2

				puts "#{files.inspect} is an array!"

				outfastq = Utility.mergeZippedFiles(files[0].to_s, files[1].to_s, prefix)
				puts outfastq.to_s
				puts "Converting from FASTQ to FASTA..."
				outfasta = Utility.convertQ2A(outfastq, prefix)
				puts "Files are reformatted to FASTA!"
				puts outfasta.to_s
				Avail.createDir("DB")
				Avail.moveFile(currentPath, "DB", outfasta)
			else

				puts "ERROR::Data files could not be found!"
			
			end	

			# if file.is_a?(Hash)
			# 	file.each do |value|
			# 		name = File.basename(value, ".*")
					
			# 		if Utility.unzip(value, "data")
						
			# 			puts "Decompressiong #{value}..."
			# 			temp = File.join(name, ".fastq")
			# 			Utility.convertQ2A(temp)

			# 		else
			# 			puts "ERROR::Oops..something terrible happened while unzipping the file!"
			# 			exit	
			# 		end							
			# 	end	
			# else
			# 	Puts "WARNING::Input file is not paired!"
			# end	
		rescue Exception => e

			e.backtrace
			puts "ERROR::Something bad happened!!"
			exit
		
		end

		Utility.navigate("DB")

		if Utility.checkFileExist(outfasta)

			name = File.basename(outfasta, ".*")
			createDbA(outfasta, prefix)
					
		else
			
			puts "ERROR::Cannot move merged fasta file to DB folder!"
			exit
		
		end

		createDbB(prefix)

	end	

	def createDbA(file, prefix)

		name = File.basename(file, ".*")
		
		if Utility.checkFileExist(file)

			puts "Creating BLAST DB of raw sequence reads...."
			puts "Check blast log file for detail BLAST DB description!"

			cmd = "makeblastdb -in #{prefix}.fasta -input_type fasta -dbtype nucl -title #{prefix} -out #{prefix} -logfile #{prefix}.log"
			Avail.executeCmd(cmd)
			#Avail.executeCmd("cat #{prefix}.log")	
		
		else
			
			puts "ERROR::Database directory does not exist!"
			exit
		
		end		

		puts " "
		puts "-------------------------------"
		puts "Database A successfully created"
		puts "-------------------------------"

	end	

	def createDbB(prefix)

		currentPath = Utility.navigate("Sample_data")

		r1 = "#{prefix}_R1_trimmed.fastq"
		r2 = "#{prefix}_R2_trimmed.fastq"

		if currentPath.include? "/Sample_data"
			
			rawReads = Dir.glob("*.gz")
			#sickle pe -f R1.fastq.gz -r R2.fastq.gz -t sanger -o R1_sickle_20.fastq -p R2_sickle_20.fastq -s lost-mate.fastq -q 20 -l 20
			#rawreads.each do |file|
			#	Utility.unzip(file)
			#end	
			puts rawReads.inspect
			cmd = "sickle pe -f #{rawReads[0].to_s} -r #{rawReads[1].to_s} -t sanger -o #{r1} -p #{r2} -s singles.fastq"
			
			puts "Creating database B from raw sequence reads...."
			puts "Processing reads by trimming low quality bases (Threashold for quality is 20 and for length is 20bp)..."
			Avail.executeCmd(cmd)
			#FileUtils.rm "singles.fastq"

			puts " "
			puts "-------------------------------"
			puts "Database B successfully created"
			puts "-------------------------------"
			puts " "
		
		else
			
			puts "ERROR::Possibly wandering in wrong directory!"
			exit
		
		end		

	end

	def useBlast(dir, partial ,threads, prefix)
		
		count = 0
		if threads == 0
			count = Utility.countSeq(partial)
		else
			count = threads
		end	

		#blastn -query partials.fasta -db R1 -out blastn_out_R1.tab -outfmt 6 -max_target_seqs 50000 -num_threads 25
		if Utility.directory_exists?("DB")
			
			cmd = "blastn -query #{partial}.fasta -db #{prefix} -out #{prefix}.out -outfmt 6 -max_target_seqs 50000 -num_threads #{count}"
			
			puts "Executing blast search...."
			
			Utility.executeCmd(cmd)	
		else
			puts "ERROR::Database directory does not exist!"
			exit
		end	
		
	end	

	def useNhmmer(partial, prefix, threads)
		
		count = 0
		if threads == 0
			count = Utility.countSeq(partial)
		else
			count = threads
		end	

		#nhmmer --tblout pacbio_nillu_hmmer.blastout --noali --incE 0.0001 -E 0.001 --max --dna --cpu 16 ER1_vF_CYP6ER1vR1.fasta 8830_1_2_3_4_5_merged.fasta
		if Utility.directory_exists?("DB")
			
			cmd = "nhmmer --tblout #{prefix}.out --noali --incE 0.001 -E 0.001 --max --dna --cpu #{count} #{partial}.fasta #{prefix}.fasta"
			
			puts "Executing nhmmer search...."
			
			Utility.executeCmd(cmd)	
		
		else
			
			puts "ERROR::Database directory does not exist!"
			exit
		
		end	


	end

	def fetchPairs(dbdir, datadir, searchAlgo)

		if Utility.directory_exists?(dbdir) && Utility.directory_exists?(datadir)

			Utility.navigate(dbdir)
			searchResultOutput = Dir.glob("*.out")
			Utility.moveFilesToTmp(searchResultOutput)
			
			Utility.navigate(datadir)
			rawDataFiles = Dir.glob("*.fastq")
			Utility.moveFilesToTmp(rawDataFiles)
			Utility.navigate("tmp")

						
				if searchAlgo.eql?("blast")		
					
					cmd_1 = "awk '{print $1}' #{searchResultOutput.to_s} > #{searchResultOutput.to_s}_filtered.out"
					Utility.executeCmd(cmd)
				
				else
					
					cmd_1 = "awk '{print $2}' #{searchResultOutput.to_s} > #{searchResultOutput.to_s}_filtered.out"
					Utility.executeCmd(cmd)
				
				end	

				rawDataFiles.each_with_index do |value, index| 
					
					cmd_2 = "seqtk subseq  #{value} #{searchResultOutput.to_s}_filtered.out > targeted_R#{index}.fasta"	
					
					Utility.executeCmd(cmd)
				
				end	
		else
			
			puts "ERROR::Either DB or Data directory is not available!!"		
		
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
