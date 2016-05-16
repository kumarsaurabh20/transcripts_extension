
require 'utility'

class Createdbs

	BadRunError=Class.new(Exception)

	def formatReadFiles(file, prefix)
		begin
			if file.is_a?(Hash)
				file.each do |value|
					name = File.basename(value, ".*")
					
					if Utility.unzip(value, "data")
						
						puts "Decompressiong #{value}..."
						temp = File.join(name, ".fastq")
						Utility.convertQ2A(temp)

					else
						puts "ERROR::Oops..something terrible happened while unzipping the file!"
						exit	
					end							
				end	
			else
				Puts "WARNING::Input file is not paired!"
			end	
		rescue Exception => e
			e.backtrace
			puts "Oopss! Faulty file conversion. Check the raw read files!"
		end

		puts "Files are reformatted to FASTA!"
		puts "Merging the Fasta files for their reuse in BLAST/HMMER searches!"

		faFiles = Dir.glob("*.fasta")
		
		if faFiles.is_a?(Array) && faFiles.size == 2
			Utility.createDbFasta(prefix, faFiles)
		else
			puts "ERROR::Cannot Create merge fasta file."
			exit
		end

	end	

	def createDbA(dir, prefix)
		
		if Utility.directory_exists?("DB")
			cmd = "makeblastdb -in #{prefix}.fasta -input_type fasta -title #{prefix} -parse_seqids -out #{prefix} -logfile #{prefix}.log"
			puts "Creating BLAST DB of raw sequence reads...."
			`#{cmd} 1> /dev/null 2> /dev/null`
			if $?.exitstatus == 0
				puts "Blast DB created!!"
			else
				raise BadRunError, "Error: Oops..Bad command execution!"
			end		
		else
			puts "ERROR::Database directory does not exist!"
			exit
		end		

	end	

	def useBlast(dir, partials ,threads, prefix)
		
		count = 0
		if threads == 0
			count = Utility.countSeq(partials)
		else
			count = threads
		end	

		#blastn -query partials.fasta -db R1 -out blastn_out_R1.tab -outfmt 6 -max_target_seqs 50000 -num_threads 25
		if Utility.directory_exists?("DB")
			cmd = "blastn -query #{partials}.fasta -db #{prefix} -out #{prefix}.out -outfmt 6 -max_target_seqs 50000 -num_threads #{count}"
			puts "Executing blast search...."
			`#{cmd} 1> /dev/null 2> /dev/null`
			if $?.exitstatus == 0
				puts "Blast search is now completed!!"
			else
				raise BadRunError, "Error: Oops..Bad command execution!"
			end		
		else
			puts "ERROR::Database directory does not exist!"
			exit
		end	
		
	end	

	def useNhmmer(partials, prefix, threads)
		
		count = 0
		if threads == 0
			count = Utility.countSeq(partials)
		else
			count = threads
		end	

		#nhmmer --tblout pacbio_nillu_hmmer.blastout --noali --incE 0.0001 -E 0.001 --max --dna --cpu 16 ER1_vF_CYP6ER1vR1.fasta 8830_1_2_3_4_5_merged.fasta
		if Utility.directory_exists?("DB")
			cmd = "nhmmer --tblout #{prefix}.out --noali --incE 0.001 -E 0.001 --max --dna --cpu #{count} #{partials}.fasta #{prefix}.fasta"
			puts "Executing nhmmer search...."
			`#{cmd} 1> /dev/null 2> /dev/null`
			if $?.exitstatus == 0
				puts "nhmmer search is now completed!!"
			else
				raise BadRunError, "Error: Oops..Bad command execution!"
			end		
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
					`#{cmd_1} 1> /dev/null 2> /dev/null`
				else
					cmd_1 = "awk '{print $2}' #{searchResultOutput.to_s} > #{searchResultOutput.to_s}_filtered.out"
					`#{cmd_1} 1> /dev/null 2> /dev/null`
				end	

				rawDataFiles.each_with_index do |value, index| 
					cmd_2 = "seqtk subseq  #{value} #{searchResultOutput.to_s}_filtered.out > targeted_R#{index}.fasta"	
					`#{cmd_1} 1> /dev/null 2> /dev/null`
				end	

		else
			puts "ERROR::Either DB or Data directory is not available!!"		
		end		


	end	

	def createDbB(dir, prefix)

		if Utility.directory_exists?("Data")
			
			rawReads = Dir.glob("*.gz")
			#sickle pe -f R1.fastq.gz -r R2.fastq.gz -t sanger -o R1_sickle_20.fastq -p R2_sickle_20.fastq -s lost-mate.fastq -q 20 -l 20
			cmd = "sickle pe -f #{rawReads[0]} -r #{rawReads[1]} -t sanger -o R1_sickle_20.fastq -p R2_sickle_20.fastq -q 20 -l 20"
			puts "Creating database B from raw sequence reads...."
			`#{cmd} 1> /dev/null 2> /dev/null`
			if $?.exitstatus == 0
				puts "Database B created!!"
			else
				raise BadRunError, "Error: Oops..Bad command execution!"
			end		
		else
			puts "ERROR::Database directory does not exist!"
			exit
		end	

		#############################################################################
		#	Check if you need to move file to tmp or you can work from Data directory
		#############################################################################	

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
