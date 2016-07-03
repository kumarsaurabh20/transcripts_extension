#require 'testdir'

class Testup

	def check(input, size)

		filterList = {}
		header = []
		sequence = ""

		File.open(input.to_s, "r")  do |file|
			
			file.each_line do |line|

				if header.empty? && line.start_with?(">")
					header = line.split(" ")[0]
					#puts "this is line 1"
					#puts line.to_s
					
				elsif !header.empty? && line.start_with?(">")

						if !sequence.empty? 
								
								#puts sequence
								
								filterList["#{header}"] = sequence  if !file.eof? && sequence.size <= size
								#puts filterList.to_s
								#puts "hash is just filled with #{line.chomp}"

								header = line.split(" ")[0]
								#puts "new header #{header} is created"
								
								sequence = ""				

						else
								next
						end	
				
				else 					
					sequence.concat(line.chomp)
					#puts "sequence is just filled with #{line.chomp}"
					filterList["#{header}"] = sequence if file.eof? && sequence.size <= size			
				end		

			end	

		end		

		#puts filterList.to_s
		writeFile("FilteredByLength.fasta", filterList)	

	end	


	def writeFile(outputFile, fileObject, headerPrefix = nil)
		

		if fileObject.is_a? Array 

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

			raise "Something terrible has happened. Try again!!"

		end	

		
	end


end	

file = Testup.new
file.check(ARGV[0].to_s, ARGV[1].to_i)