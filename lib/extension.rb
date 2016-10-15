#!/usr/bin/ruby

require 'utility'
require 'optparse'
require 'optparse/time'
require 'ostruct'
require 'createAndQueryDb'
require 'querydb'

class Extension

	VERSION = 1.0

	BadRunError=Class.new(Exception)
	ArgumentError=Class.new(StandardError)
	NoMethodError=Class.new(NameError)

	o = OpenStruct.new
        o.transcripts = nil
        o.algorithm = nil
        o.prefix = nil

	ARGV << '-h' if ARGV.size==0
	OptionParser.new do |opts|
		opts.banner = "Usage: itransmap [options]"
		opts.separator ""
      		opts.separator "Required options:"
      		opts.separator "-----------------"

		opts.on('-s', '--transcripts <string>', 'incompletely assembled or partial transcripts. Accepts a fasta file with single partial transcript or a consensus sequence.') do |transcripts|
			#check = Utility.fileType(transcripts)
			o.transcripts = transcripts
		end

		opts.on('-a', '--algorithm <blast|nhmmer>', 'Search algorithm (blast for faster searches and hmmer for a sensitive search).') do |algorithm|
			o.algorithm = algorithm	
		end

		opts.on('-p', '--prefix <String>', 'Database Prefix.') do |prefix|
			o.prefix = prefix
		end

		opts.separator ""
      		opts.separator "Common options:"
      		opts.separator "---------------"

		opts.on_tail('-k', '--precheck', 'Use this option to check the availability of required tools in the system') do |precheck|
			Utility.checkAndSetTools
			exit
		end

		opts.on_tail('-d', '--database <database_prefix>', 'Use this option to create blast and hmmer database') do |prefix|
			create = CreateAndQueryDb.new(prefix)
			create.formatReadFiles("partials.fasta", "blast")
			exit
		end
			
		opts.on_tail('-v', '--version', 'Show version') do |version|
			puts "Extension workflow: version #{VERSION}"
		        exit
		end

		opts.on_tail('-i', '--info', 'Display tool information') do |info|
			Utility.intro
		        exit
		end

		opts.on_tail('-h', '--help', 'Displays all options') do
			puts opts
			exit
	    	end
	
	end.parse!
				
	abort "-s is mandatory" if o.transcripts.nil?
	abort "-a is mandatory" if o.algorithm.nil?
	abort "-p is mandatory"	if o.prefix.nil?
		
	#puts o.algorithm.to_s
	#puts o.transcripts.to_s
	#puts o.prefix.to_s

	begin
		query = Querydb.new(o.prefix.to_s)

                if o.algorithm.eql?("blast")
                        query.useBlast("DB", o.transcripts.to_s)
                elsif o.algorithm.eql?("nhmmer")
                        query.useNhmmer("DB", o.transcripts.to_s)
                else
                        puts "[#{@step.strftime("%d%m%Y-%H:%M:%S")}]  Please enter the search algorithm!"
                        raise ArgumentError
                end

	rescue => err
		$stderr.puts "Exception: #{err}\n\n"
		err.backtrace.each { |l| $stderr.puts l + "\n" }
		err	
	end	
end
