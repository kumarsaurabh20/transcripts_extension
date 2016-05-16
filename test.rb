#!/usr/bin/ruby

require 'bio'

class Test

	def countSeq(file)
		partials = Bio::FlatFile.new(Bio::FastaFormat,file)
		partials.each_entry do |record|
  			puts [record.definition, record.nalen.to_s ].join(" ")
		end
	end
end

pp = Test.new
pp.countSeq(ARGV.to_s)
