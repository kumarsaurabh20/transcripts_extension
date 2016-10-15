#!/usr/bin/ruby

##Transcript Extension###

##Launch this file from the command line to start the application and get started###

APP_ROOT = File.dirname(__FILE__)
$:.unshift(File.join(APP_ROOT, 'lib'))

require 'extension'
#require "Testup"


map = Extension.new
#map.options(ARGV)
map
#test = Testdir.new
#test.runCommand("cat", "file1", "file2", "file3", "file4")

#test = Testup.new
#test.check
