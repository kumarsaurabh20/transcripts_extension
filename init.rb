#!/usr/bin/ruby

##iTransMap###

##Launch this file from the command line to start the application and get started###

APP_ROOT = File.dirname(__FILE__)
$:.unshift(File.join(APP_ROOT, 'lib'))

require "itransmap"

map = Itransmap.new
map.options(ARGV)