#!/usr/bin/ruby
##Transcript Extension###
##Launch this file from the command line to start the application and get started###

APP_ROOT = File.dirname(__FILE__)
$:.unshift(File.join(APP_ROOT, 'lib'))
require 'extension'
map = Extension.new
map