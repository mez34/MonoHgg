#!/usr/bin/python

import sys, getopt, re, os

def main(argv):
   inputfile = ''
   outputfile = ''
   outputdir = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:d:",["ifile=","ofile=","odir"])
   except getopt.GetoptError:
      print 'extractJSONS.py -i <inputfile> -o <outputfile> -d <outputdir>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'extractJSONS.py -i <inputfile> -o <outputfile> -d <outputdir>'
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      elif opt in ("-d", "--odir"):
         outputdir = arg
   print 'Input file is ', inputfile
   hand = open(inputfile) # open catalogue from Flashgg

 
   outputname = outputdir+'/'+outputfile+'.json'
   #check that output directory exists and if not create it
   dir = os.path.dirname(outputname)
   if not os.path.exists(dir):
      os.makedirs(dir)   

   target = open(outputname,'w') #output file with the EOS file of dataset to be analyzed

          
   lines = hand.readlines()
   hand.close()

   #print out the line+next 2 (so name,nentries,weight)
   for i in range(0, len(lines)):
      line = lines[i]
      if outputfile in line and i+2 < len(lines):
         target.write(lines[i])
         target.write(lines[i+1])
         target.write(lines[i+2])
	 target.write(lines[i+3])

   target.close()        
   print 'Output file is ', outputname

if __name__ == "__main__":
   main(sys.argv[1:])
