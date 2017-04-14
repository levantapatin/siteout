#!/usr/bin/python

'''
Takes a directory of patser results files and makes them csv's for use with inSite

Outputs one file per DNA sequence

Usage:
$ patesr2csv.py fastaFile patserResultDirectory outputDirectory 

'''

import os
import re
import sys
import fasta

seqFile = sys.argv[1]
resultsDir = sys.argv[2]
outputDir = sys.argv[3]

seqs = {}

#make header of each csv file based on fasta sequences
fh = open(seqFile, 'r')
for fs in fasta.seq.tokenize(fh):
	seqs[fs.fid] = ['#NAME: %s' % (fs.fid), '#START: 0', '#LENGTH: %i' % (len(fs)), '##class,start,length,motif_type, strength']
fh.close()

#for each result file extract list of binding sites and convert to InSite format
resultFiles = os.listdir(resultsDir)
for elem in resultFiles:
	fh = open(os.path.join(resultsDir, elem), 'r')
	
	nameMatch = re.search('(\w+)\_', elem) #extract binding site name from file name
	pwmName = nameMatch.group(1)
	
	pwmWidth = 0 #default
	
	for line in fh:
		pwmWidthMatch = re.search('^width of the matrix.+?(\d+)', line.strip()) #check if line has width info
		if (pwmWidthMatch):
			pwmWidth = int(pwmWidthMatch.group(1)) #set actual pwm width
		
		lineRe = '(\w+?)\.txt\s+position=\s+(\d+).+score=\s+([\d\.]+)' #regex of patser output format
		lineMatch = re.search(lineRe, line.strip()) #check if line has binding site info
		if (lineMatch):  #if it has binding site
			if(seqs.has_key(lineMatch.group(1))):  #check what sequence binding site comes from
				#now put binding site info in InSite format and save in 'seqs' holder
				seqs[lineMatch.group(1)].append('binding_site,%i,%i,%s,%.2f' % (int(lineMatch.group(2)), pwmWidth, pwmName, float(lineMatch.group(3))))
	fh.close()

#go through 'seqs' and write a .csv file for each sequence
for elem in seqs.keys():
	fh = open(os.path.join(outputDir, '%s.csv' % elem), 'w')
	fh.write('\n'.join('%s' % (x) for x in seqs[elem]))
	fh.close()



		
			