#!/usr/bin/env python
# Zeba Wunderlich

import re
import sys
import math
import types
import os.path

import fasta

class fm(object):

	'''
	A single frequency matrix object.  
	
	Matrix rows correspond to binding site positions.  
	Matrix columns are in the order A C G T.

	>>> import pwm

	>>> f = pwm.fm(matrix=[[0.1788, 0.3692, 0.1000, 0.3519], [0.2173, 0.0615, 0.0231, 0.6981], [0.9481, 0.0038, 0.0423, 0.0058]])	
	
	>>> print f
	0.1788	0.3692	0.1000	0.3519	
	0.2173	0.0615	0.0231	0.6981	
	0.9481	0.0038	0.0423	0.0058	
	
	>>> fstr = str(f)
	
	>>> print fstr
	0.1788	0.3692	0.1000	0.3519	
	0.2173	0.0615	0.0231	0.6981	
	0.9481	0.0038	0.0423	0.0058	
		
	>>> f2 = pwm.fm.parse(fstr)
	
	>>> f2.matrix
	[[0.1788, 0.3692, 0.1000, 0.3519], [0.2173, 0.0615, 0.0231, 0.6981], [0.9481, 0.0038, 0.0423, 0.0058]]
	
	>>> fstr2 = '0.25	0.5 0.1 0.3 0.9\n0.25	0.5 0.2 0.3 0\n0.25 0	0.6 0.3 0.05\n0.25	0	0.1 0.1 0.05'
	
	>>> f3 = pwm.fm.parse(fstr2, 'v')
	
	>>> len(f3)
	5
	
	>>> print f3
	0.2500	0.2500	0.2500	0.2500
	0.5000	0.5000	0.0000	0.0000
	0.1000	0.2000	0.6000	0.1000
	0.3000	0.3000	0.3000	0.1000
	0.9000	0.0000	0.0500	0.0500

	'''

	def __init__(self, matrix):

		self.matrix = matrix
		
	def __str__(self):
		''' Pretty print the matrix.'''
		matrixString = '\n'.join('\t'.join(map(str, row)) for row in self.matrix)
		return matrixString

	def __len__(self):	 
		'''Return number of positions in PWM.'''
		return len(self.matrix)
	
	def __getitem__(self, k):
		''' Support slicing.

		>>> f = pwm.fm(matrix=[[0.1788, 0.3692, 0.1000, 0.3519], [0.2173, 0.0615, 0.0231, 0.6981], [0.9481, 0.0038, 0.0423, 0.0058]])	

		>>> len(f)
		3

		>>> f2 = f[0:2]
		
		>>> len(f2)
		2
		
		>>> print f2
		0.1788	0.3692	0.1000	0.3519	
		0.2173	0.0615	0.0231	0.6981	

		'''		  
		sliced = self.matrix[k]
		return self.__class__(sliced)
	
	@classmethod
	def parse(self, fm_string, orientation='h'):
		''' Parse a string of a frequency matrix into a pwm.fm object. 
		
		Keyword arguments:
		orientation -- 'v' if matrix is vertical
		'''
		lines = fm_string.strip().split('\n')
		matrix = []
	   
		for x in range(len(lines)):
			row_match = re.split('\s+', lines[x].strip())			 
			matrix.append(map(float, row_match))
			
		if (orientation == 'v'):
			matrix = zip(*matrix)
			
		for x in range(len(matrix)):
			if len(matrix[x]) != 4:
				raise ValueError("Position %i of matrix isn't 4 columns" % (x+1))
			if (sum(matrix[x]) < 0.99 or sum(matrix[x]) > 1.01):
				raise ValueError("Frequencies at position %i do not sum to 1" % (x+1))
			
		return self(matrix)
				
	def weightMatrix(self, gcContent=0.5):
		'''Calculate weight matrix from a frequency matrix.
		
		Keyword arguments:
		gcContent -- fractional GC content, must be between 0 and 1
		'''
		if (gcContent <= 0 or gcContent >= 1):
			raise ValueError("gcContent must be between 0 and 1")
					
		bg = [0.5 - gcContent/2, gcContent/2, gcContent/2, 0.5 - gcContent/2]
		wm = []
		
		for i in range(len(self)):
			wmRow = [0, 0, 0, 0]
			for j in range(4):
				if self.matrix[i][j] != 0:
					wmRow[j] =	math.log(self.matrix[i][j]/bg[j], 2)
			wm.append(wmRow)
		return wm

	def klEntropy(self, gcContent=0.50):
		'''Calculate KL entropy for a frequency matrix.
				
		Keyword arguments:
		gcContent -- fractional GC content, must be between 0 and 1
		'''
		if (gcContent <= 0 or gcContent >= 1):
			raise ValueError("gcContent must be between 0 and 1")
		
		bg = [0.5 - gcContent/2, gcContent/2, gcContent/2, 0.5 - gcContent/2]
		kl = 0
		
		for i in range(len(self)):
			for j in range(4):
				if self.matrix[i][j] != 0:
					kl += self.matrix[i][j] * math.log(self.matrix[i][j]/bg[j], 2) 
		return kl
		
	def seqWeight(self, sequence, gcContent=0.50, concentration=1):
		'''Calculate the statistical weight of a given sequence.
		
		Sequence can be given as a string of a fasta object
		
		Keyword arguments:
		gcContent -- fractional GC content, must be between 0 and 1
		concentration -- number of TFs
		'''
		
		if (gcContent <= 0 or gcContent >= 1):
			raise ValueError("gcContent must be between 0 and 1")
		if (concentration < 0):
			raise ValueError("concentration must be positive")
		if (type(sequence) is not types.StringType): #fix: find better way
			sequence = sequence.sequence
		
		bg = [0.5 - gcContent/2, gcContent/2, gcContent/2, 0.5 - gcContent/2]
		weight = concentration
		
		sequence = sequence.upper()
		sequence = re.sub('A', '0', sequence)
		sequence = re.sub('C', '1', sequence)
		sequence = re.sub('G', '2', sequence)
		sequence = re.sub('T', '3', sequence)
		
		for i in range(len(self)):
			weight = weight * self.matrix[i][int(sequence[i])]/bg[int(sequence[i])]
		return weight

	def patser(self, gcContent = 0.5, filename = 'fileName', cutoff = 0, minPValue = 1):
		'''Submit a weight matrix to patser to search the sequences in fileName.
		Returns a dict where keynames are sequence names and values are a 2D matrix
		where the first entry is the match location and the second is log(pValue)
		
		Keyword arguments:
		gcContent -- fractional GC content, must be between 0 and 1
		filename -- name of file containing names of sequence files
		cutoff -- lower score threshold used by patser
		minPValue -- minimum p-value of hit to be included in output matrix
		
		Note that temporary files alphabet, matrix and temp.out are written out.
		'''
		
		if (gcContent <= 0 or gcContent >= 1):
			raise ValueError("gcContent must be between 0 and 1")
		if (os.path.isfile(filename) is False):
			raise ValueError("%s does not exist" % fileName)
		if (cutoff < 0):
			raise ValueError("cutoff must be greater than 0")
		if (minPValue < 0 or minPValue > 1):
			raise ValueError("minPValue must be between 0 and 1")
			
		tempFile = 'temp.out'
		hitLocations = {}
		logPValue = math.log(minPValue)
					
		
		#Write out alphabet and matrix files for patser
		fh = open('alphabet', 'w')
		fh.write('A:T\nC:G\n\n')
		fh.close()
		
		wm = self.weightMatrix(gcContent)
		wmString = '\n'.join(["%.4f\t%.4f\t%.4f\t%.4f" % (row[0], row[1], row[2], row[3]) for row in wm])
		
		fh = open('matrix', 'w')
		fh.write('A\tC\tG\tT\n')
		fh.write("%s\n" % wmString)
		fh.close()
		
		'''
		submit to patser
		-w = weight matrix, -v = vertical matrix, -p = print matrix, -f = filename, -c = complementary, -l = lower score threshold -d2 = ignore unrecognized characters
		'''
		os.system("patser-v3b -w -v -p -f %s -c -l %.4f -d2 > %s " % (filename, cutoff, tempFile))
		
		#filter through patser results		
		fh = open(tempFile, 'r')
		for line in fh:
			lineRe = '([\w\.\-]+)\s+position=\s+(\d+).+ln\(p-value\)=\s+(.+)'
			lineMatch = re.search(lineRe, line.strip())
			if (lineMatch):
				if (float(lineMatch.group(3)) < logPValue):
					if (lineMatch.group(1) in hitLocations):
						hitLocations[lineMatch.group(1)].append([int(lineMatch.group(2)), float(lineMatch.group(3))])
					else:
						hitLocations[lineMatch.group(1)] = [[int(lineMatch.group(2)), float(lineMatch.group(3))]]
		fh.close()
		
		return hitLocations
	
class cm(object):

	'''
	A single count matrix object.  
	
	Matrix rows correspond to binding site positions.  
	Matrix columns are in the order A C G T.

	>>> import pwm

	>>> c = pwm.cm(matrix=[[1, 5, 10, 12], [12, 5, 1, 10], [6, 0, 11, 11]]) 
	
	>>> print c
	1	5	10	12
	12	5	1	10
	6	0	11	11
	
	>>> len(c)
	3	
		
	>>> c2 = pwm.cm.parse(str(c))
	
	>>> c2.matrix
	[[1, 5, 10, 12], [12, 5, 1, 10], [6, 0, 11, 11]]
	
	
	>>> f = c.cm2fm()
	
	>>> print f
	0.0357142857143 0.178571428571	0.357142857143	0.428571428571
	0.428571428571	0.178571428571	0.0357142857143 0.357142857143
	0.214285714286	0.0 0.392857142857	0.392857142857
	'''

	def __init__(self, matrix):

		self.matrix = matrix
		
	def __str__(self):
		''' Pretty print the matrix. '''
		matrixString = '\n'.join('\t'.join(map(str, row)) for row in self.matrix)
		return matrixString

	def __len__(self):	 
		'''Return number of positions in PWM.'''
		return len(self.matrix)
		
	def __getitem__(self, k):

		''' Support slicing.

		>>> c = pwm.cm(matrix=[[1, 5, 10, 12], [12, 5, 1, 10], [6, 0, 11, 11]]) 

		>>> len(f)
		3

		>>> c2 = c[0:2]
		
		>>> len(c2)
		2
		
		>>> print c2
		1	5	10	12
		12	5	1	10	

		'''
		
		sliced = self.matrix[k]
		return self.__class__(sliced)
	
	@classmethod
	def parse(self, cm_string, orientation='h'):
		''' Parse a string of a count matrix into a pwm.cm object. 
		
		Keyword arguments:
		orientation -- 'v' is matrix is vertical		
		'''

		lines = cm_string.strip().split('\n')
		matrix = []
	   
		for x in range(len(lines)):
			row_match = re.split('\s+', lines[x].strip())			 
			matrix.append(map(int, row_match))
			
		if (orientation == 'v'):
			matrix = zip(*matrix)
			
		for x in range(len(matrix)):
			if len(matrix[x]) != 4:
				raise ValueError("Position %i of matrix isn't 4 columns" % (x+1))
			
		return self(matrix)
	
	@classmethod	
	def sites2cm(self, siteMatrix):
		'''Parse a matrix that contains one binding site/element into a pwm.cm object.  All sites must be the same length.
		'''
		matrix = []
		letterMapping = {'A':0, 'C':1, 'G':2, 'T':3}
		
		for i in range(len(siteMatrix[1])):
			matrix.append([0, 0, 0, 0])
		
		for i in range(len(siteMatrix)):
			siteMatrix[i] = siteMatrix[i].upper()
			
			for j in range(len(siteMatrix[i])):
				matrix[j][letterMapping[siteMatrix[i][j]]] = matrix[j][letterMapping[siteMatrix[i][j]]]+1
		
		return self(matrix)
	
	def cm2fm(self, psuedoCount=0.0):
		'''Convert a count matrix into a frequency matrix, using psuedoCount.
		
		Keyword arguments:
		psuedoCount -- added to each count in the count matrix
		'''
		cMatrix = self.matrix
		fMatrix = []
		
		for x in range(len(cMatrix)):
			rowTotal = float(sum(cMatrix[x])) + float(psuedoCount)*4.0
			fMatrix.append(map(lambda y: (float(y) + float(psuedoCount))/rowTotal, cMatrix[x]))
		
		return fm(matrix=fMatrix)
