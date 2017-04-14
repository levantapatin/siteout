#!/usr/bin/env python
# Venky Iyer
# Time-stamp: <2008-09-23 15:45:13 viyer>

import re
import sys
import textwrap

class seq(object):

    '''
    A single FASTA sequence object.

    >>> import fasta

    >>> f = fasta.seq(fid='CG2328', sequence='ATCG', description='eve')
    
    >>> print f
    >CG2328 eve
    ATCG
    
    >>> fstr = str(f)
    
    >>> print fstr
    >CG2328 eve
    ATCG
    
    >>> f2 = fasta.seq.parse(fstr)
    
    >>> f2.fid
    'CG2328'
    
    >>> f2.description
    'eve'
    
    >>> f2.sequence
    'ATCG'

    '''

    def __init__(self, fid, sequence, description=None):

        # We use fid because id is a reserved word in Python
        self.fid = fid
        
        self.sequence = sequence
        self.description = description

    def __str__(self, width=80):

        ''' Pretty print the sequence '''

        header = '>%s' % self.fid

        if self.description is not None:
            header += ' %s' % self.description

        sequence = '\n'.join(textwrap.wrap(self.sequence, width))
        
        return '%s\n%s' % (header, sequence)

    def __len__(self):
        return len(self.sequence)

    @classmethod
    def parse(self, fasta_string):

        ''' Parse a string representation of a FASTA sequence into a fasta.seq
        object. 
        '''

        lines = fasta_string.strip().split('\n')

        (header, sequence) = (lines[0], lines[1:])

        # combine lines of sequence
        sequence = ''.join((x.strip() for x in sequence))
        
        # Parse bits out of header
        header = header.strip()
        
        header_regexp = '>(?P<fid>\S+)(?:\s+(?P<description>.+))?'
        header_match = re.match( header_regexp, header )

        if not header_match:
            raise ValueError( 'FASTA header %s did not match regexp %s' % (header, header_regexp) )

        fid = header_match.group('fid')
        description = header_match.group('description')

        return self(fid, sequence, description)

    @classmethod
    def tokenize(self, stream):

        '''
        Tokenize a stream (files, stdin etc) into fasta.seq objects.

        >>> fh = open('test.fasta', 'r')

       >>> for fs in fasta.seq.tokenize(fh):
       ...     print fs.fid
       ...     print len(fs)
       ...
       FBtr0071764
       5180
       FBtr0100521
       4644
       FBtr0071763
       4831
       FBtr0083388
       3946
 
        
       '''

        header = None
        sequence = None
        
        for line in stream:

            line = line.strip()
            
            if not line or line.startswith('#'):
                continue
            elif line.startswith('>'):
                
                if sequence:
                    yield self.parse('\n'.join((header, sequence)))

                header = line
                sequence = ''
                
            elif line is not None:
                sequence += line

        if sequence:
            yield self.parse('\n'.join((header, sequence)))

    def __getitem__(self, k):

        ''' Support slicing.

        >>> f = fasta.seq(fid='CG2328', sequence='ATCG'*20, description='eve')

        >>> len(f)
        80

        >>> print f
        >CG2328 eve
        ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
        

        >>> f2 = f[0:50]
        
        >>> len(f2)
        50
        
        >>> print f2
        >CG2328 eve
        ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT

        '''
        
        sliced = self.sequence[k]
        return self.__class__(self.fid, sliced, self.description)
    
    @classmethod
    def pairwisePercentID(self, *seqList):
    	'''
    	Determines the pairwise percent identity of several sequences.
    	
    	Returns a vector of the pairwise percent identities.
    	
    	Ignores '-' positions.
    	
    	'''
    	
    	id = []
    	bl = seq(fid = 'blank', sequence = '-')
    	for i in range(len(seqList)):
    		for j in range(i+1, len(seqList)):
    			idCount = 0.0
    			seq1 = seqList[i]
    			seq2 = seqList[j]
    			minLength = min(len(seq1), len(seq2))
    			for k in range(minLength):
    				if (seq1[k] == seq2[k] or seq1[k] == bl or seq2[k] == bl):
    					idCount = idCount + 1
    			id.append(idCount/minLength)
    	return id
    	
    @classmethod
    def globalPercentID(self, *seqList):
    	'''
    	Determines the global percent identity of several sequences.
    	
    	Returns the percent of sites which are completely identical.
    	
    	Optionally, the last argument can be the length over which the sequences
    	are compared.
    	
    	Ignores '-' positions.
    	
    	>>> f = fasta.seq(fid='CG2328', sequence='ATCG')
    	>>> f1 = fasta.seq(fid='CG2329', sequence='ATCGG')
    	>>> id = fasta.seq.globalPercentID(f, f1)
    	>>> print id
    	1.0
    	
    	>>> id = fasta.seq.globalPercentID(f, f1, 5)
    	>>> print id
    	0.8
    	
    	>>> seqList = [f, f1]
    	>>> id = fasta.seq.globalPercentID(*seqList)
    	
    	'''
    	
    	id = 0.0
    	bl = seq(fid = 'blank', sequence = '-')
    	seq1 = seqList[0]
    	
    	if type(seqList[len(seqList)-1]).__name__ == 'int':
    		seqList = list(seqList)
    		length = seqList.pop()	
    	else:
    		length = len(seq1)
    	
    	for i in range(length):
    		posConserved = 1
    		for j in range(1, len(seqList)):
    			seq2 = seqList[j]
    			if (min(len(seq1), len(seq2)) < (i + 1)):
    				posConserved = 0
    				break
    			if (seq1[i] == bl or seq2[i] == bl):
    				continue
    			elif (seq1[i] != seq2[i]):
    				posConserved = 0
    				break
    		if posConserved:
    			id = id + 1
    	return id/length
    	
    def __eq__(self, other):
    	if (self.sequence == other.sequence):
    		return 1
    	else:
    		return 0
    		
    def __ne__(self, other):
    	if (self.sequence != other.sequence):
    		return 1
    	else:
    		return 0
   
class file(object):

    ''' This will create and interface transparently to a SQLite database
    containing pickled versions of fasta.seq objects, in order to efficiently
    query and work with a very large FASTA file, such as all binding sites in a
    genome.'''
    
    pass
 

