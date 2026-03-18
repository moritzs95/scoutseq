#!/usr/bin/python
"""
Filter interleaved paired-end FASTQ records by sequence motifs in read 2 to remove e.g. pseudogene reads of off-target reads.
"""

import sys

#get fname from parameter
patternFile=sys.argv[1]

#load patterns
patterns = set( x.strip() for x in open(patternFile) )

#read interleaved paired-end fastq from stdid
handle=sys.stdin

#parse fastq
while True:
  name=handle.readline()
  if not name:
    break
  seq     =handle.readline()
  spacer  =handle.readline()
  quals   =handle.readline()
  nameR2=handle.readline()
  seqR2   =handle.readline()
  spacerR2=handle.readline()
  qualsR2 =handle.readline()
  #print('%s%s%s%s%s%s%s%s' % ( idline, seq, spacer, quals, idlineR2, seqR2, spacerR2, qualsR2 ))
  #check
  if all(x not in seqR2 for x in patterns):
    #print fastq
    sys.stdout.write( '%s%s%s%s%s%s%s%s' % ( name, seq, spacer, quals, nameR2, seqR2, spacerR2, qualsR2 ) )
