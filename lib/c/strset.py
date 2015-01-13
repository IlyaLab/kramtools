
"""
This script pounds on an executable generated from strset.c
that, in turn, exercises the hash set implementation.
1. A random set of ASCII strings of random lengths is generated, and
	a random sampling with replacement of this set is written to a
	tmpfile (so the tmpfile *will* have duplicates).
2. The test executable is run which inhales this tmpfile, adds all its
	strings to a hash set, then dumps the contents of the hashset back
	out to stdout.
3. This script reads back the dumped hashset and verifies that the set
	read back is 1) a set and 2) equal to the original set.

It immediately aborts on any failure printing the name of a
tmpfile in which the failure manifests.
"""

import sys
import random
import subprocess
import tempfile
import os

assert sys.version_info.major >= 3

def _randUpper( length ):
	return ''.join([chr(65+random.randint(0,25)) for i in range(length)])

def _randAscii( length ):
	return ''.join([chr(33+random.randint(0,93)) for i in range(length)])

def _test( executable, string_count, max_strlen=32 ):
	strings = [ _randAscii( random.randint(1,max_strlen) )
		for x in range(string_count) ]
	strings = frozenset(strings)
	sample  = []

	# Create an array with random counts of duplicates of each string
	# in the set.

	for s in strings:
		for i in range(random.randint(1,10)):
			sample.append(s)

	# Mix it up.

	random.shuffle( sample )

	# Write it to a tmp file

	tmpfile_name = None
	with tempfile.NamedTemporaryFile( mode="w+", delete=False ) as f:
		tmpfile_name = f.name
		for s in sample:
			print( s, file=f )

	# Let the test executable inhale the file from which it will derive
	# and write back out a non-redundant set.
	# "2" is the starting capacity of the set (forces testing of resizing)
	# "1" tells it to strdup entries (tests memory mgmt)

	roundtrip_result = set()
	count = 0
	with subprocess.Popen([executable,"2","1",tmpfile_name ], stdout=subprocess.PIPE ) as proc:
		while True:
			line = proc.stdout.readline()
			if len(line) < 1:
				break
			roundtrip_result.add( line.decode().rstrip() )

	# Verify that the set read back from the test executable is identical
	# to the set that generated the array (with duplicates) sent the
	# executable.

	if roundtrip_result == strings:
		os.unlink( tmpfile_name )
		return None
	else:
		return tmpfile_name

if len(sys.argv) < 2:
	print( "{} <executable> [<max strings per test>] [<# tests>]".format( sys.argv[0] ) )
	sys.exit(-1)

M = int(sys.argv[2]) if len(sys.argv) > 2 else 10
N = int(sys.argv[3]) if len(sys.argv) > 3 else 1
while N > 0:
	N -= 1
	n = random.randint( min(2,M), M )
	result = _test( sys.argv[1], n )
	if result:
		print( "Failure on {}".format( result ) )
		break
	else:
		print( "OK {} strings".format( n ) )

