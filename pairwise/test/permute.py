
"""
Generate a permutation of length <length>, possibly leaving the first
<fix> items in place. That is, if the second optional argument is present
then a permutation is generated that only permutes the last <length>-<fix>
items.
"""

import sys
import random

assert sys.version_info.major >= 3

def permutation( n, fix=0 ):
	v = [ i for i in range(n) ]
	for i in range(fix,n):
		j = random.randint( fix, i ) # j is in the CLOSED interval [0,i]
		v[i] = v[j]
		v[j] = i
	return [ n for n in v ]

if __name__=="__main__":
	if len(sys.argv) < 2:
		print( "{0} <length> <fix>".format( sys.argv[0] ) )
		sys.exit(0)
	P = permutation( int(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 0 )
	print( P )
	#print( '\n'.join( [ str(v) for v in P ] ) )

