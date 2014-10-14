
"""
Create categorical feature pairs that (depending on the argument) violate
cardinality limits for the purpose of checking pairwise' response to this
type of input.
"""

import sys
import permute
import random

def _randStr( length ):
	return ''.join([chr(65+random.randint(0,25)) for i in range(length)])

if len(sys.argv) < 2:
	print( sys.argv[0], "<# discrete values>", "[ <string label length>]", file=sys.stderr )
	sys.exit(-1)

R = 10
N = int(sys.argv[1])
L = int(sys.argv[2]) if len(sys.argv) > 2 else 0

# Create a base sequence consisting of either N random strings
# or the sequence [0..N-1] repeated R times.
C = [ _randStr(L) if L > 0 else i for i in range(N) ]*R

# Permute it once to create the 1st vector...
P = permute.permutation(N*R)
print( "C:x", *[ C[i] for i in P ], sep="\t" )

# ...and again to create the 2nd vector.
P = permute.permutation(N*R)
print( "C:y", *[ C[i] for i in P ], sep="\t" )

