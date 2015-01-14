
"""
It is intended to provide sample input files primarily for regression
testing pairwise.
This is NOT intended to validate the statistical tests, which is
why there are no statistical parameters to tweak the random output.
"""

import sys
import math
import random

_NAN_STRINGS  = ('na','nA','Na','NA')

Ncont = int(sys.argv[1])
Ndisc = int(sys.argv[2])
Nsamp = int(sys.argv[3])

def _randStr( length ):
	return ''.join([chr(65+random.randint(0,25)) for i in range(length)])

# Write a header...

print( "TCGA", *[ "{}{}".format(_randStr(3),i) for i in range(Nsamp) ], sep="\t" )

# ...and Ncont+Ndisc data lines.

while Ncont+Ndisc > 0:
	n = random.random()
	if n < 0.2:    # 20% of features have possibly large NA counts
		Nemp = random.randint( 0, Nsamp )
	elif n < 0.95: # 75% have small NA counts
		Nemp = random.randint( 0, Nsamp // 2 )
	else:
		Nemp = 0   # 05% have no missing data
	data = [ random.choice( _NAN_STRINGS ) for i in range(Nemp) ]
	if random.random() > ( Ncont / (Ncont+Ndisc) ):
		####################################################################
		# Categorical
		# Choose a random number of categories...
		Ncat = random.randint( 2, math.ceil( math.sqrt(Nsamp)) )
		# ...with random names, insuring no duplicates (so the actual
		# category count my by < Ncats).
		cats = list( frozenset( [ _randStr(3) for i in range(Ncat) ] ) )
		# Build the data 
		data += [ random.choice(cats)
			for i in range(Nsamp-Nemp) ]
		random.shuffle( data )
		print( "{}:{}".format( "B" if len(cats) == 2 else "C", _randStr(3) ),
			*data, sep="\t" )
		Ndisc -= 1
	else:
		####################################################################
		# Numeric/continuous
		data += [ "{0:.1f}".format(random.uniform(0,100))
			for i in range(Nsamp-Nemp) ]
		random.shuffle( data )
		print( "N:{}".format( _randStr(3) ), *data, sep="\t" )
		Ncont -= 1
	assert len(data) == Nsamp

