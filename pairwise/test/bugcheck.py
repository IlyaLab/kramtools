
"""
This script scans a TCGA feature matrix for features that might
have been affected by a bug discovered in pairwise.

It expects the matrix on stdin.
It assumes a header row is present (and skips it).
It emits lines of the form...

	|<featurename>|=C, range=[min,max]

...for "bad" features only--that is, features that might have been
affected by this bug. C is the cardinality of the feature's levels.
"""

import sys

assert sys.version_info.major >= 3

def _bad_line( values ):
	if all([ v == "NA" or v.isdigit() for v in values ]):
		values.discard("NA")
		cardinality = len(values)
		values = [ int(v) for v in values ]
		mn = min(values)
		mx = max(values)
		if mn < 0 or cardinality <= mx:
			return ( cardinality, mn, mx )
	return None


def _proc( fp ):
	line = fp.readline()
	while len(line) > 0:
		if line.startswith("C:"):
			fields = line.rstrip().split('\t')
			result = _bad_line( set(fields[1:]) )
			if result:
				print( "|{}|={}, range=[{},{}]".format( fields[0], *result ) )
		line = fp.readline()

_proc( sys.stdin )

