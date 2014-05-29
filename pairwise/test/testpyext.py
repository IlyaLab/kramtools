
"""
This script merely tests that the Python extension pairwise.so
is loadable and runnable.
This script should 
"""

import sys
import tempfile
import io

try:
	import pairwise
except ImportError:
	print """Add the directory containing pairwise.so to your PYTHONPATH environment 
variable. e.g. PYTHONPATH=<whatever> python %s""" % sys.argv[0]
	sys.exit(-1)

_RANDOM_MATRIX="""HEADER	A	B	C	D	E	F	G	H
ehYSX	5.4578334e+08	-2.8507611e+08	-4.3018952e+08	7.1989919e+08	9.1774179e+08	-3.4384220e+08	NA	-1.5198611e+08
EPftK	zxBqp	zxBqp		lIOWU	lIOWU	MqHit	TMRPm	lIOWU
pSbjR	+490206467.09833	+564394923.03903	+771901139.75079	-225271250.65280	+929190595.73979	-695924195.60089	-710747222.12050	+742906094.84894
myihR	+237082349.308	-286646212.095	-363622130.949	-463150322.513	-591150767.727	-338147341.062	-807692895.817	
SzgAR	ItCUG	ItCUG	NA	ItCUG	GekWw	NA	bjvHC	GekWw
VPHLQ	DXJZC	wirWa	UwQpz	DXJZC	wirWa	wirWa	hCKZA	
GDpjt	8.153e+08	-8.851e+08	8.493e+08	-5.547e+07	1.803e+08		9.815e+05	3.661e+08
"""

# Write the dummy data above to a tmpfile.
i = tempfile.TemporaryFile()
i.write( _RANDOM_MATRIX )
# ...then use it as input to pairwise.
i.seek( 0, io.SEEK_SET )
pairwise.run( i, sys.stdout )
i.close()

