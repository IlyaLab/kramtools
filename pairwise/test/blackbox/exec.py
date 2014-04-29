
import sys
import random
import subprocess
import os.path
import math
import os

from covpair import CovariatePair
from randfeat import Feature,Num,Cat

_PAIRWISE_INPUT = '/dev/shm/z.tab'
_R_INPUT = os.path.splitext(_PAIRWISE_INPUT)[0]+'.ref'
_MAX_K = 8

_MAX_SAMPLES     = 1000
# User must explicitly set environment variables to turn on 
# degeneracy generation.
_PROB_UNI_DEGEN = float( os.environ.get( "UNIDEGEN", 0.0 ) ) #-1.00
_PROB_COV_DEGEN = float( os.environ.get( "COVDEGEN", 0.0 ) ) #-1.00
_REQUIRED_EXECUTABLES = ( 'pairwise', )

# Verify user's directory is set up to run.
# _REQUIRED_EXECUTABLES may require user setup (symlinking)...
for exe in _REQUIRED_EXECUTABLES:
	if not os.path.exists( exe ):
		print( "Either copy the relevant program or symlink to",exe,"in the local dir", file=sys.stderr )
		sys.exit(-1)
# ...but we can create missing directories ourselves.
if not os.path.isdir('failures'):
	os.mkdir( 'failures' )

"""
Two types of degeneracies:
	1) those inherent in one (or both) univariate features
	2) those that only emerge in the covariate pair because
	   NAs in one feature exclude some values in the other
I need to simulate both.

Covariate degeneracies happen when the missing values in one
feature effectively decimate present elements in the other.
Degeneracy types:
	Content-independent (low counts)
		1) 0 surviving pair
		2) 1 surviving pair
		3) 2 surviving pairs
	Content-dependent
		1) categorical: single surviving category
		2) numeric:     all surviving values identical, i.e. constant
"""

def _generateTestCase():
	"""
	Generate a covariate pair of random length, optionally with a random
	type and amount of degeneracy.
	"""
	N = random.randint(3,_MAX_SAMPLES)
	########################################################################
	# Choose covariate degeneracy first as it will constrain the types
	# of covariates required to manifest it.
	degen = None
	count = None
	const = None # ...unless set otherwise below.
	if random.random() < _PROB_COV_DEGEN:
		degen = random.choice(["count","const"])
		const = random.choice(["num","cat"])
		if degen=="count":
			count = random.randint(0,2)
			types = random.randint(0,3) # any features will do
		elif const=="num":
			types = random.randint(0,2) # pair must have a numerical feature
		else:
			types = random.randint(1,3) # pair must have a categorical feature
	else: # no degeneracy in this case
		types = random.randint(0,3) # any features will do
	########################################################################
	# Now generate covariates.
	# Both can contain univariate degeneracy.
	k = [ random.randint( 1 if random.random() < _PROB_UNI_DEGEN else 2, _MAX_K ) for n in range(2) ]
	r = [ random.randint( 1, N ) if random.random() < _PROB_UNI_DEGEN else 0, 0 ]
	i = random.randint(0,1) # to insure both parameters not same
	# Generate the covariate pair. 
	# It might contain univariate degeneracy, but it does not yet contain
	# any (forced) covariate degeneracy...though some could have occurred
	# entirely by chance.
	if types == 0:   # NN
		p = CovariatePair( Num(N,rep=r[i]),  Num(N,rep=r[1-i]) )
	elif types == 1: # CN
		p = CovariatePair( Cat(N,k[0]),      Num(N,rep=r[0])   )
	elif types == 2: # NC
		p = CovariatePair( Num(N,rep=r[0]),  Cat(N,k[0])       )
	else:        # CC
		p = CovariatePair( Cat(N,k[i])     , Cat(N,k[1-i])     )
	########################################################################
	# Now muck up the covariates if degeneracy was selected above.
	if degen:
		if degen=="count":
			p.fewSurvivors( count )
		elif const=="num":
			p.constScalar()
		else:
			p.constCategory()
	# Return a "specification" that contains the covariate data and a
	# qualitative description of it with respect to presence/absence of 
	# degeneracy.
	return {'n':N, 
		'degen':degen,
		'param':count if degen=="count" else const, 
		'types':types, 
		'covar':p }

def _closeEnough( f1, f2 ):
	return (math.isnan(f1) == math.isnan(f2)) or (abs(f1-f2) < 1e6)

def _verify( pout, rout, case ):
	"""
	Returns a brief string describing a discrepancy or None if there
	is no discrepancy.
	"""
	assert case in frozenset([0,1,2,3])
	pfields = pout.rstrip().split('\t')
	# pairwise' output fields (when standard format active)
	#  5: p-value
	# 12: Spearman rho
	# 13: log
	pval = float(pfields[5])
	prho = float(pfields[3])
	plog =       pfields[10]
	rfields = rout.rstrip().split('\t')
	if case <= 2: # NN, CN, NC
		# R output should have p-value, rho
		if not _closeEnough( pval, float(rfields[0]) ):
			return "p-value"
		if case == 0 and not _closeEnough( prho, float(rfields[1]) ):
			return "rho"
	else: # CC
		if not _closeEnough( pval, float(rfields[0]) ):
			return "p-value"
		if plog != rfields[1]:
			return "cull-log"
	return None


def _saveFailedData( testnum, data ):
	"""
	"""
	assert isinstance( data , CovariatePair )
	# Just re-writing since trying to rename from /dev/shm generated
	# "cross-devic link" error
	with open("failures/{0:04x}.tab".format(testnum),"w") as fp:
		data.writeForPairwise( fp )
	with open("failures/{0:04x}.ref".format(testnum),"w") as fp:
		data.writeForR( fp )


def _report( testnum, status, spec ):
	print( "{0}\t{1}\t{n}\t{degen}\t{param}\t{types}".format( status, testnum, **spec ) ) 

############################################################################

TEST_COUNT = int(sys.argv[1]) if len(sys.argv) > 1 else 1

R_SCRIPTS = (  # these MUST correspond to the types field in the dict 
	'num.R', # returned by _generateTestCase.
	'mix.R',
	'mix.R',
	'cat.R' )

good = 0
for tnum in range(TEST_COUNT):
	# Generate a test case
	spec = _generateTestCase()
	# Generate files
	with open(_PAIRWISE_INPUT,"w") as fp:
		spec['covar'].writeForPairwise( fp )
	with open(_R_INPUT,"w") as fp:
		spec['covar'].writeForR( fp )
	# Run the executables
	try:
		p_out = subprocess.check_output(
			['./pairwise','-h','-f','std',_PAIRWISE_INPUT], 
			universal_newlines=True )
		r_out = subprocess.check_output(
			['R','--slave','-f',R_SCRIPTS[spec['types']],'--args',_R_INPUT ], 
			universal_newlines=True )
		try:
			discrep = _verify( p_out, r_out, spec['types'] )
		except IndexError:
			# when pairwise does not produce output
			discrep = "silent"
		if discrep:
			_saveFailedData( tnum, spec['covar'] )
		# Always write test/reference output...
		sys.stdout.write( '#'+p_out )
		sys.stdout.write( '#'+r_out )
		# ...then a report.
		if not discrep:
			_report( tnum, "ok", spec )
			good += 1
		else:
			_report( tnum, discrep, spec )
	except subprocess.CalledProcessError as e:
		sys.stdout.write( '#\n#\n' ) # indicate ABSENCE of output for consistency
		_report( tnum, "procerror({0})".format(e.cmd), spec )
		_saveFailedData( tnum, spec['covar'] )

print( "{} ok, {} possible failures".format( good, TEST_COUNT-good ) )
