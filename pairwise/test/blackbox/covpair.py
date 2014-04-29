
import sys
import random
from randfeat import Feature,Num,Cat

def _randomperm(n):
	l = list(range(n))
	random.shuffle(l)
	return l

class CovariatePair:
	def __init__(self,f1,f2):
		self.f1 = f1
		self.f2 = f2

	def _names( self ):
		if type(self.f1) == type(self.f2):
			n1 = "{}1".format( self.f1.typePrefix() )
			n2 = "{}2".format( self.f2.typePrefix() )
			return (n1,n2)
		else:
			return ( self.f1.typePrefix(), self.f2.typePrefix() )

	def anyCategorical( self ):
		"""
		Returns a randomly selected categorical feature or
		None if neither is categorical.
		"""
		return random.choice( [
			self.f1 if self.f1.isCategorical() else None,
			self.f2 if self.f2.isCategorical() else None ] )

	def anyNumeric( self ):
		"""
		Returns a randomly selected numeric feature or
		None if neither is numeric.
		"""
		return random.choice( [
			self.f1 if self.f1.isNumeric() else None,
			self.f2 if self.f2.isNumeric() else None ] )

	def other( self, f ):
		return self.f1 if self.f2 is f else self.f2

	# Next few methods are corrupters, introducers of degeneracy!

	def fewSurvivors( self, n ):
		"""
		Produce covariate pair with n survivors (pairs in which
		neither member is NA).
		"""
		N = len(self.f1)
		pair = (self.f1,self.f2)
		# Leave n intact at the beginning, then...
		for i in range(n,N):
			random.choice(pair).decimate(i)
		# ...shuffle them in among the rest.
		perm = _randomperm(N)
		self.f1.permute( perm )
		self.f2.permute( perm )

	def constCategory( self ):
		"""
		Produce a covariate pair in which (one of) the categorical
		feature(s) is reduced to one category by NAs in the other
		feature.
		"""
		degen = self.anyCategorical()
		if degen:
			other = self.other(degen)
			for i in degen.complementRandomCategory():
				other.decimate(i)
		elif __debug__:
			print( "warning: no categorical feature", file=sys.stderr )

	def constScalar( self ):
		"""
		Produce a covariate pair in which (one of) the categorical
		feature(s) is reduced to one category by NAs in the other
		feature.
		"""
		degen = self.anyNumeric()
		if degen and degen.hasRepeats():
			other = self.other(degen)
			for i in degen.complementRepeats():
				other.decimate(i)
		elif __debug__:
			print( "warning: no numeric feature or no repeats in it", file=sys.stderr )

	def writeForPairwise( self, fp ):
		name = self._names()
		print( name[0]+":", self.f1, sep="\t", file=fp )
		print( name[1]+":", self.f2, sep="\t", file=fp )

	def writeForR( self, fp ):
		print( *self._names(), sep="\t", file=fp )
		for i in range(len(self.f1)):
			print( self.f1[i], self.f2[i], sep="\t", file=fp )

############################################################################
# Unit testing

if __name__=="__main__":
	import sys
	import random
	N = int(sys.argv[1])
	DEGEN = sys.argv[2] if len(sys.argv) > 2 else None
	x = random.randrange(4)
	k = random.randint(2,4)
	r = random.randint(2,8)
	# Generate a normal covariate pair. "Normal" means it may
	# well have univariate degeneracies like too few categories,
	# but no (explicitly-introduced) covariate degeneracies.
	if x == 0:   # NN
		rr = random.randint(2,4)
		p = CovariatePair( Num(N,rep=r),  Num(N,rep=rr) )
	elif x == 1: # CN
		p = CovariatePair( Cat(N,k),      Num(N,rep=r)  )
	elif x == 2: # NC
		p = CovariatePair( Num(N,rep=r),  Cat(N,k)      )
	else:        # CC
		kk = random.randint(2,4)
		p = CovariatePair( Cat(N,k)     , Cat(N,kk)     )
	# If a degeneracy argument was provided, muck up the covariate pair.
	if DEGEN:
		if DEGEN.isdigit(): # Decimate the 
			p.fewSurvivors( int(DEGEN) )
		elif DEGEN.lower().startswith("const"):
			if DEGEN.lower().endswith("cat"):
				if x > 0: # if any feature is categorical
					p.constCategory()
				else:
					print( "warning: no categorical features", file=sys.stderr )
			elif DEGEN.lower().endswith("num"):
				if x < 3: # if any feature is numeric
					p.constScalar()
				else:
					print( "warning: no numeric features", file=sys.stderr )
			else:
				raise RuntimeError("unknown option:"+DEGEN)
		else:
			 raise RuntimeError("unknown option:"+DEGEN)
	p.writeForPairwise( sys.stdout )

