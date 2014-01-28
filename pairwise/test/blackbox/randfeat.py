
import random

"""
Covariate degeneracies happen when the missing values in one
feature effectively decimate present elements in the other.
Degeneracy types:
	Content-independent (low counts)
		1) 0 surviving pair
		2) 1 surviving pair
		3) 2 surviving pairs
	Content-dependent
		Categorical
			1) single surviving category
		Numeric
			1) all surviving values identical

Start by creating "good" covariate pair and decimate accordingly.
Run prep.py -H -O on it
Run ./pairwise on result, reading res
Run R on the data
Compare the results.
"""

def _randStr( length ):
	return ''.join([chr(65+random.randint(0,25)) for i in range(length)])

class Feature:
	def __init__( self ):
		pass
	def __len__( self ):
		return len(self.data)
	def __str__( self ):
		return '\t'.join(self.data)

	def __getitem__( self, i ):
		return self.data[i]

	def decimate( self, i ):
		self.data[i] = 'NA'

	def permute( self, p ):
		self.data = [ self.data[i] for i in p ]

	def isNumeric( self ):
		return not self.isCategorical()

	def isCategorical( self ):
		raise NotImplementedError("subclass should define this")

class Num(Feature):
	"""
	Generates and holds a random numeric feature with n samples in [0,Max).
	Optionally, generates 
	"""
	def __init__(self,n,rep=0,Max=100):
		if not n >= 2:
			raise RuntimeError("counts cannot be < 2")
		assert rep >= 0
		self.repeats = None
		# Create a random vector of n-rep
		v = [ random.uniform(0.0,Max) for i in range(n-rep) ]
		# ...and add some repeated elements, if requested.
		if rep > 0:
			r = random.uniform(0.0,Max)
			v.extend([ r for i in range(rep) ])
			random.shuffle( v )
			self.repeats = set( filter(lambda i:r == v[i],list(range(n))) )
		# Convert to strings; we don't need float values anymore.
		self.data = [ "{0:.1f}".format(f) for f in v ]

	def isCategorical( self ):
		return False

	def typePrefix(self):
		return 'N'

	def hasRepeats( self ):
		return self.repeats != None

#	def repeatedValue( self ):
#		if not self.hasRepeats():
#			raise RuntimeError("repeatValue requested on feature without repeats")
#		return self.data[ self.repeats[0] ]

	def complementRepeats( self ):
		"""
		Note that this does not imply that none of the elements at the
		returned indices are duplicates! ...only that these are not the
		ones explicitly generated in the constructor.
		"""
		universe = list(range(len(self.data)))
		return list(filter(lambda i: not (i in self.repeats), universe ))


class Cat(Feature):

	def __init__(self, n,k=3,lablen=3):
		self.labels = [ _randStr(lablen) for i in range(k) ]
		self.data = [ random.choice( self.labels ) for i in range(n) ]

	def isCategorical( self ):
		return True

	def typePrefix(self):
		return 'C'

	def complementRandomCategory( self ):
		"""
		Returns the indices of 
		"""
		x = random.choice(self.labels)
		universe = list(range(len(self.data)))
		return list(filter(lambda i: self.data[i] != x, universe))


############################################################################
# Unit testing

if __name__=="__main__":
	import sys
	N = int(sys.argv[2])
	if sys.argv[1].upper().startswith("N"):
		R = int(sys.argv[3]) if len(sys.argv) > 3 else 0
		print( str( Num(N,rep=R) ) )
	else:
		K = int(sys.argv[3]) if len(sys.argv) > 3 else 3
		print( str( Cat(N,K) ) )

