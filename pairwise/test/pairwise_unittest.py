import unittest
import doctest

class DeviceTest( unittest.TestCase ):
	# This is a simple test that just tries to load the module
	def runTest( self ):
		try:
			import pairwise
			try:
				f = open("test/pairwise_unittest.py")
				pairwise.run( f )
				f.close()
			except IOError, e:
				self.fail( str( e ) )
		except ImportError, e:
			self.fail( str( e ) )

