from setuptools import setup, Extension

setup(
    # Name of this package
    name="pairwise",
    
    # Package version
    version=0.1,
    
    # This tells setup how to find our unit tests.
    test_suite = "test.pairwise_unittest",
    
    # Describes how to build the actual extension module from C source files.
    ext_modules = [
        Extension(
          'pairwise',
          ['src/py.c']
        )]
    )
