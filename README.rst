============================================================================
ABOUT
============================================================================

This repository actually contains multiple *independent* projects.
Each is contained in its own directory.

The ./lib directory contains generic code shared between one or more project 
subfolders (and organized by language). 
Each subfolder is self-contained and independently buildable with the
sole exception that any project may depend on the contents and relative
location of lib.
In other words, **the directory hierarchy must be preserved**.

Because I actively try to identify 
general-purpose routines within projects and factor them out the lib 
directory contains a very *ad hoc* mix of things.

