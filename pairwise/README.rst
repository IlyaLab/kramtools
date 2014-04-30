============================================================================
ABOUT
============================================================================

This is command-line software for performing (relatively) rapid pairwise 
analysis of *mixed-type* covariate data to identify *associations* between 
covariates.

The prototypical kind of association is (Pearson) correlation, but Pearson 
correlation is, in fact, only one of a wide variety of statistical tests
that can be used to establish that two variables ("features") are related.
This program implements a selection of tests and automatically chooses
the appropriate test for a given pair of covariates. It also carries out
a number of other processes on data to increase the efficacy of the
statistical tests.

This software has three high-level functions:

	1. (Optionally) select pairs of features (rows from the input matrix)
	   specified in one of several ways by the command line.
	   By default it iterates through *all* pairs.
	2. For every pair of features it chooses and executes a statistical test 
	   according to the features' statistical classes (boolean, categorical, 
	   ordinal, continuous).
	3. Report a configurable amount of information on each test in either
	   tabular or (soon) JSON format subject to configurable filters.
	   Output filters include the option of Benjamini-Hochberg FDR conrol.

Every run of pairwise involves these function, but *many* options control 
the exact behavior.

Input consists of textual matrices of tab-separated values in which
rows represent features and columns represent samples. See the libmtm
project in a neighboring directory tree for thorough documentation of 
input options and constraints as well as rationales behind the "statistical
classes" mentioned above.

**In the interest of minimizing stale and out-of-sync documentation functionality
that is thoroughly documented in the command-line executable's help text is not
repeated here.**

^^^^

============================================================================
STATISTICS
============================================================================

Three separate statistics are computed for each pair of features.
One of the following is calculated *between* the features using only 
those samples present in both features.

======================= ================================
Features                Statistic
======================= ================================
categorical/categorical Chisquare or Fisher Exact
categorical/numerical   Kruskal-Wallis
numerical/numerical     Spearman correlation
======================= ================================

Kruskal-Wallis is a non-parametric analog of 1-way ANOVA, 
which is approximately Chi-square distributed. The Spearman
correlation is the same as Pearson correlation but is calculated
on the ranks, not actual values, of the data.

Additionally, *within* each feature (F) the Kruskal-Wallis statistic is calculated 
between the samples of F present in both and samples present in F but missing in "the other."
For example, in the figure below where 'X' indicates missing data, 
the Kruskal-Wallis statistic is calculated in Feature 1 
between the light green and dark green samples and in Feature 2 between the light
blue and dark blue samples. The primary statistic (whichever of the above is used)
is calculated between the dark blue and dark green.

.. image:: ./doc/subsets.png

The K-W statistic *within* each feature provides some indication of whether 
the statistical result *between* the features was skewed by missing data.

The Spearman correlation is calculated for all pairs involving ordered
data. This obviously includes continuous pairs, but it also includes
categorical data with only two categories (i.e. boolean) since two categories
are always orderable (implicitly or otherwise). In other words, the Spearman
correlation is computed as long as one of the features is *not* strictly 
categorical with more than 2 categories. The absolute value of the Spearman 
correlation in some cases may not be meaningful, but its sign is.

The case where both features are categorical is slightly more
complex. In the interest of having statistically "good" tables
a preprocessing step is performed in which rows and/or columns
with empty cells are removed in such a way as to minimize the
total loss of samples. The approach, however, is greedy and not
guaranteed to be optimal. If after this "culling" step either of
the table's dimensions exceeds 2, a Chi-Square test is performed. 
Otherwise, if the table is 2x2, the Fisher exact test is performed.

^^^^

============================================================================
OUTPUT (FILTERING AND FORMATTING)
============================================================================

The *content* and *format* of output are *independently* configured.
Ouput is provided in one of two overarching formats: tabular and JSON_.
The "micro" formatting of tabular output is also tweakable.

.. _JSON : http://json.org

----------------------------
Content filtering
----------------------------

The following tables summarize the information available for every pair
of features analyzed. The *content* of the output can be filtered by
providing a format string that consists of the regular expressions listed
below separated by *any kind* of whitespace.

=========================== ==================================
Covariate information       Format regular expression
=========================== ==================================
sample count                c(ount)? [1]_
statistic name              st(atistic)?
test error                  e(rror)?[2]_
sign of correlation         si(gn)?("[3]_")?
statistic value             v(alue)?("[3]_")?
p-value                     p(robability)?("[3]_")?
-log10(p-value)             P(robability)?("[3]_")?
extra info                  e?x(tra)? [4]_
=========================== ==================================

.. [1] Count of the number of *pairs* in which neither value is missing.
.. [2] This is a hexadecimal value reporting the bit flags of the test.
       See online help.
.. [3] A printf-style format string valid for a floating-point value 
       delimited by double quotes.
.. [4] This is a semi-structures string consisting of any additional
       information made available by the statistical test (e.g. number of
       ties in rank data, or number of empty cells in a contingency table).

Each univariate feature is potentially subjected to a transformation
(e.g. outlier culling) as well as "internal" statistical tests as
described above. 

=========================== ================================
Univariate information      Format regular expression 
=========================== ================================
feature names               [<&>]f(eature)? [1]_
feature offsets             [<&>]o(ffset)? [1]_
feature statistical classes [<&>]cl(ass)? [1]_
transformation              [<&>]t(ransform)?
count of unused values      [<&>]u(nused)?
statistic                   [<&>]s(tatistic)?
statistic value             [<&>]v(alue)?
p-value                     [<&>]p(robability)?
extra info                  [<&>]e?x(tra)?
=========================== ================================

The left or right feature's univariate information is specified with "<"
and ">", respectively. "Left" and "right" corresponds to the ordering
of the covariate
"&" requests reporting of the information for both univariate features
as adjacent pairs.

----------------------------
False discovery rate control
----------------------------

This is implemented but untested.

^^^^

============================================================================
BUILDING
============================================================================

Dependencies:
	1. GNU Scientific Library (GSL_) 
	2. Lua_ is an optional dependency

Update the Makefile's GSLINC and GSLLIB variables with the location of
GSL's headers and libraries on your system.

Running::

	make
	
...in the src directory on a Linux machine with a suitable 
development environment builds the pairwise executable. 

.. _GSL: http://www.gnu.org/software/gsl
.. _Lua: http://www.lua.org

^^^^

============================================================================
TESTING
============================================================================

C code is structured as much as possible in unit-testable modules.
All unit testable C code has a conditionally compiled main() in its tail
to facilitate unit testing. See the relevant files and Makefile.
The unit tests are not automated; they are command line apps that are
intended for manual verification.

A blackbox testing framework is also provided in test/blackbox.
The exec.py script in test/blackbox:

	1. generates random data in R 
	2. analyzes the data in R
	3. compares the results of the two data paths.

It is run simply as::

	python3 exec.py <# of tests to perform>.

Note that for the case of categorical/categorical covariates
this involved a non-trivial reimplemention in R of the C++ "culling" code
in cat.cpp

^^^^

============================================================================
USAGE
============================================================================

Run the pairwise as::

	./pairwise yourdata.tab 

Many command line options are available. Run pairwise with no arguments.

^^^^

============================================================================
LIMITATIONS
============================================================================

1. No categorical feature may have more than 32 categories.

^^^^

============================================================================
OPEN ISSUES/TODO/WISHLIST
============================================================================

----------------------------------------------------------------------------
Reporting/filtering/counting of tests
----------------------------------------------------------------------------

Computation of statistics can fail for a variety of reasons related to
degeneracies in the input data. Some of these can be detected early (before
computation begins); some only become apparent in the coarse of computation.

Handling of degeneracies is furthermore bound up with three different
requirements that are somewhat at odds:

	1. maximizing useful output, "useful" being context-dependent.
	2. the need to filter output (i.e. to avoid a combinatoric explosion
	   of "uninteresting" results)
	3. the need to count *actually-performed* tests (for FDR control)

The typical sparsity of "interesting" results in the N-choose-2 possibilities
demands some in-program filtering (as opposed to simply piping the output
through a shell filter like awk), but there are multiple ways one
might filter. Only p-value filtering is currently available.

The need to support FDR control requires clear distinction between tests
*not performed* (e.g. because of early degeneracy detection) and failed
tests.

Proper handling of these issues is not fully resolved in the current 
implementation. 

Currently, NaNs are *intentionally* allowed to propagate to output; this is 
not a bug.

============================================================================
DESIGN
============================================================================

----------------------------------------------------------------------------
Goals
----------------------------------------------------------------------------

This software began as a reimplementation of an existing pipeline.
The requirement for compatibility with the prototype drove much
of its design.

The original program was motivated by one goal: speed. ..specifically, 
fast calculation of the several standard statistics describe above
on input with significant amounts of missing data. 

It was originally intended strictly for *exhaustive pairwise batch 
processing*.  Everything that deviates from this, e.g. row subselection, is 
an afterthought/add-on.

The goal of speed is approached in two ways:

	1. No memory allocation within loops; all memory is allocated before iteration commences.
	2. Implementation in a compiled language

Speed is no longer the *primary* driver of design. 

----------------------------------------------------------------------------
Degeneracy handling
----------------------------------------------------------------------------

Two types of degeneracies occur:
	1) those inherent in one (or both) *univariate* features
		a) categorical data with < 2 categories
		b) numerical data that is constant
	2) those that only emerge in the covariate pair because missing
	   data in one feature forces exclusion of values in the other

Univariate degeneracies are detected and flagged by the input parser (libmtm).

