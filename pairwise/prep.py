
"""
Converts a "TCGA format" tab-delimited ASCII feature matrix into
the binary format consumed by the pairwise executable.

Basically this just carries out a number of inspections, validations, 
and conversions on each row of the data that would otherwise 
1. be done redundantly in the C executable 
2. and would be more error-prone to implement in C.

The most important aspects of this format are:
1. Every value in the original matrix is converted to a 4-byte
   float or integer value in the result
2. Every row in the binary result is prefixed by a 32-bit word 
   that conveys attributes of that row to the pairwise 
   executable (the compiled C code, that is).

Row prefixes have the following format:

From the least to the most significant bit we have:
- The low-order 16 bits are the count of NAs
- The third byte is the category count for categorical variables.
- The high byte is unused

Data with a category count of 'one' is implicitly numeric.

Pictorially...
---------------------------------------
3322 2222 2222 1111 1111 1100 0000 0000
1098 7654 3210 9876 5432 1098 7654 3210
uuuu uuuD cccc cccC nnnn nnnn nnnn nnnn
---------------------------------------
u   unused
D   discrete flag
    1 => binary/categorical variable
	0 => numeric/scalar variable
c|C count of categories OR pseudo-count of numeric values
    if D==0, then
		C==0 => "OK"
		C==1 => numeric variable is a constant vector (excluding NA's)
    if D==1, then 'cccc cccC' is category count
n   count of missing values

Notice that an intentional implication of the above definitions is
that if the 'cccc cccC' field == 1, the feature is degenerate
regardless of its type.

Note that the column count refers to DATA columns, excluding the
(currently) 32-bit row prefixes. The row headers are excluded
because they might change size in the future in which case
treating the row headers as a column would invalidate calculations
in associated C code since byte widths might not agree.

The overall format of the resulting file is:
1) HEADER
2) MATRIX
3) STRINGS
4) ROWMAP
...and STRINGS and ROWMAP are optional

Just FYI, in the outer 'with' clause at the bottom, the overall
call structure of the program is:
	_verifyAbsence
	_process
	   _processNumeric
	   _processCategorical
	_appendStringTable (optionally)
"""

import sys
from struct import pack
import os.path
import subprocess
import optparse

assert sys.version_info.major >= 3

_SIG = "PAIRWISE"
_NAN_AS_UINT  = 0x7FC00000
_NAN_AS_FLOAT = float('nan')
_NAN_STRINGS  = frozenset(['na','nA','Na','NA'])
_SECTOR_SIZE  = 512
_MAX_LEVELS   = 32

parser = optparse.OptionParser( usage="usage: %prog [options]")
parser.add_option("-n","--no-names", action="store_false", dest="include_string_table",
		default=True,
		help="""Do not include row names in the binary file.
		Excluding row names precludes feature lookup by name.
		""")
parser.add_option("-H","--noheader", action="store_true", dest="noheader",
		default=False,
		help="""Do not expect a header line in the input file.
		""")
parser.add_option("-O","--overwrite", action="store_true", dest="overwrite",
		default=False,
		help="""Overwrite previous existing matrix/map files (if any).
		""")

parser.add_option("-C","--cat2num", action="store_true", dest="implicitNumeric",
		default=False,
		help="""Implicitly convert categorical features with >{} levels to numeric.
		""".format( _MAX_LEVELS ) )
def _verifyAbsence( fname ):
	if os.path.exists( fname ):
		print( "{} exists, and I won't overwrite it!".format( fname ) )
		sys.exit(-1)


def _processNumeric( fields, NC, fp, fnum, met ):
	FORMAT = "I{}f".format(NC)
	# Analyze the content of fields as strings before transforming it.
	# Regarding the non-missing values...
	PRESENT = [ float(v) for v in fields if v.upper() != "NA" ]
	# ...are they all identical, meaning degenerate? And if there's only 
	# 1 non-missing, yes, that is considered degenerate, too!
	isdegen = len(PRESENT) < 2 or all( [ PRESENT[0]==v for v in PRESENT[1:]] )
	NA = len(fields) - len(PRESENT)
	print( "NA", fnum, NA, sep="\t", file=met )
	f = [ _NAN_AS_FLOAT if v.upper()=="NA" else float(v) for v in fields ]
	b = pack( FORMAT, ( (0x00010000 if isdegen else 0) | (0x0000FFFF & NA) ), *f )
	if fp.write( b ) != (NC+1)*4:
		sys.exit( b )
	return len(f)


def _processCategorical( fields, NC, fp, fnum, met, convertCatToNum=False ):
	"""
	Be aware of categorical features that actually are not!...
	"""
	FORMAT = "I{}I".format(NC)
	NA = len([ 1 for v in fields if v.upper()=="NA" ])
	print( "NA", fnum, NA, sep="\t", file=met )
	S = set(fields).difference( _NAN_STRINGS )
	# Special-case genuinely boolean categories (without labels)...
	if S == frozenset(['0','1']):
		f = [ _NAN_AS_UINT if v.upper()=="NA" else int(v) for v in fields ]
		b = pack( FORMAT, ( 0x01020000 | (0x0000FFFF & NA) ), *f )
		if fp.write( b ) != (NC+1)*4:
			sys.exit( b )
		NF = len(f)
	elif len(S) <= _MAX_LEVELS:
		S = sorted(list(S))
		f = [ _NAN_AS_UINT if v.upper()=="NA" else S.index(v) for v in fields ]
		b = pack( FORMAT,
				( 0x01000000 | (0x00FF0000 & (len(S)<<16)) | (0x0000FFFF & NA) ), *f )
		if fp.write( b ) != (NC+1)*4:
			sys.exit( b )
		# Write out the encoding
		print( "FE\t{}\t{}".format( fnum, len(S) ), file=met )
		for i in range(len(S)):
			print( "CA\t{}\t{}".format(i,S[i]), file=met )
		NF = len(f)
	elif convertCatToNum:
		# There are too many categories, and implicit conversion has been
		# requested, so delegate this to _processNumeric with the data...
		if all([ "NA"==v.upper() or v.isdigit() for v in fields ]):
			# ...unchanged, if the labels are all integers (or NA), or...
			NF = _processNumeric( fields, NC, fp, fnum, met )
			conversionType = "natural"
		else:
			# ...converted to ARBITRARY integers if ANY of the labels 
			# are not already integers.
			S = sorted(list(S))
			NF = _processNumeric(
					[ "NA" if v.upper()=="NA" else str(S.index(v)) for v in fields ], 
					NC, fp, fnum, met )
			conversionType = "arbitrary"
			# ...and warn in metadata that implicit conversion has occurred.
		print( "IC\t{}\tcategorical -> numeric ({} order)".format( fnum, conversionType ), file=met )
	else:
		raise RuntimeError(
			"error: factor with too many levels ({}).\n\tMaybe try the --cat2num option?".format(len(S)) )
	return NF


def _process( idat, odat, ometa, expectheader=True, implicitNumeric=False ):
	"""
	Binarized matrix is written to odat.
	Textual category encodings are written to fp_enc.
	Returns number of DATA rows and DATA columns.
	"""
	loff = 0 # Feature offset and line offset need not correspond...
	feat = 0 # ...because of comment lines.
	rows = []
	NC = 0
	# Assume the first line is a header, but check for expected content.
	if expectheader:
		line = idat.readline()
		fields = line.strip().split('\t') # both leading AND trailing ws.
		NC = len(fields)-1
		loff += 1
	# Remaining lines are data.
	line = idat.readline()
	while len(line) > 0:
		fields = line.strip().split('\t') # both leading AND trailing ws.
		ty = line[0].upper()
		rowname = fields[0]  # Strip the TCGA label off the front...
		datcols = fields[1:] # ...so NC remainder are data.
		nc = len(datcols)
		if NC == 0:
			NC = nc # if a header missing, 1st data line sets column count.
		elif NC != nc: # if column count doesn't match that from header.
			raise RuntimeError(
				"{} data columns on line {}, expected".format(
				nc, loff+1, NC ) )
		if 'N' == ty: # check most common first!
			_processNumeric(     datcols, nc, odat, feat, ometa )
		elif 'B' == ty or 'C' == ty:
			try:
				_processCategorical( datcols, nc, odat, feat, ometa, implicitNumeric )
			except Exception as x:
				raise RuntimeError( "{} in row {}".format( x, rowname ) )
		else:
			raise RuntimeError("unexpected type at line "+str(loff+1) )
		rows.append( [ rowname, feat ] )
		feat += 1
		loff += 1
		line = idat.readline()
	return ( nc, rows )


def _appendStringTable( strings, fp ):
	"""
	This appends 
	1) a packed string table (ragged array) sorted by string value and
	2) a map of (string_offset,matrix_row_offset) pairs
	...to the end of the file.
	"""
	# Pad the file to a 16-byte boundary (just to simplify debugging in
	# hexdump).
	NUL = bytes(1)
	while fp.tell() % 16:
		fp.write( NUL )
	# This is where the string section will reside, so note the offset...
	OFF_STR = fp.tell()
	# Store the strings, noting their offsets FROM THE BASE of the section
	# and after each is stored, REPLACE it in the pair with its offset.
	ordered = sorted( strings, key=lambda x:x[0] )
	for pair in ordered:
		OFF = fp.tell() - OFF_STR
		s = pair[0].encode('ascii')
		fp.write( pack( "{}sc".format(len(s)), s, NUL ) )
		pair[0] = OFF # REPLACE the string with its offset.
	# ...pad to a disk sector size boundary. (See notes at top for why.)
	while fp.tell() % _SECTOR_SIZE:
		fp.write( NUL )
	# This is where the row map will reside, so note its file offset.
	OFF_MAP = fp.tell()
	for pair in ordered:
		fp.write( pack( "II", pair[0], pair[1] ) ) # string offset, matrix row
	return ( OFF_STR, OFF_MAP )

############################################################################

opts,args = parser.parse_args()

INPUT  = args[0]
OUTPUT = os.path.splitext(INPUT)[0] + ".bin"
CATEGS = os.path.splitext(INPUT)[0] + ".cat"

if not opts.overwrite:
	_verifyAbsence( OUTPUT )
	_verifyAbsence( CATEGS )

with open( INPUT ) as fp_in:
	# Calculate the input matrix' MD5 hash to provide a hard tie to the
	# resulting output matrix' source.
	p = subprocess.Popen( [ "md5sum", INPUT ], stdout=subprocess.PIPE )
	if p:
		md5 = p.communicate()[0][:32].decode()
	else:
		raise RuntimeError("Is md5sum on your system?")
	try:
		with open( OUTPUT, "wb" ) as fp_out:
			# Write out the header used by the C allpairs. Most of the values
			# are unknown at this point, so we write 0's as placeholders then
			# backup and fill them in at the very end when values are known.
			SIGNATURE_BYTES = _SIG.encode('ascii')
			HDR_SIZE = fp_out.write( pack(
				"8sIIIIII32s",
				SIGNATURE_BYTES,
				0, # header size
				0, # rows
				0, # columns (of data)
				0, # Offset to strings (0 => no row map/string table)
				0, # Offset to row map (0 => no row map/string table)
				0, # Unused
				md5.encode('ascii') ) )
			assert HDR_SIZE == 64
			# Open another file to receive (textual) metadata.
			with open( CATEGS, "w" ) as fp_cat:
				print( "# for matrix {} with MD5 {}".format( INPUT, md5 ), file=fp_cat )
				# ...and here begins the processing in earnest.
				data_column_count,row_names = \
						_process( fp_in, fp_out, fp_cat, not opts.noheader, opts.implicitNumeric )
			if opts.include_string_table:
				off_str, off_map = _appendStringTable( row_names, fp_out )
			else:
				off_str, off_map = 0,0
			if off_str > 0xFFFFFFFF or off_map > 0xFFFFFFFF:
				raise RuntimeError( "File too large" )
			# Now back up to header and fill in header size, rows, and columns.
			fp_out.seek( len(SIGNATURE_BYTES) )
			fp_out.write( pack( "IIIII",
				HDR_SIZE,
				len(row_names),
				data_column_count,
				off_str,
				off_map ) )
	except Exception as x:
		print( "Failed processing file: {}".format(INPUT), file=sys.stderr )
		print( sys.exc_info() )
		# In the event of failure, don't leave ANYTHING behind that might
		# allow a parent script to continue.
		if os.path.isfile( OUTPUT ):
			print( "Removing", OUTPUT, file=sys.stderr )
			os.remove( OUTPUT )
		if os.path.isfile( CATEGS ):
			print( "Removing", CATEGS, file=sys.stderr )
			os.remove( CATEGS )

if __debug__:
	print( HDR_SIZE, "header bytes" ) # just informational/sanity check

