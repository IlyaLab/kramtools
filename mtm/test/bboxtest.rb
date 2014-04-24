
require 'set'

=begin ---------------------------------------------------------------------

Testing in two general modes:
1. random generation of all valid formats
   Ideally variation on valid formats ought to be generated.
2. targeted insertion of point flaws

Header in wrong (not first) row generally cannot be detected except that
it might appear as a categorical variable with too many categories and
it might have empty first field (an aberration that is only permitted to
the header).

Capabilities requiring testing:
1. type inference
2. correct parsing

Independent options
1. header
2. row labels
   2a. with type info
   2b. without type info
3. All valid representations of all class types
   3a. boolean as
       string
       0,1
       +,1
       true, false
   3b. categorical as
       random strings
       integers
   3c. ordinal
   3d. floating point
       decimal
       exponential
       both signs
       integral early values

Include aberration
1. wrong column count in a row
2. empty cells (adjacent separators)
3. incompatible mixed token types (e.g. string, float)
4. empty cells at end-of-line
=end


DEFAULT_PROB_NA=0.1
SINGLE_PREC_MAX=3.4e38
ABS_TEST_BOUNDS=1e9
MAX_CATEGORY_COUNT=32

OUTPUT_FLOAT_FORMAT = "%.1e" # same as the executable under test
OUTPUT_NA_MARKER    = "NA"

def randomAsciiString( len )
	len = rand(64) + 1 if len == 0
	(0x21..126).to_a.sample(len).collect {|i|i.chr}.join("")
end

def randomAlphaString( len )
	len = rand(64) + 1 if len == 0
	((65...91).to_a + (97...123).to_a).sample(len).collect {|i|i.chr}.join("")
end

def randomChoice()
	return rand() < 0.5
end

=begin ---------------------------------------------------------------------
Each subclass of feature needs to be able to emit two different 
representations of itself:
1. random formats with occasional violations of formatting rules
2. deterministic format mirroring what mtproc is expected to emit.

The randomly formatted version is consumed by mtm_parse, and the 
deterministically formatted version is used for round-trip testing.
=end

class Feature
	def initialize(args={})
		# @na controls whether missing data is represented explicitly or 
		# implicitly ("implicit" meaning adjacent separators with nothing
		# between them).
		@na   = args.has_key?(:na) ? args[:na] : ( randomChoice() ? "NA" : "" )
		@name = args[:name] # may be nil
	end

	def degenerate?
		( (@vector.length - countMissing() ) < 2) or Set.new(@vector.select {|x| not x.nil? }).length == 1
	end

	def countMissing
		@vector.count {|x| x.nil? }
	end
end

=begin ---------------------------------------------------------------------
Each continuous feature generates its own random float format.
=end
class RandomContinuousFeature < Feature

	def initialize( n, min, max, args={} )
		super(args)
		# Create a random floating point format...
		@precision = 3+rand(6)
		@notation  = randomChoice() ? 'f' : 'e'
		@flag      = randomChoice() ? ""  : "+"
		@format = "%%%s.%d%c" % [ @flag, @precision, @notation ]
		# ...and a random floating-point array with a random amount of 
		# missing data.
		p = args.has_key?(:prob_na) ? args[:prob_na].to_f : DEFAULT_PROB_NA
		@vector = Array.new(n) do |i|
			rand() < p ?  nil : min + rand()*(max-min)
		end
	end

	def printTestInput( fp )
		fp.printf( "%s\t", @name ) if @name
		fp.puts( @vector.collect { |f| f.nil? ? @na : @format % f }.join("\t") )
	end

	#                         !!! WARNING !!!
	#     Parser can't produce more precision than was in its input.
	#                      Format must respect this.
	#                         !!! WARNING !!!


	def rnd( val )
		raise 'Trouble' if val == 0.0
		p = 10**(Math.log10(val.abs).floor)
		p = val < 0 ? -p : p
		return (val / p).round(1) * p
	end

	def printTestExpectedOutput( fp, format )
		fp.printf( "%s\t", @name ) if @name
		missing = countMissing()
		fp.printf( "%s:%s:%d:%d\t",
				missing < @vector.length ? "F" : "?",
				degenerate? ? "!" : "-", 
				0, 
				missing )
		fp.puts( @vector.collect { |f| f.nil? ? OUTPUT_NA_MARKER : format % f }.join("\t") )
	end
end # RandomContinuousFeature

=begin ---------------------------------------------------------------------
Each categorical feature generates random labels
=end
class RandomCategoricalFeature < Feature

	def initialize( n, k, args={} )
		super(args)
		largest = 0
		p = args.has_key?(:prob_na) ? args[:prob_na].to_f : DEFAULT_PROB_NA
		@vector = Array.new(n) do |i|
			r = nil
			if rand() > p 
				# This creates a random array of [0..k-1] such that 
				# integers' first occurrences are ordered. This 
				# mirrors the way strset.c maps symbols to integers.
				r = rand(largest+1)
				if r == largest and largest+1 < k then
					largest += 1
				end
			end
			r
		end
		@cats = Set.new(@vector.select {|x| not x.nil? }).length 
		# Now that we know how many categories -actually- exist in the
		# feature (which will be <= k), create some labels for them.
		@labels = Array.new(@cats) {|i| randomAlphaString(args[:l])}
	end

	def printTestInput( fp )
		fp.printf( "%s\t", @name ) if @name
		fp.puts( @vector.collect { |x| x.nil? ? @na : @labels[x] }.join("\t") )
	end

	# Categorical values become integers in the output.
	def printTestExpectedOutput( fp )
		fp.printf( "%s\t", @name ) if @name
		missing = countMissing()
		fp.printf( "%s:%s:%d:%d\t",
				missing < @vector.length ? "I" : "?",
				degenerate? ? "!" : "-", 
				@cats, 
				missing )
		fp.puts( @vector.collect { |x| x.nil? ? OUTPUT_NA_MARKER : x  }.join("\t") )
	end
end # RandomCategoricalFeature


=begin ---------------------------------------------------------------------
 Generate the matrix then emit it in two forms, one for mtproc's consumption
 and the other should be the output mtproc will produce from round-tripping
 the matrix.
=end
class RandomMatrix

	def initialize( nrows, ncols )
		@has_row_names = randomChoice()
		# Use 'H' suffix just to (visually) distinguish header from possible
		# first line categorical.
		@col_names = randomChoice() ? Array.new(ncols) {|c| randomAlphaString(5)+'H' } : nil
		@rows = Array.new( nrows ) do |r|
			# Insure that no row label starts with a comment prefix ('#')!
			# Keep row labels short to simplify running manual unit tests
			# that refer to rows by name.
			rowlabel = @has_row_names ? randomAlphaString(5) : nil
			# Randomly produce either continuous or categorical features.
			if randomChoice() then
				RandomContinuousFeature.new( 
					ncols, -ABS_TEST_BOUNDS, ABS_TEST_BOUNDS, name:rowlabel )
			else
				RandomCategoricalFeature.new(
					ncols, 1+rand(MAX_CATEGORY_COUNT), name:rowlabel, l:5 )
			end
		end
	end

	# Print the matrix in the format for the mtm_parse API's consumption,
	# the format that serves as -input- to mtm_parse
	def printTestInput( fp )
		if @col_names then
			# Need a place holder (or empty field) 
			if @has_row_names then
				fp.printf( "%s\t%s", randomChoice() ? "X" : "", @col_names[0] )
			else
				fp.print( @col_names[0] )
			end
			@col_names.each_index do |i|
				fp.printf( "\t%s", @col_names[i] ) if i > 0
			end
			fp.print( "\n" )
		end
		@rows.each do |row|
			row.printTestInput( fp )
		end
	end

	# Print the matrix in the format in which mtproc is expected to echo
	# a parsed matrix.
	def printTestExpectedOutput( fp )
		# mtproc does not re-emit the header (column names)
		@rows.each do |row|
			if row.instance_of?( RandomContinuousFeature ) then
				row.printTestExpectedOutput( fp, OUTPUT_FLOAT_FORMAT )
			else
				row.printTestExpectedOutput( fp )
			end
		end
	end

	def has_col_names
		not @col_names.nil?
	end

	attr_reader :has_row_names
end

=begin ---------------------------------------------------------------------
If the first token in the first argument refers to an executable file, run 
in batch test mode.  Otherwise, just generate a pair of matrices to stdout.
The remaining tokens in the first argument are assumed to be switches that
are passes verbatim to the program.
=end

DEFAULT_INPUT="i"
DEFAULT_EXPECT="e"
DEFAULT_OBSERVE="o"

if File.executable?( ARGV[0].split(' ')[0] ) then
	n = ARGV.length > 3 ? ARGV[3].to_i : 1
	executed = 0
	while n > 0 do
		n = n - 1
		m = RandomMatrix.new( ARGV[1].to_i, ARGV[2].to_i )
		File.open( DEFAULT_INPUT, "w") do |fp|
			m.printTestInput( fp )
		end
		File.open( DEFAULT_EXPECT, "w") do |fp|
			m.printTestExpectedOutput( fp )
		end
		# Build a command line to run the matrix parser on the generated input
		# and re-emit it to a file that will constitute the "observed."
		cmd = sprintf( "%s %s %s %s > %s", 
			ARGV[0], 
			m.has_col_names ? "" : "-h",
			m.has_row_names ? "" : "-r",
			DEFAULT_INPUT,
			DEFAULT_OBSERVE )
		break if not system( cmd )
		if not system( "diff e o" ) then
			STDOUT.printf( "-------------------\n" )
			system('cat i')
			STDOUT.printf( "-------------------\n" )
			STDOUT.printf( "%d tests passed\n", executed )
			break
		end
		executed += 1
	end
else
	m = RandomMatrix.new( ARGV[0].to_i, ARGV[1].to_i )
	m.printTestInput( STDOUT )
	STDOUT.puts ""
	m.printTestExpectedOutput( STDOUT )
end

