
=begin
This script verifies that pairwise' results are independent of the order
of feature pairs.
It generates random feature pairs (random in type and content) and executes
pairwise twice on each pair, reversing the order of the pairs.

If the computed p-values are identical it emits nothing.
If the computed p-values are different it emits the offending 2-row table
to stdout, and halts testing immediately.
=end

require 'set'
require 'tempfile'

DEFAULT_COLUMN_COUNT = 100
DEFAULT_PROB_NA      = 0.1
ABS_TEST_BOUNDS      = 1e9
MAX_CATEGORY_COUNT   = 32
OUTPUT_NA_MARKER     = "NA"

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

=begin
=end

if ARGV.length < 1 or /--?h(elp)?/.match( ARGV[0] ) then
	printf( "commute.rb <executable> [ column_count [ test_count ] ]\n" )
	exit(0)
end

executable = ARGV[0]
ncols      = ARGV.length > 1 ? ARGV[1].to_i : DEFAULT_COLUMN_COUNT
ntests     = ARGV.length > 2 ? ARGV[2].to_i : 1

while ntests > 0 do
	# Generate a pair of random features
	feature = Array.new( 2 ) do |r|
		if randomChoice() then
			RandomContinuousFeature.new( 
				ncols, -ABS_TEST_BOUNDS, ABS_TEST_BOUNDS )
		else
			RandomCategoricalFeature.new(
				ncols, 1+rand(MAX_CATEGORY_COUNT), l:5 )
		end
	end
	input = Tempfile.new('pairwise')
	begin
		feature[0].printTestInput( input )
		feature[1].printTestInput( input )
		input.close
		results = []
		# Execute pairwise on them both directions
		IO.popen( "#{executable} -v 0 -h -r -P 0,1 #{input.path}", 'r' ) { |fd|
			results << fd.gets.chomp.split("\t")
		}	
		IO.popen( "#{executable} -v 0 -h -r -P 1,0 #{input.path}", 'r' ) { |fd|
			results << fd.gets.chomp.split("\t")
		}
		if (results[0][4].to_i != results[1][4].to_i) or (results[0][5].to_f != results[1][5].to_f) then
			system( "cat #{input.path}" )
			break
		end
	ensure
		input.unlink
	end
	# Verify that the results are identical
	ntests -= 1
end

