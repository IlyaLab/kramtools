
require 'bigdecimal'

=begin
0. Repeat
1. Generate random data with given attributes.
2. Feed it to ut-feature, and
3. Confirm ut-feature's output is correct.

Priorities:
	1. Confirm feature.c correctly handles all correct data.
		?? How to confirm random generation covers all valid
		formats?
		Independently specifiable/encodeable data elements:
		(prefix)?label data data data ...
		If prefix is present
			Constraints exist on data.
		else
			Data is interpreted in most restrictive way possible
			but parsing effectively cannot fail.
	2. Confirm feature.c is non-silent (where possible) on
		incorrect input and fails gracefully in all failure cases.
Corner cases:
	1. Hit the three types of degeneracy
		a. entirely missing data
		b. categorical data with only one category represented (in
			multiple values)
		c. continuous data with only one non-missing value
	2. boolean data with >2 categories
=end

=begin
:class one-of { :bool, :cat, :ord, :cont }
:rep
	if cat/bool one-of { str, int }
	if cont one-of { int, fp }
:card
	if cat random int
	if bool 2
	if cont nil
:labeled {T,F}
:prefixed {T,F}
:length 
:missing
:na []
=end

STATCLASS = [ :bool, :cat, :cont ] # "ordinal"
REP_CAT   = [ :str, :int ]
REP_NUM   = [ :fp, :int ]
MAXLEN    = 100
MAXLABLEN = 32
MAXCARD   = 32
MAXFPVAL  = 2**23 # Within
MAXINTVAL = 2**23

def prefix_for( sym )
	case sym
	when :bool then "B:"
	when :cat  then "C:"
	when :cont then "N:"
	else throw "invalid statistical class: " + sym.id2name
	end
end

# Generate a random string of *non-control character* ASCII text
# of a specified length.
def random_ascii( len )
	len = rand(64) + 1 if len == 0
	(0x21..126).to_a.sample(len).collect {|i|i.chr}.join("")
end

# Generate a random string of *strictly alphabetic*, upper- and
# lower-case ASCII text of a specified length.
def random_alpha( len )
	len = rand(64) + 1 if len == 0
	((65...91).to_a + (97...123).to_a).sample(len).collect {|i|i.chr}.join("")
end

def random_float( prec=6 )
	digits = (0x30..0x39).to_a.sample(prec).collect {|i|i.chr}.join("")
	return digits.insert( rand(digits.length+1),'.')
end

# Generate an ordered random sequence of n integers 0...k.
# The result might not exhaust [0,k-1]
def random_ordered_array( n, k )
	a = Array.new(n) {0}
	for i in 1...n do
		cmax = a[0...i].max
		a[i] = rand( [ cmax+2, k].min )
	end
	return a
end

# Generate a random permutation of the contiguous integers 0...n,
# optionally holding the first k fixed. 
def random_permutation( n, fix=0 )
	v = (0...n).to_a
	for i in fix...n do
		# j <= random from [fix,i]
		j = fix + rand(i-fix+1)
		v[i] = v[j]
		v[j] = i
	end
	return v
end


# Create a random data specification, or finish filling out an existing but
# incomplete one. By default, this generates "sensible" specs--that is,
# specs that are internally consistent and represent generally "good" data.
# For example, a spec with boolean class will have cardinality 2;
# categorical class will have cardinality >2.
# ALL attributes will be defined one way or another when this method
# returns.
# 
# For error injection, the caller can pre-specify certain attributes, and
# these will be left intact. This allows forcing "bad" specs, in particular,
# to generate data that will test corner cases.
def random_spec( spec = Hash.new  )
	spec[ :class ] = STATCLASS[ rand( STATCLASS.length ) ] \
		if not spec.has_key?( :class )
	spec[ :rep ] = ( spec[ :class ] == :cont \
			? REP_NUM[ rand( REP_NUM.length ) ] \
			: REP_CAT[ rand( REP_CAT.length ) ] ) \
		if not spec.has_key?( :rep )
	spec[ :labeled ] = rand() > 0.5 ? TRUE : FALSE \
		if not spec.has_key?( :labeled )
	spec[ :prefixed ] = spec[ :labeled ] and rand() > 0.5 ? TRUE : FALSE \
		if not spec.has_key?( :prefixed )
	spec[ :length ] = rand( MAXLEN ) \
		if not spec.has_key?( :length )
	spec[ :missing ] = rand( spec[ :length ]/2 ) \
		if not spec.has_key?( :missing )
	spec[ :na ] = ["NA",""] \
		if not spec.has_key?( :na )
	# Choose cardinality sensibly according to :class and number of
	# of non-missing data values.
	if not spec.has_key?( :card ) then
		if spec[ :class ] == :cont then
			card = nil # TODO: Could be overridden in random_data?
		elsif spec[ :class ] == :bool then
			card = 2
		else
			card = rand( MAXCARD-3 )+3
		end
		spec[ :card ] = card
	end
	return spec
end


# Create a random Array according to a "spec" expressed in a Hash.
def random_data( spec )
	n = spec[:length] - spec[:missing]
	integral = (spec[:rep] == :int)
	a = case spec[:class]
	when :bool,:cat then \
		random_ordered_array( n, spec[ :card ] )
	when :cont then Array.new(n) { |i|
			integral ? (MAXINTVAL*rand()).round() : random_float.to_f
		}
	end
	# Finally insert missing values (nil for now)
	missing = spec[:missing]
	# Now insert missing data at random positions (preserving the
	# existing ordering of data).
	while missing > 0 do
		a = a.insert( rand( a.length+1 ), nil )
		missing -= 1
	end
	throw "Crap" if a.length != spec[:length]
	return a
end


# Transform a data vector by
# 1. Converting nil's to one of allowed missing data representation(s)
# 2. Mapping integers to string labels, 
# Combine a data vector with the prefix/label given in a feature
# spec to produce a single string constituting one *line* of input.
def feature_string( spec, data )
	na = spec[ :na ]
 	# Force na to be an array to simplify below.
	na = [na,] if not na.is_a?(Array)
	# Convert missing data (nil's) to output format.
	if spec[:class] == :cat then
		# Determine the largest integer
		m = data.reject {|v| v.nil? }.max
		labels = Array.new(m+1) {|i| random_alpha(5) }
		d = data.map {|v| v.nil? ? na[rand(na.length)] : labels[v] }
	else
		d = data.map {|v| v.nil? ? na[rand(na.length)] : v }
	end
	# Create and prepend a (prefix)?label string.
	head = ""
	if spec[:labeled] then
		head = random_alpha(rand(MAXLABLEN)+1)
		if spec[:prefixed] then
			head = prefix_for( spec[:class] ) + head
		end
		d.insert(0,head)
	end
	return d.join("\t") + "\n"
end

def exec( spec, idata )
	# Convert the data into single string which might constitute one
	# line of a matrix/table of data...
	input = feature_string( spec, idata )
	p input
	# ...and send it to the unit test.
	cmd = ['./ut-feature']
	if not spec[:labeled] then
		cmd << "-r"
	end
	output = IO.popen( cmd, "r+" ) { |f|
		f.write( input )
		f.readlines
	}
	# Output now contains all lines emitted from ut-feature.
	# Convert them to a Hash containing both the data and metadata.
	result = { :odata => []}
	for l in output do
		l.chomp!
		if l[0] == '#' then
			# Turn the metadata into hash key:value pairs...
			pair = l.slice(1,l.length).split(':',2)
			result[ pair[0].to_sym ] = eval( pair[1] )
		else
			# ...and append data to the array.
			result[:odata] << eval(l)
		end
	end
	# Count the number of data values that did NOT survive the round trip.
	d = result[:odata]
	if result[:int]
		return (0...idata.length).count { |i| idata[i] != d[i] }
	else # Allow for rounding error in feature.
		count = 0
		for i in 0...idata.length
			if (idata[i].nil? != d[i].nil?) and (idata[i] - d[i]).abs/idata[i] > 1e-7 then
				puts "#{idata[i]} <> #{d[i]}"
				count += 1
			end
		end
		return count
	end
end

# Verify that the unit test code for feature.c output
# 1. a correct interpretation of the input and
#	1) class
#	2) class specified/inferred
#	2) representation
#	3) actual cardinality
#	4) actual length
#	5) missing count
#	6) label
# 2. all non-missing data output is identical to input (within bounds
#    of floating-point precision)
def verify( input, output )
	# Class was interpreted correctly
	# Representation was interpreted correctly.
end

s = {:length => 10, :class => :cont, :rep => :fp, :missing => 0 }
s = random_spec({:length => 10 })
d = random_data(s)

p s
p exec(s, d)
#p random_float.to_f

