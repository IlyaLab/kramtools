
# Turn any line containing "USAGE" into a C variable definition.
# There should be only 1, and it should be the first line.
/^USAGE/ {
s/$/ =/
s/^/const char */
n
}

# Escape any and all EXISTING quotes
s/"/\\"/g
# Convert any and all EXISTING to fixed whitespace.
s/\t/    /g

# Surround every line NOT prefixed with # or U with quotes and an escaped newline.
/^[^#]/s/^/"/
/^[^#]/s/$/\\n"/

# Last line only, append a semicolon, too.
$s/$/;/

# Insure empty lines are reflected in C string result
/^$/c\
"\\n"
