#!/usr/bin/awk -f

# This script finds all line numbers matching a given regular expression
# Usage: awk -f scan_line_idxs.awk "pattern" filename

# BEGIN block to process the command line arguments
BEGIN {
    if (ARGC != 3) {
        print "Usage: " ARGV[0] " \"pattern\" filename"
        exit 1
    }
    pattern = ARGV[1]
    filename = ARGV[2]
    ARGV[1] = filename  # Make the filename the second argument for awk to process
}

# Main block to process each line in the file
{
    if ($0 ~ pattern) {
        print FNR  # Print the line number if it matches the pattern
    }
}