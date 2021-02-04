# this awkwardly converts a tree in NCBI's "new" text format to "old" text format, for backward-compatibility with other code
# all of this is just for visualisation and it is recommended to use a neater pipeline with e.g. phytools or ETE3

sed -E 's/\+\-/\+ /g' $1 | sed -E 's/\|/\+/g' | sed -E 's/\\\-/\+ /g' | sed -E 's/^ /\+/g' | sed -E 's/\+  /\+ \+/g' | sed -E 's/\+  /\+ \+/g' | sed -E 's/\+  /\+ \+/g' | sed -E 's/\+  /\+ \+/g' | sed -E 's/\+  /\+ \+/g' | sed -E 's/\+  /\+ \+/g' |  sed -E 's/\+  /\+ \+/g' | sed -E 's/\+\+/\+ /g' | sed -E 's/ \\\+/ \+ /g'
