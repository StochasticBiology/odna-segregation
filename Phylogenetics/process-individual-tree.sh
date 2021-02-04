./tree-to-old-format.sh $1 > $1-old.txt
./pad-tree.sh $1-old.txt $2
gcc translate-tree.c -o translate-tree.ce
./translate-tree.ce $1-old.txt.padded
