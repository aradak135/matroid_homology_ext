#!/usr/bin/env sh
# First time setup script to import the extension and generate resource files

echo "Importing matroid_homology_ext to polymake...";
BASEDIR=$(dirname "$0")
echo "obliterate_extension(\"$BASEDIR\");exit;" | polymake -
ninja -C $BASEDIR/build/Opt
echo "import_extension(\"$BASEDIR\");exit;" | polymake -

echo "Generating chainfiles (this may take a couple minutes)...";
mkdir -p "$BASEDIR/chainfiles"
for n in $(seq 1 8)
do
	for r in $(seq 1 $n)
	do
		polymake --script "$BASEDIR/scripts/generate_chain_files.pl" "$BASEDIR/chainfiles/chainIndices_r${r}_n$n.txt" $r $n
	done
done

echo "Done! To try out the extension run:"
echo "\tpolymake --script \"$BASEDIR/scripts/matroid_homology_fast.pl\""
