#!/usr/bin/env sh
# Repeatable reload script to update the extension (e.g. with newer code)

echo "Destroying and rebuilding matroid_homology_ext in polymake...";
BASEDIR=$(dirname "$0")

ninja -C $BASEDIR/build/Opt
echo "reconfigure_extension(\"$BASEDIR\");exit;" | polymake -

echo "Done! To try out the extension run:"
echo "\tpolymake --script \"$BASEDIR/scripts/matroid_homology_fast.pl\""
