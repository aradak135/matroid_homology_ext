# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

This is a working repository hosting the polymake extension used to compute the matroid homology introduced by Simon Hampe. Applications include the computing of a single indicator vector from a lone matroid and computing basic properties of the intersection ring of matroids endowed with either the deletion or contraction induced boundary operator.

### How do I get set up? ###

*Dependency: you must be running a Linux machine with Polymake 4.X or later installed.*

First, clone (or simply download and unzip) this repository to a folder on your machine.
```
git clone 

https://github.com/aradak135/matroid_homology_ext.git
```

Then, invoke "setup.sh" to load the Polymake extension and populate resource (chain files). Chain files are included with the download, but may be omitted.
```
cd matroid_homology_ext && sh setup.sh
```

When that completes, you're done! You can run any of the scripts built on the extension, for instance:
```
polymake --script "scripts/matroid_homology_fast.pl"
```

### Who do I talk to? ###

* Austin Alderete aaldere@math.utexas.edu