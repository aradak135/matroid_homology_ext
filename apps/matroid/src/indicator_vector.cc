#include "polymake/client.h"
#include "polymake/Vector.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/graph/Lattice.h"
#include "polymake/graph/Decoration.h"
#include "polymake/internal/sparse2d.h"
#include <list>
#include <algorithm>
#include <stdio.h>
#include <iostream>
#include <bitset>
#include <math.h>
#include <vector>
#include <typeinfo>

namespace polymake { namespace matroid {

Vector<Int> indicator_vector_full(BigObject m, Array<Matrix<Int>> allChains, Array<Matrix<Int>> matroidChains) {
	using graph::Lattice;
	using graph::lattice::BasicDecoration;

	int dim = (int)allChains.size();
	Vector<Int> indicatorVector(dim, 0);

	// Lookup the index of each matroid chain in the global list and create a [0,1] incidence array.
	// One perforamnce improvement might be to stash the indexes in a cache for reuse across multiple matroids. 
	for (Matrix<Int> chain : matroidChains) {
		// Use find to get the index of the found chain (scans the full list)
		Int idx = std::find(allChains.begin(), allChains.end(), chain) - allChains.begin();
		indicatorVector [idx] = 1;
	}
	
	return indicatorVector;
}

Vector<Int> indicator_vector(BigObject m){
	Int n = m.give("N_ELEMENTS");
	Int r = m.give("RANK");
	Array<Matrix<Int>> allChains = call_function("generate_chain_matrices",r,n);
	Array<Matrix<Int>> matroidChains = call_function("chains_of_flats",m);

	return indicator_vector_full(m,allChains,matroidChains);
}

UserFunctionTemplate4perl("# @category Other"
                  "# Generates an array of chains representing a depth first search across the power set 2^n."
                  "# @param Matroid m the original matroid"
                  "# @return Vector<Int> an indicator vector [0/1] noting which of allChains are included in the matroid",
                  "indicator_vector(Matroid)");


UserFunctionTemplate4perl("# @category Other"
                  "# Generates an array of chains representing a depth first search across the power set 2^n."
                  "# @param Matroid m the original matroid (optional if allChains/matroidChains are provided)"
                  "# @param Array<Matrix<Int>> allChains the full list of chains available to the matroid"
                  "# @param Array<Matrix<Int>> matroidChains the list of maximal chains which define the matroid"
                  "# @return Vector<Int> an indicator vector [0/1] noting which of allChains are included in the matroid",
                  "indicator_vector_full(Matroid, Array<Matrix<Int>>, Array<Matrix<Int>>)");
} }
