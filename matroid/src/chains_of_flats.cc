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
#include<iostream>
#include<bitset>
#include <math.h>
#include <vector>
#include <typeinfo>

namespace polymake { namespace matroid {


// Custom comparator for comparing maximal chain indexes by order (# of elements)
struct lattice_face_idx_comparator
{
	// Constructor
	std::vector<std::vector<Int>> _facesOfLattice;
	lattice_face_idx_comparator(std::vector<std::vector<Int>> facesOfLattice) {
		_facesOfLattice = facesOfLattice;
	}
	
	// Comparator
    inline bool operator() (Int idx1, Int idx2) {
        return _facesOfLattice[idx1].size() < _facesOfLattice[idx2].size();
    }
};


std::vector<Matrix<Int>> chains_of_flats(BigObject m) {
	using graph::Lattice;
	using graph::lattice::BasicDecoration;
	
	// We will a vector of (r+1) x (n) matrices, each one representing a chain
	std::vector<Matrix<Int>> listOfChainMatrices;
	int n = m.give("N_ELEMENTS");
	int r = m.give("RANK");

	// LATTICE_OF_FLATS->FACES returns all possible chain matrix rows (chain elements) in the matroid
	// e.g. facesOfLattice[3] might return {2} corresponding to a chainMatrix row of [0,0,1,0]
	BigObject latticeOfFlats = m.give("LATTICE_OF_FLATS");
	std::vector<std::vector<Int>> facesOfLattice = latticeOfFlats.give("FACES");

	// maximal_chains_of_lattice returns a list of maximal chains which define the matroid
	// e.g. maxChainIndexes.row(0) might return {0,2,5} indicating a chain of latticeOfFlats[0], latticeOfFlats[2], and latticeOfFlats[5]
	IncidenceMatrix<NonSymmetric> maxChainIndices = call_function("maximal_chains_of_lattice", latticeOfFlats);
	
	for (int chainNum = 0; chainNum < maxChainIndices.rows(); ++chainNum) {
		Set<Int> chain = maxChainIndices.row(chainNum);
		Matrix<Int> chainMatrix(r+1, n);

		// lexicological order for chain sets is sometimes ascending (0000 -> 0100 --> 1111) and other times descending (1111 -> 0100 --> 0000).
		// Sets are orderless so it could just as easily be randomly orderered. So to be safe, sort on face size to ensure ascending order.
		// Maybe this is overkill. But computationally it's quite cheap, so who cares?
		std::vector<Int> ascOrderChain;
		for (Int faceOfChain : chain) {
			ascOrderChain.push_back(faceOfChain);
		}
		std::sort(ascOrderChain.begin(), ascOrderChain.end(), lattice_face_idx_comparator(facesOfLattice));

		// build matrix for chain
		int rowNum = 0;
		for (Int faceOfChain : ascOrderChain) {
			std::vector<Int> setInChain = facesOfLattice[faceOfChain];
			for (Int elem : setInChain) {
				chainMatrix(rowNum, elem) = 1;
			}
			++rowNum;
		}
		listOfChainMatrices.push_back(chainMatrix);
	}

	return listOfChainMatrices;
}

UserFunctionTemplate4perl("# @category Other"
                  "# Generates an array of chain matrices corresponding to the maximal chains of a matroid."
                  "# @param Matroid M a matroid"
                  "# @return Array<Matrix<Int>> each maximal chain matrix of the matroid",
                  "chains_of_flats(Matroid)");

} }
