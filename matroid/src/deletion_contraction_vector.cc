#include "polymake/client.h"
#include "polymake/Vector.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/graph/Lattice.h"
#include "polymake/graph/Decoration.h"
#include "polymake/internal/sparse2d.h"
#include "polymake/GenericSet.h"
#include <list>
#include <algorithm>
#include <stdio.h>
#include<iostream>
#include<bitset>
#include <math.h>
#include <vector>
#include <typeinfo>

namespace polymake { namespace matroid {
	
Int chain_length (Int r, Int n) {
	// Hard code length of chain here for all reasonable n,r pairs.
	std::vector<Matrix<Int>> allChains = call_function("generate_chain_matrices",r,n);
	Int lenOfVec = allChains.size();
	return lenOfVec;
}

Vector<Int> deletion_vector(BigObject m) {
	
	Int n = m.give("N_ELEMENTS");
	Int r = m.give("RANK");
	
	Int lenOfVec = chain_length(r,n-1);
	Vector<Int> delVector (lenOfVec,0);
	Vector<Int> tempVec (lenOfVec,0);

	// TODO: There needs to be a way to obtain COLOOPS directly, but give("COLOOPS") isn't recognized by polymake. Taking the dual is slow, but works.
	BigObject dualM = call_function("dual",m);
	Vector<Int> coloops = dualM.give("LOOPS");
	BigObject currentDelM;

	Int coloopFlag = 0;

	for(Int elem; elem<n; ++elem){
		coloopFlag =0;
		for(Int coloop: coloops){
			if(coloop == elem){
				coloopFlag = 1;
			}
		}
		if (coloopFlag == 0){
			currentDelM = call_function("deletion",m,elem);
			tempVec = call_function("indicator_vector",currentDelM);
			delVector = delVector+pow(-1,elem)*tempVec;
			//delVector = delVector + call_function("indicator_vector",currentDelM);
		}
	}
/*
	for(Int elem;elem<n;++elem){
		Int coloopFlag = 0;
		for (Int coloop : coloops){
			if (coloop == elem){
				coloopFlag =1;
				printf("%d",elem);
				printf("\n");
			}
		}
	}
*/
	return delVector;
}	

Vector<Int> contraction_vector(BigObject m) {
	Int n = m.give("N_ELEMENTS");
	Int r = m.give("RANK");
	
	Int lenOfVec = chain_length(r-1,n-1);
	Vector<Int> conVector (lenOfVec,0);
	Vector<Int> tempVec (lenOfVec,0);

	Array<Set<Int>> nodesOfFlatsRankOne = m.give("MATROID_HYPERPLANES");

	BigObject currentConM;
	Int flatFlag = 0;

	for (Int elem; elem<n;++elem){
		flatFlag = 0;
		for(Set<Int> flat : nodesOfFlatsRankOne){
			if(call_function("contains",flat,elem)){
				flatFlag = 1;
			}
		}
		if (flatFlag ==1){
			currentConM = call_function("contraction",m,elem);
			tempVec = call_function("indicator_vector",currentConM);
			conVector = conVector+ pow(-1,elem)*tempVec;
		}
	}
	return conVector;
}


UserFunctionTemplate4perl("# @category Other"
                  "# Check if a given list of sets satisfies the axioms to be the flats of a matroid."
                  "# @param Matroid a matroid whose indicator vector is to be computed"
                  "# @return std::vector<Int>",
                  "contraction_vector(Matroid)");

UserFunctionTemplate4perl("# @category Other"
                  "# Check if a given list of sets satisfies the axioms to be the flats of a matroid."
                  "# @param Matroid a matroid whose indicator vector is to be computed"
                  "# @return Vector<Int>",
                  "deletion_vector(Matroid)");
} }