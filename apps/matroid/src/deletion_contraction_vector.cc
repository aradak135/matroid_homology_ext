#include "polymake/client.h"
#include "polymake/Vector.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/GenericSet.h"
#include "polymake/graph/Lattice.h"
#include "polymake/graph/Decoration.h"

#include <list>
#include <algorithm>
#include <stdio.h>
#include<iostream>
#include<bitset>
#include <math.h>
#include <vector>

namespace polymake { namespace matroid {
	
Int chain_length (Int r, Int n) {
	Array<Matrix<Int>> allChains = call_function("generate_chain_matrices",r,n);
	Int lenOfVec = allChains.size();
	return lenOfVec;
}

Vector<Int> deletion_vector(BigObject m) {
	
	Int n = m.give("N_ELEMENTS");
	Int r = m.give("RANK");
	
	Int lenOfVec = chain_length(r,n-1);
	Vector<Int> delVector (lenOfVec,0);
	Vector<Int> tempVec (lenOfVec,0);

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
		}
	}
	return delVector;
}

Vector<Int> deletion_vector_full(BigObject m, Array<Matrix<Int>> rangeChains) {
	Int n = m.give("N_ELEMENTS");
	Int lenOfVec = rangeChains.size();
	Vector<Int> delVector (lenOfVec,0);
	Vector<Int> tempVec (lenOfVec,0);

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
			Array<Matrix<Int>> matroidDelChains = call_function("chains_of_flats",currentDelM);
			tempVec = call_function("indicator_vector_full",currentDelM,rangeChains,matroidDelChains);
			delVector = delVector+pow(-1,elem)*tempVec;
		}
	}
	return delVector;
}

Vector<Int> contraction_vector(BigObject m) {
	using polymake::graph::Lattice;
	using polymake::graph::lattice::BasicDecoration;
	using polymake::graph::lattice::Sequential;
	
	Int n = m.give("N_ELEMENTS");
	Int r = m.give("RANK");
	
	Int lenOfVec = chain_length(r-1,n-1);
	Vector<Int> conVector (lenOfVec,0);
	Vector<Int> tempVec (lenOfVec,0);
	
	Lattice<BasicDecoration, Sequential> latticeOfFlats = m.give("LATTICE_OF_FLATS");
	Set<Int> nodesOfFlatsRankOne = latticeOfFlats.nodes_of_rank(1);
	BigObject bigFaces = m.give("LATTICE_OF_FLATS");
	Array<Set<Int>> faces = bigFaces.give("FACES");
	
	BigObject currentConM;
	Int flag = 0;
	
	for (Int elem; elem <n; ++elem) {
		flag = 0;
		for (Int node : nodesOfFlatsRankOne ) {
			Set<Int> face = faces[node];
			if (call_function("contains",face,elem) and face.size() == 1) {
				flag = 1;
				break;
			}
		}
		
		if (flag == 1) {
			currentConM = call_function("contraction",m,elem);
			tempVec = call_function("indicator_vector",currentConM);
			conVector = conVector+ pow(-1,elem)*tempVec;
		}
	}
	return conVector;
}

Vector<Int> contraction_vector_full(BigObject m, Array<Matrix<Int>> rangeChains) {
	using polymake::graph::Lattice;
	using polymake::graph::lattice::BasicDecoration;
	using polymake::graph::lattice::Sequential;
	
	Int n = m.give("N_ELEMENTS");
	Int lenOfVec = rangeChains.size();
	Vector<Int> conVector (lenOfVec,0);
	Vector<Int> tempVec (lenOfVec,0);

	Lattice<BasicDecoration, Sequential> latticeOfFlats = m.give("LATTICE_OF_FLATS");
	Set<Int> nodesOfFlatsRankOne = latticeOfFlats.nodes_of_rank(1);
	BigObject bigFaces = m.give("LATTICE_OF_FLATS");
	Array<Set<Int>> faces = bigFaces.give("FACES");
	
	BigObject currentConM;
	Int flag = 0;
	
	for (Int elem; elem <n; ++elem) {
		flag = 0;
		for (Int node : nodesOfFlatsRankOne ) {
			Set<Int> face = faces[node];
			if (call_function("contains",face,elem) and face.size() == 1) {
				flag = 1;
				break;
			}
		}
		
		if (flag ==1){
			currentConM = call_function("contraction",m,elem);
			Array<Matrix<Int>> matroidConChains = call_function("chains_of_flats",currentConM);
			tempVec = call_function("indicator_vector_full", currentConM, rangeChains, matroidConChains );
			conVector = conVector+ pow(-1,elem)*tempVec;
		}
	}
	return conVector;
}

UserFunctionTemplate4perl("# @category Other"
                  "# Produces a vector corresponding to the alternating sum of the indicator vectorf of each single-element contraction"
                  "# @param Matroid a matroid whose indicator vector is to be computed"
                  "# @return Vector<Int>",
                  "contraction_vector(Matroid)");

UserFunctionTemplate4perl("# @category Other"
                  "# Produces a vector corresponding to the alternating sum of the indicator vectors of each single-element deletion"
                  "# @param Matroid a matroid whose indicator vector is to be computed"
                  "# @return Vector<Int>",
                  "deletion_vector(Matroid)");

FunctionTemplate4perl("contraction_vector_full(Matroid, *)");


FunctionTemplate4perl("deletion_vector_full(Matroid, *)");
} }
