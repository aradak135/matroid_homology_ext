#include "polymake/client.h"
#include "polymake/Vector.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"

#include <list>
#include <algorithm>
#include <stdio.h>
#include<iostream>
#include<bitset>
#include <math.h>
#include <vector>

namespace polymake { namespace matroid {

BigObject matroid_product(BigObject m, BigObject n) {
	
	BigObject intersectionProduct = call_function("intersection",m,n);
	Int nLoops = intersectionProduct.give("N_LOOPS");
	if (nLoops > 0) {
		Int nElem = m.give("N_ELEMENTS");
		intersectionProduct = call_function("uniform_matroid",0,nElem);
		cout << "Product is zero. \n";
	}

	return intersectionProduct;
}

Vector<Int> matroid_addition(Array<BigObject> arrayMatroids, Vector<Int> coefficients) {
	Int numCoeff = coefficients.size();
	
	BigObject currentMatroid = arrayMatroids[0];
	Vector<Int> currentVector = call_function("indicator_vector",currentMatroid);
	Int currentCoeff = coefficients[0];
	Vector<Int> indicatorVectorSum = currentCoeff*currentVector;
	
	for (int count =1; count < numCoeff; count++) {
		currentMatroid = arrayMatroids[count];
		
		if ( currentMatroid.give("N_LOOPS") == 0){
			currentVector = call_function("indicator_vector",currentMatroid);
			currentCoeff = coefficients[count];
			indicatorVectorSum = indicatorVectorSum+currentCoeff*currentVector;
		}
	}
	
	return indicatorVectorSum;
}

UserFunctionTemplate4perl("# @category Other"
                  "# Generates a new matroid which is the intersection product of two matroids."
                  "# @param Matroid M a matroid"
				  "# @param Matroid N a matroid"
                  "# @return Matroid intersection product",
                  "matroid_product(Matroid,Matroid)");

UserFunctionTemplate4perl("# @category Other"
                  "# Takes in a linear combination of matroids and outputs sum of indicator vectors."
                  "# @param Array<Matroid> matroids"
				  "# @param Vector<Int> coefficients"
                  "# @return Vector<Int> sum of indicator vectors.",
                  "matroid_addition(Array<Matroid>,Vector<Int>)");
}
}
