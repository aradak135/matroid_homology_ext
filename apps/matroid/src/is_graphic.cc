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
	
//Forward declaration.
bool binary_graphic(BigObject m);

//This even partition problem is equivalent to finding a 2-coloring on the intersection graph of bridgesY.

bool isYeven (BigObject m, Set<Int> deleteCols, Vector<Set<Int>> bridgesY, Set<Int> groundSet){
	//0 unallocated, 1 allocated
	Vector<Int> colors (bridgesY.size(),0);
	Int numColored = 0;
	
	while(numColored < bridgesY.size()){
		Int nextUncolored;
		for(Int i=0; i<bridgesY.size();++i){
			if(colors[i] == 0){
				nextUncolored = i;
				break;
			}
		}
		Set<Int> currentBridge = bridgesY[nextUncolored];
		colors[nextUncolored] = 1;
		numColored +=1;
		Set<Int> adjSetIndex;
		adjSetIndex += nextUncolored;
		
		while(adjSetIndex.size()>0){
			Int adjIndex = adjSetIndex.front();
			adjSetIndex-= adjIndex;
			for(Int j =0; j< bridgesY.size();++j){
				Set<Int> tempInter = bridgesY[adjIndex]*bridgesY[j];
				if(tempInter.size()>0 and adjIndex != j){
					if(colors[adjIndex] == colors[j]){
						return false;
					}
					else if(colors[j] == 0){
						colors[j] = -1*colors[adjIndex];
						numColored+=1;
						adjSetIndex += j;
					}
				}
			}
		}
	}
	return true;
}
	
bool bridges(BigObject m, Matrix<Int> repMatrix, Int colIndex, Int rowIndex) {
	Int colCount = repMatrix.col(0).size();
	Int rowCount = repMatrix.row(0).size();
	
	Set<Int> deleteCols;
	
	for(Int j=0; j<colCount ;++j){
		if(repMatrix(rowIndex,j) ==1){
			deleteCols += j;
		}
	}
	
	Int n = m.give("N_ELEMENTS");
	Set<Int> groundSet;
	Vector<Int> groundLabels;
	Vector<Int> deletionLabels;
	
	for (Int elem = 0; elem<n;++elem) {
		groundSet += elem;
		groundLabels = groundLabels| elem;
		bool containsBool = call_function("contains",deleteCols,elem);
		if(containsBool == false){
			deletionLabels = deletionLabels | elem;
		}
	}
	
	BigObject mY = call_function("deletion",m,deleteCols);
	Int numComMY = mY.give("N_CONNECTED_COMPONENTS");
	
	if(numComMY==1){
		return false;
	}
	Vector<Set<Int>> componentsY = mY.give("CONNECTED_COMPONENTS");
	Vector<Set<Int>> bridgesY;
	
	for(Int comp;comp<componentsY.size();++comp){
		Set<Int> bridgeY;
		Set<Int> currentComp = componentsY[comp];
		//Vector<Int> currentComp = currentCompSet;
		for (Int index =0; index <deletionLabels.size(); ++index){
			if(call_function("contains",currentComp,index)){
				bridgeY += deletionLabels[index];
			}
		}
		bridgesY = bridgesY | bridgeY;
	}
	
	if(isYeven(m,deleteCols,bridgesY,groundSet) == false){
		return false;
	}
	
	for (Int bridgeNum=0; bridgeNum < bridgesY.size();++bridgeNum){
		Set<Int> bridgeAndY = deleteCols+bridgesY[bridgeNum];
		Set<Int> complementbridgeY = groundSet-bridgeAndY;
		BigObject bridgeMatroid = call_function("contraction",m,complementbridgeY);
		if(binary_graphic(bridgeMatroid) == false){
			return false;
		}
	}
	return true;
}

bool binary_graphic(BigObject m) {
	Matrix<Int> repMatrixTrans = m.give("BINARY_VECTORS");
	Matrix<Int> repMatrix = call_function("transpose",repMatrixTrans);
	
	Int colCount = repMatrix.col(0).size();
	Int rowCount = repMatrix.row(0).size();
	
	bool mod2Check = true;
	Int colIndex;
	
	for (Int colNum = 0;colNum < colCount; ++colNum) {
		Vector<Int> currentColumn = repMatrix.col(colNum);
		Int colSum = 0;
		for(Int i =0; i<rowCount;++i){
			colSum += currentColumn[i];
		}
		if(colSum >2){
			mod2Check = false;
			colIndex = colNum;
			break;
		}
	}
	
	if(mod2Check){
		return true;
	}
	
	//The following implementation is inefficient, but was written like this to resolve a few issues.
	Set<Int> rowIndices;
	for (Int i=0;i<rowCount;++i){
		if(repMatrix(i,colIndex)==1){
			rowIndices += i;
		}
	}
	Int componentCounter = 0;
	//Only need to check three, editted out for debugging.
	while(rowIndices.size()>0){
		Int rowIndex = rowIndices.front();
		rowIndices -= rowIndex;
		if(bridges(m,repMatrix,colIndex,rowIndex)){
			return true;
		}
	}
	return false;
}
	
bool connected_component_breakdown(BigObject m) {
	
	if(m.give("CONNECTED")){
		return binary_graphic(m);
	}
	Int n = m.give("N_ELEMENTS");
	Set<Int> groundSet;
	for (Int elem = 0; elem<n;++elem) {
		groundSet += elem;
	}
	Vector<Set<Int>> components = m.give("CONNECTED_COMPONENT");
	
	for (Int componentNumber =0 ; componentNumber < m.give("N_CONNECTED_COMPONENTS"); ++componentNumber) {
		Set<Int> component = components[componentNumber];
		BigObject restrictionMatroid = call_function("deletion",m,groundSet-component);
		if(binary_graphic(restrictionMatroid)==false){
			return false;
		}
	}
	return true;
}

bool is_graphic(BigObject m) {
	if (m.give("REGULAR")){
		return connected_component_breakdown(m);
	}
	else {
		return false;
	}
}

UserFunctionTemplate4perl("# @category Other"
                  "# Check if a matroid is a graphic matroid using Tutte's algorithm"
                  "# @param Matroid a matroid to be checked"
                  "# @return bool",
                  "is_graphic(Matroid)");
} }