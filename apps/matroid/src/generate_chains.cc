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

namespace polymake { namespace matroid {

// Enable this code to perform print-debugging.
void print_chain(std::vector<std::vector<Int>> chain) {
	int first = 1;
	for (std::vector<Int> link : chain) {
		if (!first) {		
			printf(" -> ");
		}
		for (Int elem : link) {
			printf("%ld", elem);
		}
		first = 0;
	}		
	printf("\n");
}

// recursive method to build up a chain one link at a time to a fixed length r
void recursive_build_chains(std::vector<std::vector<std::vector<Int>>> & global_chain_list, std::vector<std::vector<Int>> chain_so_far, int remaining_positions, int r, int n) {
	//print_chain(chain_so_far);
	std::vector<Int> current_chain_value(n,0);
	if (chain_so_far.size()) {
		current_chain_value = chain_so_far.back();
	}
	
	// base success case: exactly one jump from final value
	if (chain_so_far.size() == r && remaining_positions > 0) {
		// TODO: appending all chains with the 1 vector should not be necessary
		std::vector<Int> all_ones(n, 1);
		chain_so_far.push_back(all_ones);
		global_chain_list.push_back(chain_so_far);
		return;
	}

	// base failure case: no more jumps (or, too many for remaining positions)
	if (chain_so_far.size() >= r || r - chain_so_far.size() >= remaining_positions) {
		return;
	}

	// recursion case: attempt each additional jump which flips one or more elements to 1
	// do so via binary incrementat over the remaining slots (e.g. 001, 010, 011, 100, 101, 110, 111)
	int target_inc = 2 << remaining_positions - 1;
	for (int inc = 1; inc < target_inc; ++inc) {
		std::vector<Int> next_chain_value(current_chain_value);
		int zero_elems_checked = 0;
		int zero_elems_set = 0;
		for (int i = 0; i < n; ++i) {
			// only examine 0 values left in the chain
			if (next_chain_value[i] == 0) {
				// extract the increment for this 0 value to determine whether to flip it to 1
				if (inc & (1 << (remaining_positions - zero_elems_checked - 1))) { 
					next_chain_value[i] = 1;
					++zero_elems_set;
				}
				++zero_elems_checked;
			}
		}

		// add to the chain_so_far and recurse, if sufficient slots remaining
		if (remaining_positions - zero_elems_set >= r - chain_so_far.size() - 1) {
			chain_so_far.push_back(next_chain_value);
			recursive_build_chains(global_chain_list, chain_so_far, remaining_positions - zero_elems_set, r, n);
			chain_so_far.pop_back(); // reset for next
		}
	}
}

std::vector<std::vector<std::vector<Int>>> generate_chains(Int r, Int n) {
	// Initialize data structure to store all chains & initial chain
	std::vector<std::vector<std::vector<Int>>> all_chains;
	std::vector<std::vector<Int>> base_chain;
	// TODO: prepending all chains with the 0 vector should not be necessary
	std::vector<Int> base_chain_starting_value(n,0);
	base_chain.push_back(base_chain_starting_value);

	// Enter recursive call to generate all possible paths from [0,0,...0] to [1,1,...1] of lengh r
	recursive_build_chains(all_chains, base_chain, (int)n, (int)r, (int)n);

	return all_chains;
}


// Variant of generate_chains which returns chains as matrices rather than a 2D array.
// This is possibly a bit slower because of abstraction, but makes permutation and comparison easier.
std::vector<Matrix<Int>> generate_chain_matrices(Int r, Int n) {
	std::vector<std::vector<std::vector<Int>>> all_chains = generate_chains(r, n);
	std::vector<Matrix<Int>> all_chain_matrices;
	for (std::vector<std::vector<Int>> chain : all_chains) {
		Matrix<Int> chainMatrix(r+1, n);
		for (int rowNum = 0; rowNum < r+1; ++rowNum)
			for (int colNum = 0; colNum < n; ++colNum)
				chainMatrix(rowNum, colNum) = chain[rowNum][colNum];
		all_chain_matrices.push_back(chainMatrix);
	}
	return all_chain_matrices;
}

UserFunctionTemplate4perl("# @category Other"
                  "# Generates an array of chains representing a depth first search across the power set 2^n."
				  "# @param Int R the rank (# of links) of each chain"
                  "# @param Int N the size (# of elements) of each link in each chain"
                  "# @return Array<Array<Array<Int>>> An array of chains each representing a distinct path from [0,...0] to [1,...1]",
                  "generate_chains(Int, Int)");

UserFunctionTemplate4perl("# @category Other"
                  "# Generates an array of chain matrices representing a depth first search across the power set 2^n."
                  "# @param Int R the rank (# of links) of each chain"
                  "# @param Int N the size (# of elements) of each link in each chain"
                  "# @return Array<Matrix<Int>> An array of chain matrices each representing a distinct path from [0,...0] to [1,...1]",
                  "generate_chain_matrices(Int, Int)");
} }
