use application 'matroid';
use Data::Dump qw(dump);
#$Polymake::User::Verbose::cpp = 1;

# Obtains size of the ground set, rank, and the class of matroids to be conisdered.
sub getUserInfo {
	print "This script will compute the dimension of one of the (r,n)-homology groups. \n";
	print "How many elements in the matroids: ";
	my $n = <>;
	chomp($n);
	print "What is the rank of the matroids: ";
	my $r = <>;
	chomp($r);
	
	my $oper = 0;
	while($oper != 1 and $oper != 2){
		print "Contraction (1) or deletion (2)? ";
		$oper = <>;
		chomp($oper);
	}
	
	my $matroidType = -1;
	my $allowedTypes = new Set(0,1,2,3);
	
	while($allowedTypes->contains($matroidType) != 1){
		print "What class of matroids are considered? \n";
		print "ALL = NESTED = 0 , REGULAR = 1, BINARY = 2 , TERNARY = 3: \n";
		$matroidType = <>;
		chomp($matroidType);
	}
	return ($n,$r,$oper,$matroidType);
}

sub loadChains {
	my $r = $_[0];
	my $n = $_[1];
	my $path = "chainfiles/chainIndices_r${r}_n${n}.txt";
	print "Loading chains from $path...\n";
	my $fileExists = open (FH,'<',$path);
	my $chains;
	if ($fileExists) {
		my $contents = do {local $/; <FH> };
		$chains = eval($contents);
		close(FH);
	} else {
		print "\tNot found, generating in-memory. For faster execution, ensure you have run setup.sh and are providing a reasonable (r,n).\n";
		$chains = generate_chain_matrices($r,$n);
	}

	my $numChains = scalar(@$chains);
	print "\tGot $numChains chains!\n";
	return $chains;
}

#Open polyDB and create proper cursor
#Currently uses ifelifelse structure rather than passing a string. Should be changed when added to extension.
sub collectionFromPolyDB {
	my $n = $_[0];
	my $r = $_[1];
	my $matroidType = $_[2];
	print "Loading matroids from PolyDB...\n";
	polyDB;
	$PolyDB::default::db_section_name = "Matroids";
	$PolyDB::default::db_collection_name = "Small";
	my $collection = polyDB->get_collection();
	#For the dual case ask that it have no coloops
	my $matroidsRN;
	my $matroidsCount;
	
	#"ALL = NESTED = 0 , REGULAR = 1, BINARY = 2 , TERNARY = 3
	if($matroidType == 0){
		if($n/2>=$r){
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,"NESTED"=>true});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,"NESTED"=>true});
		}
		else{
			# We use the fact that matorid classes of a fixed cyclic width are dual-closed, in particular cyclic width of one.
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0,"NESTED"=>true});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0,"NESTED"=>true});
		}
	}
	
	if($matroidType == 1){
		if($n/2>=$r){
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,"REGULAR"=>true});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,"REGULAR"=>true});
		}
		else{
			# Regular matroids are dual-closed
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0,"REGULAR"=>true});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0,"REGULAR"=>true});
		}
	}
	
	if($matroidType == 2){
		if($n/2>=$r){
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,"BINARY"=>true});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,"BINARY"=>true});
		}
		else{
			# Binary matroids are dual-closed
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0,"BINARY"=>true});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0,"BINARY"=>true});
		}
	}
	if($matroidType == 3){
		if($n/2>=$r){
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,"TERNARY"=>true});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,"TERNARY"=>true});
		}
		else{
			# Ternary matroids are dual-closed
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0,"TERNARY"=>true});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0,"TERNARY"=>true});
		}
	}
	
	return ($matroidsRN,$matroidsCount);
}

#This computes the S_n representation given by the action on the indicator vectors.
# This is the slower of the two Perl implementations, in future this will be implemented in C++.
sub representationOfSn {
	my $n = $_[0];
	my $r = $_[1];
	my $allChains = $_[2];
	my $length = $_[3];
	print "Precomputing regular permutations on matroids...\n";
	my $listOfMatrixReps = [];
	my $group = group::symmetric_group($n);
	for my $perm (@{group::all_group_elements($group->REGULAR_REPRESENTATION)}){
		my $matrix = new Matrix($length,$length);
		for my $originalIndex (0..($length-1)){
			my $chainMatrix = new Matrix($allChains->[$originalIndex]);
			my $permChainMatrix = $chainMatrix*$perm;
			for my $finalIndex (0..($length-1)){
				if($permChainMatrix == new Matrix($allChains->[$finalIndex])){
					$matrix->elem($originalIndex,$finalIndex) = 1;
				}
			}
		}
		push @$listOfMatrixReps, $matrix;
	}
	return $listOfMatrixReps;
}

# What it says on the tin. Takes in the list of chain matrices corresponding to a matroid and computes an indicator vector to describe the matroid.
sub indicatorVectorProducer {
	my $allChains = $_[0];
	my $dim = scalar(@$allChains);
	my $listOfChainMatrices = $_[1];
	
	my $indicatorVector = zero_vector($dim);
	for my $index (0..($dim-1)){
		#print "$index\n";
		my $indexChain = new Matrix($allChains->[$index]);
		my $check =0;
		for my $chainMatrix (@{$listOfChainMatrices}){
			#print $indexChain;
			#print "\n = \n";
			#print $chainMatrix;
			#print "\n \n";
			#<>;
			if($indexChain == $chainMatrix){
				$check =1;
			}
		}
		if($check ==1){
			$indicatorVector = $indicatorVector + unit_vector($dim,$index);
		}
	}
	return $indicatorVector;
}

# This subroutine converts a matroid into a set of matrices, each matrix representing a chain of flats. It is inefficient and a better implementation is known and applied in C++ using polymake methods.
sub chainsOfFlats {
	my $matroid = $_[0];
	my $n = $matroid->N_ELEMENTS;
	my $r = $matroid->RANK;
	my $latticeOfFlats = $matroid->LATTICE_OF_FLATS;
	my $chains = maximal_chains_of_lattice($latticeOfFlats);
	my $listOfChainMatrices = [];
	
	for my $chain (@{$chains}){
		#print $chain;
		#print "\n";
		#<>;
		my $chainMatrix = zero_matrix(0,$n);
		my $chainMatrix = new Matrix($r+1,$n);
		for my $rowCounter (0..$r){
			my $faceOfChain = $chain->[-(1+$rowCounter)];
			my $setInChain = $latticeOfFlats->FACES->[$faceOfChain];
			for my $elem (@{$setInChain}){
				$chainMatrix->elem($rowCounter,$elem) = 1;
			}
		}
		#This line is needed due to the way the lattice of flats is transversed. Sometimes it will began at the emptyset and othertimes at the entire set. This check determines if it's transversed it backward and flips the matrix if so.
		if($chainMatrix->row(0) != zero_vector($n)){
			my $flipChainMatrix = new Matrix($r+1,$n);
			for my $rowCounter (0..$r){
				$flipChainMatrix->row($r-$rowCounter)=$chainMatrix->row($rowCounter);
			}
			$chainMatrix = $flipChainMatrix;
		}
		
		push @$listOfChainMatrices,$chainMatrix;
	}
	return $listOfChainMatrices;
}

sub imageKernelDeletion {
	my $matroid = $_[0];
	my $n = $matroid->N_ELEMENTS;
	my $rangeChains = $_[1];
	my $lengthOfRangeVector = scalar(@$rangeChains);
	
	my $imageOfDeletion = zero_vector($lengthOfRangeVector);
	
	for my $i (0..$n-1){
		if($matroid->COLOOPS->contains($i) == false){
			my $matroidMinor = deletion($matroid,$i);
			my $listOfChainMatrices = chainsOfFlats($matroidMinor);
			my $indicatorVector = indicatorVectorProducer($rangeChains, $listOfChainMatrices);
			$imageOfDeletion = $imageOfDeletion+((-1)**$i)*$indicatorVector;
		}
	}
	
	return $imageOfDeletion;
}

sub imageContractionDeletion {
	my $matroid = $_[0];
	my $n = $matroid->N_ELEMENTS;
	my $rangeChains = $_[1];
	my $lengthOfRangeVector = scalar(@$rangeChains);
	
	my $imageOfContraction = zero_vector($lengthOfRangeVector);
	
	for my $i (0..$n-1){
		#The condition here is that closure(i) does not equal itself.
		my $possibleFlat = new Set($i);
		my $containsI = false;
		for my $flat (@{$matroid->LATTICE_OF_FLATS->FACES}){
			if($possibleFlat == $flat){
				$containsI = 1;
			}
		}
		if($containsI){
			my $matroidMinor = contraction($matroid,$i);
			my $listOfChainMatrices = chainsOfFlats($matroidMinor);
			my $indicatorVector = indicatorVectorProducer($rangeChains, $listOfChainMatrices);
			$imageOfContraction = $imageOfContraction+((-1)**$i)*$indicatorVector;
		}
	}
	
	return $imageOfContraction;
}

sub main {
	my ($n,$r,$oper,$matroidType) = getUserInfo();
	my $nMinusOne = $n-1;
	my $rMinusOne = $r-1;
	my $dualCheck = 0;
		
	#Obtains cursor from polyDB
	my ($matroidsRN,$matroidsCount) = collectionFromPolyDB($n,$r,$matroidType);
	
	if($r>$n/2){
		$dualCheck =1;
	}
	
	# Obtains chainIndices from correct file. All chains produced by chain_generator_python.py
	# Maps are considered d:(M_{r,n}->M_{r,n-1}) or c:(M_{r,n}->M_{r,n-1}
	my $domainChains = loadChains($r,$n);

	# $oper = 1 is contraction and 2 is deletion.
	my $rangeChains;
	if($oper == 1){
		$rangeChains = loadChains($rMinusOne, $nMinusOne);
	}
	elsif($oper == 2){
		$rangeChains = loadChains($r, $nMinusOne);
	}
	
	my $lengthOfDomainVector = scalar(@$domainChains);
	my $lengthOfRangeVector = scalar(@$rangeChains);
	
	#We initialize the matrix of domain and range indictor vectors, including some redundancies due to relabeling.
	my $facN = new Int(fac($n));
	my $group = group::symmetric_group($n);
	my $listOfMatrixRepOfSn = representationOfSn($n,$r,$domainChains,$lengthOfDomainVector);
	my $matrixOfDomainIndicatorVectors = new Matrix($matroidsCount*$facN,$lengthOfDomainVector);
	my $matrixOfRangeIndicatorVectors = new Matrix($matroidsCount*$facN,$lengthOfRangeVector);
	my $indVecRowCount = 0;
	
	while($matroidsRN->PolyDB::Cursor::has_next){
		my $matroid = $matroidsRN->PolyDB::Cursor::next;
		if($dualCheck){
			$matroid = dual($matroid);
		}
		
		#The following computations relate only to the domain and are not affected by the operations, so we can later use permutation matrices to obtain all other indicator vectors just by having one in the isomorphism class
		my $listOfChainMatrices = chainsOfFlats($matroid);
		my $domainIndicatorVector = indicatorVectorProducer($domainChains, $listOfChainMatrices);
		#my $domainIndicatorVector = indicatorVectorProducer($matroid,$n,$r,$domainChains,$lengthOfDomainVector,$listOfChainMatrices);
		my $permCounter = 0;
		
		#Unfortunantly, the kernel is not generated by single matroid indicator vectors, but by some linear combinations of them (i.e. m1+m2 may be in the kernel but not m1 or m2). Because of this, we need to store the actual image rather than just check if it's zero. Moreover, relabelings change the signs of the sum on the minors. It is possible that there exists a matroid whose image is unrelated to the image of a relabeled matroid by permutations. This requires expensive computations at this time. If the action of S_n on deletion is known, this can be resolved. I.e. if one can guess \sigma \partial_d(m) when given \partial_d( \sigma m).
		for my $perm (@{group::all_group_elements($group->PERMUTATION_ACTION)}){
			my $permBases = [];
			for my $base (@{$matroid->BASES}){
				push @{$permBases}, group::action($perm,$base);
			}
			my $permMatroid = new Matroid(N_ELEMENTS => int($n),BASES=>$permBases);
			#Here we can use the nice properties of the domain to avoid calling indicatorVectorProducer.
			$matrixOfDomainIndicatorVectors->row($indVecRowCount) = $domainIndicatorVector*$listOfMatrixRepOfSn->[$permCounter];
			
			my $rangeIndicatorVector;
			if($oper == 1){
				$rangeIndicatorVector = imageContractionDeletion($permMatroid,$rangeChains);
			}
			elsif($oper == 2){
				$rangeIndicatorVector = imageKernelDeletion($permMatroid,$rangeChains);
			}
			$matrixOfRangeIndicatorVectors->row($indVecRowCount) = $rangeIndicatorVector;
			
			$permCounter = $permCounter+1;
			$indVecRowCount = $indVecRowCount+1;
		}
		
	}
	#I have it print out an extra row so that we can be sure the printed matrix contins all entries. If the final row is non-zero, something has gone wrong.
	
	#print $matrixOfDomainIndicatorVectors;#->minor([0..$indVecRowCount],All);
	#print "\n";
	#print $matrixOfRangeIndicatorVectors;
	
	print "Dim of domain M_{r,n}: ";
	print rank($matrixOfDomainIndicatorVectors);
	print "\n";
	print "Dim of image M_{r,n}->: ";
	print rank($matrixOfRangeIndicatorVectors);
	print "\n";
	
	print "Dim of kernel: M_{r,n}->:";
	print rank($matrixOfDomainIndicatorVectors)-rank($matrixOfRangeIndicatorVectors);
	print "\n";
	exit;
}

main()
