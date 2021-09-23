use application 'matroid';
use Data::Dump qw(dump);

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
	
	print "Class of matroids (use all caps and a polyDB recognized keyword e.g. NESTED, REGULAR, BINARY, TERNARY): ";
	my $matroidType = <>;
	chomp($matroidType);
	
#	my $matroidType = -1;
#	my $allowedTypes = new Set(0,1,2,3);
#	while($allowedTypes->contains($matroidType) != 1){
#		print "What class of matroids are considered? \n";
#		print "ALL = NESTED = 0 , REGULAR = 1, BINARY = 2 , TERNARY = 3: \n";
#		$matroidType = <>;
#		chomp($matroidType);
#	}
	return ($n,$r,$oper,$matroidType);
}

sub loadChainMatrices {
	my $r = $_[0];
	my $n = $_[1];
	my $path = "chainfiles/chainIndices_r${r}_n${n}.txt";
	print "Loading chain matrices from $path...\n";
	my $fileExists = open (FH,'<',$path);
	my $chainMatrices = [];
	if ($fileExists) {
		# load as 2D arrays and cast to Matrix for easier handling
		my $contents = do {local $/; <FH> };
		my $chains = eval($contents);
		for my $chain (@$chains) {
			my $matrix = new Matrix<Int>($chain);
			push @$chainMatrices, $matrix;
		}
		close(FH);
	} else {
		print "\tNot found, generating in-memory. For faster execution, ensure you have run setup.sh and are providing a reasonable (r,n).\n";
		$chainMatrices = generate_chain_matrices($r,$n);
	}

	my $numChains = scalar(@$chainMatrices);
	print "\tGot $numChains chain matrices!\n";
	return new Array<Matrix<Int>>($chainMatrices);
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
	
	if($n/2>=$r){
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,$matroidType=>true,NESTED=>true,REGULAR=>true});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,$matroidType=>true,NESTED=>true,REGULAR=>true});
		}
		else{
			# We use the fact that matorid classes of a fixed cyclic width are dual-closed, in particular cyclic width of one.
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0,$matroidType=>true,NESTED=>true,REGULAR=>true});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0, $matroidType=>true,NESTED=>true,REGULAR=>true});
		}

	if($matroidType eq "UNIFORM"){
		my $facN = new Int(fac($n));
		if($n/2>=$r){
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,"N_AUTOMORPHISMS" => int($facN)});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,"N_AUTOMORPHISMS" => int($facN)});
		}
		else{
			# We use the fact that matorid classes of a fixed cyclic width are dual-closed, in particular cyclic width of one.
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0, "N_AUTOMORPHISMS" => int($facN)});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0, "N_AUTOMORPHISMS" => int($facN)});
		}
	}
	
	if($matroidsCount == 0){
		die("Class of matroids is empty or not polyDB supported. Aborting program. \n");
	}

	return ($matroidsRN,$matroidsCount);
}

#This computes the S_n representation given by the action on the indicator vectors.
# This is the faster of the two Perl implementations, in future this will be implemented in C++.
sub representationOfSn {
	my $n = $_[0];
	my $r = $_[1];
	my $allChainMatrices = $_[2];
	my $length = $_[3];

	print "Precomputing regular permutations on matroids...\n";
	my $listOfMatrixReps = [];
	my $group = group::symmetric_group($n);
	# Pre-hash the indexes of each chain to make reverse lookups after permutation faster
	my %chainIndexes = ();
	for my $index (0..($length-1)){
		my $indexChain = $allChainMatrices->[$index];
		$chainIndexes{"${indexChain}"} = $index;
	}
	# Consider each permutation on the elements within each chain (n x n) 
	for my $perm (@{group::all_group_elements($group->REGULAR_REPRESENTATION)}){
		my $matrix = new Matrix($length,$length);
		for my $originalIndex (0..($length-1)){
			# Permute this chain matrix and lookup its index to determine the row in the indicator permutation matrix
			my $chainMatrix = $allChainMatrices->[$originalIndex];
			my $permChainMatrix = $chainMatrix*$perm;
			my $permChainIndex = $chainIndexes{"${permChainMatrix}"};
			$matrix->elem($originalIndex,$permChainIndex) = 1;
		}
		push @$listOfMatrixReps, $matrix;
	}
	return $listOfMatrixReps;
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
	
	# Operation 1 is contraction and 2 is deletion
	my $domainChainMatrices = loadChainMatrices($r,$n);
	my $rangeChainMatrices;
	if($oper ==1){
		$rangeChainMatrices = loadChainMatrices($r-1,$n-1);
	}
	else{
		$rangeChainMatrices = loadChainMatrices($r,$n-1); 
	}
	
	my $lengthOfDomainVector = scalar(@$domainChainMatrices);
	my $lengthOfRangeVector = scalar(@$rangeChainMatrices);
	
	#We initialize the matrix of domain and range indictor vectors, including some redundancies due to relabeling.
	my $facN = new Int(fac($n));
	my $group = group::symmetric_group($n);
	my $listOfMatrixRepOfSn = representationOfSn($n,$r,$domainChainMatrices,$lengthOfDomainVector);
	my $matrixOfDomainIndicatorVectors = new Matrix($matroidsCount*$facN,$lengthOfDomainVector);
	my $matrixOfRangeIndicatorVectors = new Matrix($matroidsCount*$facN,$lengthOfRangeVector);
	my $indVecRowCount = 0;
	
	while($matroidsRN->PolyDB::Cursor::has_next){
		print "Examining matroid...\n";
		my $matroid = $matroidsRN->PolyDB::Cursor::next;
		if($dualCheck){
			$matroid = dual($matroid);
		}
		
		#The following computations relate only to the domain and are not affected by the operations, so we can later use permutation matrices to obtain all other indicator vectors just by having one in the isomorphism class
		my $listOfChainMatrices = chains_of_flats($matroid);
		my $domainIndicatorVector = indicator_vector_full($matroid, $domainChainMatrices, $listOfChainMatrices);
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
				$rangeIndicatorVector = contraction_vector_full($permMatroid,$rangeChainMatrices);
			}
			elsif($oper == 2){
				$rangeIndicatorVector = deletion_vector_full($permMatroid,$rangeChainMatrices);
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
