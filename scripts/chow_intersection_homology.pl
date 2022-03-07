#From a matroid, generates the subring of the intersection ring which arises from
# the degree c component of the Chow ring of a matroid.
#We generate a monomial basis of the degree c component and then compute homology.
#Implemented by Austin Alderete

use application 'matroid';
use Data::Dump qw(dump);

#Define the collection to be pulled from polyDB.
sub collectionFromPolyDB {
	my $n = $_[0];
	my $r = $_[1];

	#print "Loading matroids from PolyDB...\n";
	polyDB;
	$PolyDB::default::db_section_name = "Matroids";
	$PolyDB::default::db_collection_name = "Small";
	my $collection = polyDB->get_collection();
	my $loopColoop = "N_LOOPS";
	
	if($n/2<$r){
		$r = $n-$r;
		$loopColoop = "DUAL.N_LOOPS";
	}
	
	#-----------------------------USER EDIT HERE-------------------------------------
	# Matroids must form a basis. Otherwise, use matroid_homology.
	#Class must be minor-closed.
	
	#Add attributes as desired here. Note that uniform matroids are obtained by setting N_AUTOMORPHISMS to n!.
	my $matroidsRN = $collection->find({
		N_ELEMENTS => int($n),
		RANK => int($r),
		$loopColoop => 0
		,NESTED=>true
		});
		
	#------------------------------END USER EDIT--------------------------------
	
	if($matroidsRN->PolyDB::Cursor::has_next){
		return $matroidsRN;
	}
	else{	
	die("Class of matroids is empty or not polyDB supported. Aborting program. \n");
	}
}

#This generates non-equal yet isomorphic matroids. It does so by computing all bases which arise from permutation and throwing out (via collect) duplicates.
sub permNaiveAllMat {
	my $matroidsRN = $_[0];
	my $n = $_[1];
	my $r = $_[2];
	my $group = group::symmetric_group($n);
	
	my $arrayAllMatroids = new Array<Matroid>(0);
	while($matroidsRN->PolyDB::Cursor::has_next){
		my $matroid = $matroidsRN->PolyDB::Cursor::next;
		
		if(2*$r>$n){
			$matroid = dual($matroid);
		}
		
		#------------------------USER EDIT-------------------------------
		#Additional properties are checked here.
		#if $matroid->N_CONNECTED_COMPONENTS-1 or $matroid->N_CONNECTED_COMPONENTS or $matroid->N_CONNECTED_COMPONENTS+1
		# if(-1 > $matroid->N_CONNECTED_COMPONENTS){
			# if($matroidsRN->PolyDB::Cursor::has_next){
				# next;
			# }
			# else{
				# last;
			# }
		# }
		#------------------END USER EDIT-----------------------

		my $allPossibleBases = new Set<Set<Set<Int>>>();
		
		for my $perm(@{group::all_group_elements($group->PERMUTATION_ACTION)}){
			my $permBases = [];
			for my $base (@{$matroid->BASES}){
				push @{$permBases}, group::action($perm,$base);
			}
			my $permBaseSet = new Set<Set<Int>>($permBases);
			$allPossibleBases->collect($permBaseSet);
		}
		
		for my $basesSet (@$allPossibleBases){
			my $permMatroid = new Matroid(BASES=>$basesSet,N_ELEMENTS=>$n,RANK=>$r);
			push @{$arrayAllMatroids}, $permMatroid;
		}
	}
	return $arrayAllMatroids;		
}

# Treats the given matroids as a basis and computes the deletion or contraction induced boundary map as a matrix with respect to this basis.
sub boundaryMap {
	my $allDomainMatroids = $_[0];
	my $allRangeMatroids = $_[1];
	my $n = $_[2];
	my $r = $_[3];
	my $domainDim = $_[4];
	my $rangeDim = $_[5];
	my $conOrdel = $_[6];
	my $boundaryMatrix = new Matrix<Int>($domainDim,$rangeDim);
	
	if ($conOrdel == 0){
		for (my $domCount = 0; $domCount <$domainDim;$domCount++){
			my $currentDomMatroid = $allDomainMatroids->[$domCount];
			for (my $elem = 0; $elem <$n; $elem++){
				my $conMat = contraction($currentDomMatroid,$elem);
				my $conMatBases = $conMat->BASES;
				
				for (my $ranCount = 0; $ranCount < $rangeDim ; $ranCount++){
					my $ranMatBases = $allRangeMatroids->[$ranCount]->BASES;
					if ($conMatBases == $ranMatBases ){
						$boundaryMatrix->elem($domCount,$ranCount) += (-1)**$elem;
					}
				}
			}
		}
	}
	
	elsif ($conOrdel ==1){
		for (my $domCount = 0; $domCount <$domainDim;$domCount++){
			my $currentDomMatroid = $allDomainMatroids->[$domCount];
			for (my $elem = 0; $elem <$n; $elem++){
				my $delMat = deletion($currentDomMatroid,$elem);
				my $delMatBases = $delMat->BASES;
				
				for (my $ranCount = 0; $ranCount < $rangeDim ; $ranCount++){
					my $ranMatBases = $allRangeMatroids->[$ranCount]->BASES;
					if ($delMatBases == $ranMatBases ){
						$boundaryMatrix->elem($domCount,$ranCount) += (-1)**$elem;
					}
				}
			}
		}
	}
	
	return $boundaryMatrix;
}

sub cleanUp {
	#This throws out duplicates as there are many.
	my $listOfMatroids = $_[0];
	
	my $cleanedUpList = new Array<Matroid>(0);
	my $setOfAllBases = new Set<Set<Set<Int>>>();
	for my $matroid (@{$listOfMatroids}){
		my $bases = new Set<Set<Int>>($matroid->BASES);
		if($setOfAllBases->contains($bases)){
			next;
		}
		else {
			push @{$cleanedUpList}, $matroid;
			$setOfAllBases->collect($bases);
		}
	}
	return $cleanedUpList;
}

sub principalChainTruncation  {
	my $m = $_[0];
	my $chainOfFlats = $_[1];
	my $trunMatroid = $m;
	
	for my $i (1..$chainOfFlats->dim){
		my $trunFlat = $chainOfFlats->[-$i];
		$trunMatroid = principal_truncation($trunMatroid,$trunFlat);
	}
	return $trunMatroid;
}

sub truncatorRecursion {
	my $m = $_[0];
	my $n = $_[1];
	my $r = $_[2];
	my $c = $_[3];
	my $lattice = $_[4];
	my $flat = $_[5];
	my $arrayMatroids = $_[6];
	my $nodes = $_[7];
	my $currentChainOfFlats = $_[8];
	my $currentRank = $_[9];
	
	if(scalar(@{$nodes}) == 0){
		return $arrayMatroids;
	}
	
	for my $node (@{$nodes}) {
		my $superFlat = $lattice->FACES->[$node];
		#This check requires that the rank(superSet)>rank(flat) otherwise we may have equality.
		if($superFlat+$flat != $superFlat){
			next;
		}
		my $rankSuper = $m->rank($superFlat);
		my $rankDiff = $rankSuper - $currentRank;
		for my $i (1..$rankDiff){
			my $localChainOfFlats = $currentChainOfFlats;
			for my $numberToAdd (1..$i){
				$localChainOfFlats = $localChainOfFlats | $superFlat;
			}
			if($localChainOfFlats->dim == $c){
				my $truncation = principalChainTruncation($m,$localChainOfFlats);
				#print $localChainOfFlats , "\n";
				if ($truncation->N_LOOPS == 0){
					#print $localChainOfFlats, "\n";
					#print $localChainOfFlats->dim, "\n";
					push @{$arrayMatroids}, $truncation;
				}
				#return $arrayMatroids;
			}
			elsif($localChainOfFlats->dim > $c){
				next;
			}
			else {
				$arrayMatroids = truncatorRecursion($m,$n,$r,$c,$lattice,$superFlat,$arrayMatroids,$lattice->nodes_of_rank_range($rankSuper,$r),$localChainOfFlats,$rankSuper)
			}
		}
	}
	return $arrayMatroids;
}

# This functions as main for deletion.
sub mainDel {
	my $relativeNestedBases = $_[0];
	my $testMatroid = $relativeNestedBases->[0];
	my $n = $testMatroid->N_ELEMENTS;
	my $r = $testMatroid->RANK;
	
	my $nmin1rMatroids = collectionFromPolyDB($n-1,$r,"NESTED");
	my $nmin2rMatroids = collectionFromPolyDB($n-2,$r,"NESTED");
	
	my $allnmin1rMatroids = permNaiveAllMat($nmin1rMatroids,$n-1,$r);
	my $allnmin2rMatroids = permNaiveAllMat($nmin2rMatroids,$n-2,$r);
	
	my $boundaryMatrix1 = boundaryMap($relativeNestedBases,$allnmin1rMatroids,$n,$r,$relativeNestedBases->size,$allnmin1rMatroids->size,1);
	
	my $boundaryMatrix0 = boundaryMap($allnmin1rMatroids,$allnmin2rMatroids,$n-2,$r,$allnmin1rMatroids->size,$allnmin2rMatroids->size,1);
	
	#Recast as sparse for input into topaz's homology.
	$boundaryMatrix0 = new SparseMatrix<Integer>($boundaryMatrix0);
	$boundaryMatrix1 = new SparseMatrix<Integer>($boundaryMatrix1);

	my $arrayBoundaryMatrices = new Array<SparseMatrix<Integer>>($boundaryMatrix0,$boundaryMatrix1);
	my $chainComplex = new topaz::ChainComplex<SparseMatrix<Integer>>($arrayBoundaryMatrices,0);
	print "\n";
	#The output can be used to find representatives of each homology class.
	print topaz::homology_and_cycles($chainComplex,0)->[1];
	print "\n";
	return;
}

sub main {
	my $m = shift @ARGV;
	my $c = shift @ARGV;
	my $n = $m->N_ELEMENTS;
	my $r = $m->RANK;
	my $currentChainOfFlats = new Vector<Set<Int>>();
	my $lattice = $m->LATTICE_OF_FLATS;
	my $flat = new Set<Int>();
	my $currentRank = 0;
	my $arrayMatroids = new Array<Matroid>(0);
	
	my $monomialBasisDegC =  truncatorRecursion($m,$n,$r,$c,$lattice,$flat,$arrayMatroids,$lattice->nodes_of_rank_range($currentRank,$r),$currentChainOfFlats,$currentRank);
	
	my $relativeNestedDegC = cleanUp($monomialBasisDegC);
	
	mainDel($relativeNestedDegC);
	return;
}

main();