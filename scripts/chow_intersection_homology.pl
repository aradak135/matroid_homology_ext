#From a matroid, generates the subring of the intersection ring graded by # of factors of corank 1 product and computes the homology imduced by deletion (which is where we've observed nontrivial homology
#
#Implemented by Austin Alderete

use application 'matroid';
use Data::Dump qw(dump);

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
	
	#We only want the nested matroids here.
	my $matroidsRN = $collection->find({
		N_ELEMENTS => int($n),
		RANK => int($r),
		$loopColoop => 0
		,NESTED=>true
		});
		
	if($matroidsRN->PolyDB::Cursor::has_next){
		return $matroidsRN;
	}
	else{	
	die("Class of matroids is empty or not polyDB supported. Aborting program. \n");
	}
}

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
		my $allPossibleBases = new Set<Set<Set<Int>>>();
		
		for my $perm(@{group::all_group_elements($group->PERMUTATION_ACTION)}){
			my $permBases = [];
			for my $base (@{$matroid->BASES}){
				push @{$permBases}, group::action($perm,$base);
			}
			my $permBaseSet = new Set<Set<Int>>($permBases);
			#This avoids duplicates.
			$allPossibleBases->collect($permBaseSet);
		}
		
		for my $basesSet (@$allPossibleBases){
			my $permMatroid = new Matroid(BASES=>$basesSet,N_ELEMENTS=>$n,RANK=>$r);
			push @{$arrayAllMatroids}, $permMatroid;
		}
	}
	return $arrayAllMatroids;		
}

sub circuitMatroidHGConstructor {
	my $m = $_[0];
	my $n = $m->N_ELEMENTS;
	# my $nSet = new Set<Int>(0..$n-1);
	my $r = $m->RANK;
	
	#Due to the fact that singletons induce loops, we cannot preallocate using numProperFlats.
	my $circuitMatroidHGList = new Array<Matroid>(0);
	for my $flat (@{$m->LATTICE_OF_FLATS->FACES}){
		if($flat->size >1 and $flat->size != $n){
			my $tempMatroidHG = new Matroid(N_ELEMENTS=>$n,RANK=>$n-1,CIRCUITS=>[$flat]);
			push @{$circuitMatroidHGList}, $tempMatroidHG;
		}
	}
	
	return $circuitMatroidHGList;
}

sub intersectionRingGenerator {
	my $m = $_[0];
	my $circuitMatroidHGList = $_[1];
	
	#We construct the array graded by CORANK and populate it inductively.
	#The syntax here is troublesome
		# @listOfGradedComponents = ();
		# @arrayOfGeneratorsForAComponent = new Array<Matroid>(size)
		# populate @arrayOfGeneratorsForAComponent
		# push @listOfGradedComponents, \@arrayOfGeneratorsForAComponent
	# Then to call do $listOfGradedComponents[corank]->[index]
	
	my @gradedComponents = ();
	#We already have the corank 1 matroids
	push @gradedComponents, \@{$circuitMatroidHGList};
	# currentCorank indicates the corank of the last completed entry.
	my $currentCorank = 1;
	
	while($currentCorank < $m->N_ELEMENTS-1){
		my $recentComponent = $gradedComponents[-1];
		my $currentComponent = new Array<Matroid>(0);
		for my $matroidCorankK (@{$recentComponent}){
			for my $matroidCorank1 (@{$circuitMatroidHGList}){
				my $tempMatroidProduct = dual(union(dual($matroidCorank1),dual($matroidCorankK)));
				if($tempMatroidProduct->N_LOOPS == 0){
					#We need to check for duplicates. We have a few different options, but as the indicator vector of a matroid completely determines it and we will need those later, we compute those as they're stored with the matroid
					
					my $duplicateCheck = 0;
					for my $tempMatroid (@{$currentComponent}){
						if($tempMatroidProduct->INDICATOR_VECTOR == $tempMatroid->INDICATOR_VECTOR){
							$duplicateCheck =1;
						}
					}
					if($duplicateCheck ==0){
						push @{$currentComponent}, $tempMatroidProduct;
					}
				}
			}
		}
		push @gradedComponents,\@{$currentComponent};
		$currentCorank +=1;
	}
	return @gradedComponents;
}

sub basisConstructor {
	my $allMatroids = $_[0];
	my $basisMatroids = new Array<Matroid>(0);
	
	my $matrixIndicators = new Matrix<Integer>();
	
	#This is a  temporary fix to avoid syntax errors
	for my $matroid (@{$allMatroids}){
		my $indicatorVector = new Vector<Integer>($matroid->INDICATOR_VECTOR);
		$matrixIndicators = new Matrix<Integer>($indicatorVector->dim,0);
		last;
	}
	
	for my $matroid (@{$allMatroids}){
		my $indicatorVector = new Vector<Integer>($matroid->INDICATOR_VECTOR);
		
		my $matrixIndicatorsAppend = new Matrix<Integer>($matrixIndicators|$indicatorVector);
		if(is_zero(null_space_integer($matrixIndicatorsAppend))){
			$matrixIndicators = $matrixIndicatorsAppend;
			push @{$basisMatroids},$matroid;
		}
	}
	return $basisMatroids;
}

sub boundaryMapAB {
	my $allDomainMatroids = $_[0];
	my $allRangeMatroids = $_[1];
	my $n = $_[2];
	my $r = $_[3];
	
	my $boundaryMatrix = new Matrix<Int>($allDomainMatroids->size,$allRangeMatroids->size);
	
	for (my $domCount = 0; $domCount < $allDomainMatroids->size; $domCount++){
		my $domMatroid = $allDomainMatroids->[$domCount];
		my $tropicalCycle = new tropical::MatroidRingCycle<Max>($domMatroid);
		
		my @arrayNestedMatroids = $tropicalCycle->nested_matroids;
		my $arrayNestedCoefficients = $tropicalCycle->NESTED_COEFFICIENTS;
		
		for (my $nestIndex =0; $nestIndex < $arrayNestedCoefficients->size; $nestIndex++){
			my $nestedMatroid = $arrayNestedMatroids[$nestIndex];
			my $nestCoef = $arrayNestedCoefficients->[$nestIndex];
			for (my $elem =0; $elem<$n;$elem++){
				my $delMatroid = deletion($nestedMatroid,$elem);
				for (my $ranCount =0; $ranCount < $allRangeMatroids->size; $ranCount++){
					if($delMatroid->BASES == $allRangeMatroids->[$ranCount]->BASES){
						$boundaryMatrix->elem($domCount,$ranCount) += $nestIndex*(-1)**$elem;
					}
				}
			}
		}
	}
	return $boundaryMatrix;
}

sub boundaryMapBC {
	my $allDomainMatroids = $_[0];
	my $allRangeMatroids = $_[1];
	my $n = $_[2];
	my $r = $_[3];
	
	my $boundaryMatrix = new Matrix<Int>($allDomainMatroids->size,$allRangeMatroids->size);
	
	for (my $domCount = 0; $domCount < $allDomainMatroids->size; $domCount++){
		my $domMatroid = $allDomainMatroids->[$domCount];
		for (my $elem =0; $elem<$n; $elem++){
			my $delMatroid = deletion($domMatroid,$elem);
			for (my $ranCount =0; $ranCount < $allRangeMatroids->size; $ranCount++){
				if($delMatroid->BASES == $allRangeMatroids->[$ranCount]->BASES){
					$boundaryMatrix->elem($domCount,$ranCount) += (-1)**$elem;
				}
			}
		}
	}
	return $boundaryMatrix;
}

sub localHomologyDeletion {
	my $A = $_[0];
	my $n = $_[1];
	my $r = $_[2];
	
	if($n-2<$r){
		return true;
	}
	
	#$A = basisConstructor($A);
		
	my $nmin1rMatroids = collectionFromPolyDB($n-1,$r);
	my $nmin2rMatroids = collectionFromPolyDB($n-2,$r);
	
	my $allnmin1rMatroids = permNaiveAllMat($nmin1rMatroids,$n-1,$r);
	my $allnmin2rMatroids = permNaiveAllMat($nmin2rMatroids,$n-2,$r);
	
	
	my $conOrdel = 1;
	
	my $boundaryMatrix0 = boundaryMapAB($A,$allnmin1rMatroids,$n,$r);
	my $boundaryMatrix1 = boundaryMapBC($allnmin1rMatroids,$allnmin2rMatroids,$n,$r);
	
	$boundaryMatrix0 = new SparseMatrix<Integer>($boundaryMatrix0);
	$boundaryMatrix1 = new SparseMatrix<Integer>($boundaryMatrix1);

	my $arrayBoundaryMatrices = new Array<SparseMatrix<Integer>>($boundaryMatrix1,$boundaryMatrix0);
	my $chainComplex = new topaz::ChainComplex<SparseMatrix<Integer>>($arrayBoundaryMatrices,1);
	
	print "\n";
	#The output can be used to find representatives of each homology class.
	print topaz::homology_and_cycles($chainComplex,0)->[1];
	print "\n";
	
	return true;
}

sub main {
	my $m = shift @ARGV;
	my $n = $m->N_ELEMENTS;
	my $r = $m->RANK;
	
	if($r <3){
		print "Rank must be at least 3. \n";
		return true;
	}
	
	my $circuitMatroidHGList = circuitMatroidHGConstructor($m);	
	my @gradedComponents = intersectionRingGenerator($m,$circuitMatroidHGList);
	
	#Testing that the output is correct
	#for my $array (@gradedComponents){
	#	for my $matroidTest (@{$array}){
	#		print $matroidTest->RANK;
	#	}
	#}
	#-------------------
	
	for my $rankIndex (1..scalar(@gradedComponents)-2){
		#It is important to note that @gradedComponents starts at corank 1. We go through the list backward as the end of the list is always rank 1 even though the lists are of different sizes. That is, -$rankIndex is actually the rank.
		# This indexing is the correct one, so any errors are elsewhere.
		#See localHomologyDeletion for an explanation of the notation.
		
		my $A = $gradedComponents[-$rankIndex];
		localHomologyDeletion($A,$n,$rankIndex);
	}
	return true;
}

main();