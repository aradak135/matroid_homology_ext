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
		#,NESTED=>true
		,REGULAR=>true
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
		if(-1 > $matroid->N_CONNECTED_COMPONENTS){
			if($matroidsRN->PolyDB::Cursor::has_next){
				next;
			}
			else{
				last;
			}
		}
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

# This functions as main for deletion.
sub mainDel {
	my $n = $_[0];
	my $r = $_[1];
	my $conOrdel = $_[2];
	my $nrMatroids = collectionFromPolyDB($n,$r,"NESTED");
	my $nmin1rMatroids = collectionFromPolyDB($n-1,$r,"NESTED");
	my $nplus1rMatroids = collectionFromPolyDB($n+1,$r,"NESTED");
	
	my $allnrMatroids = permNaiveAllMat($nrMatroids,$n,$r);
	my $allnmin1rMatroids = permNaiveAllMat($nmin1rMatroids,$n-1,$r);
	my $allnplus1rMatroids = permNaiveAllMat($nplus1rMatroids,$n+1,$r);
	
	my $boundaryMatrix0 = boundaryMap($allnrMatroids,$allnmin1rMatroids,$n,$r,$allnrMatroids->size,$allnmin1rMatroids->size,$conOrdel);
	my $boundaryMatrix1 = boundaryMap($allnplus1rMatroids,$allnrMatroids,$n+1,$r,$allnplus1rMatroids->size,$allnrMatroids->size,$conOrdel);
	
	#Recast as sparse for input into topaz's homology.
	$boundaryMatrix0 = new SparseMatrix<Integer>($boundaryMatrix0);
	$boundaryMatrix1 = new SparseMatrix<Integer>($boundaryMatrix1);

	my $arrayBoundaryMatrices = new Array<SparseMatrix<Integer>>($boundaryMatrix0,$boundaryMatrix1);
	my $chainComplex = new topaz::ChainComplex<SparseMatrix<Integer>>($arrayBoundaryMatrices,0);
	print "\n";
	#The output can be used to find representatives of each homology class.
	print topaz::homology_and_cycles($chainComplex,0)->[1];
	print "\n";
}

# This functions as main for contraction.
sub mainCon {
	my $n = $_[0];
	my $r = $_[1];
	my $conOrdel = $_[2];
	my $nrMatroids = collectionFromPolyDB($n,$r,"NESTED");
	my $nmin1rmin1Matroids = collectionFromPolyDB($n-1,$r-1,"NESTED");
	my $nplus1rplus1Matroids = collectionFromPolyDB($n+1,$r+1,"NESTED");
	
	my $allnrMatroids = permNaiveAllMat($nrMatroids,$n,$r);
	my $allnmin1rmin1Matroids = permNaiveAllMat($nmin1rmin1Matroids,$n-1,$r-1);
	my $allnplus1rplus1Matroids = permNaiveAllMat($nplus1rplus1Matroids,$n+1,$r+1);
	
	my $boundaryMatrix0 = boundaryMap($allnrMatroids,$allnmin1rmin1Matroids,$n,$r,$allnrMatroids->size,$allnmin1rmin1Matroids->size,$conOrdel);
	my $boundaryMatrix1 = boundaryMap($allnplus1rplus1Matroids,$allnrMatroids,$n+1,$r+1,$allnplus1rplus1Matroids->size,$allnrMatroids->size,$conOrdel);
	
	$boundaryMatrix0 = new SparseMatrix<Integer>($boundaryMatrix0);
	$boundaryMatrix1 = new SparseMatrix<Integer>($boundaryMatrix1);

	my $arrayBoundaryMatrices = new Array<SparseMatrix<Integer>>($boundaryMatrix0,$boundaryMatrix1);
	my $chainComplex = new topaz::ChainComplex<SparseMatrix<Integer>>($arrayBoundaryMatrices,0);
	print "\n";
	#The output can be used to find representatives of each homology class.
	print topaz::homology_and_cycles($chainComplex,0)->[1];
	print "\n";
}

sub main {
	# 0 = contraction, 1 = deletion
	polyDB;
	
	print "Lower bound for n (inclusive): ";
	my $nLow = <>;
	chomp($nLow);
	print "\n Upper bound for n (inclusive): ";
	my $nHigh = <>;
	chomp($nHigh);
	print "\n";
	
	#Iterate over all.
	for (my $n =$nLow;$n<$nHigh+1;$n++){
		for (my $r=2;$r<$n-1;$r++){
				print "\n";
				print "n = :";
				print $n;
				print "\n r =:";
				print $r;
				print "\n";
				print "Contraction: \n";
				mainCon($n,$r,0);
				print "Deletion: \n";
				mainDel($n,$r,1);
		}
	}
}

main()