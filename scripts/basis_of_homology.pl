use application 'matroid';
use Data::Dump qw(dump);

sub collectionFromPolyDB {
	my $n = $_[0];
	my $r = $_[1];
	my $matroidType = $_[2];

	#print "Loading matroids from PolyDB...\n";
	polyDB;
	$PolyDB::default::db_section_name = "Matroids";
	$PolyDB::default::db_collection_name = "Small";
	my $collection = polyDB->get_collection();
	#For the dual case ask that it have no coloops
	my $matroidsRN;
	my $matroidsCount;
	
	if($n/2>=$r){
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,$matroidType=>true});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($r), N_LOOPS => 0,$matroidType=>true});
		}
		else{
			# We use the fact that matorid classes of a fixed cyclic width are dual-closed, in particular cyclic width of one.
			$matroidsRN = $collection->find({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0,$matroidType=>true});
			$matroidsCount = $collection->count({N_ELEMENTS => int($n), RANK => int($n-$r), "DUAL.N_LOOPS" => 0, $matroidType=>true});
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
	return $matroidsRN;
}

#This generates non-equal yet isomorphic matroids. It does so by computing all bases which arise from permutation and throwing out (via collect) duplicates.
sub permNaiveAllMat {
	my $matroidsRN = $_[0];
	my $n = $_[1];
	my $r = $_[2];
	my $group = group::symmetric_group($n);
	
	my $arrayAllMatroids = ();
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
			$allPossibleBases->collect($permBaseSet);
		}
		
		for my $basesSet (@$allPossibleBases){
			my $permMatroid = new Matroid(BASES=>$basesSet,N_ELEMENTS=>$n,RANK=>$r);
			push @{$arrayAllMatroids}, $permMatroid;
		}
	}
	$arrayAllMatroids = new Array<Matroid>($arrayAllMatroids);
	return $arrayAllMatroids;		
}

sub boundaryMap {
	my $allDomainMatroids = $_[0];
	my $allRangeMatroids = $_[1];
	my $n = $_[2];
	my $r = $_[3];
	my $domainDim = $_[4];
	my $rangeDim = $_[5];
	my $conOrdel = $_[6];
	#print $domainDim;
	#print "\n";
	#print $rangeDim;
	
	my $boundaryMatrix = new Matrix<Int>($domainDim,$rangeDim);
	
	if ($conOrdel == 0){
		for (my $domCount = 0; $domCount <$domainDim;$domCount++){
			my $currentDomMatroid = $allDomainMatroids->[$domCount];
			for (my $elem = 0; $elem <$n; $elem++){
				my $conMat = contraction($currentDomMatroid,$elem);
				my $conMatBases = $conMat->BASES;
				
				for (my $ranCount = 0; $ranCount < $rangeDim ; $ranCount++){
					#print "\n";
					my $ranMatBases = $allRangeMatroids->[$ranCount]->BASES;
					if ($conMatBases == $ranMatBases ){
						$boundaryMatrix->elem($domCount,$ranCount) += (-1)**$elem;
						#print $boundaryMatrix;
						#print "\n";
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
					#print "\n";
					my $ranMatBases = $allRangeMatroids->[$ranCount]->BASES;
					if ($delMatBases == $ranMatBases ){
						$boundaryMatrix->elem($domCount,$ranCount) += (-1)**$elem;
						#print $boundaryMatrix;
						#print "\n";
					}
				}
			}
		}
	}
	
	return $boundaryMatrix;
}

sub iterateAllDel {
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
	
	#I'm not sure why this is needed. I should probably initialize it to be sparse in the first place...
	$boundaryMatrix0 = new SparseMatrix<Integer>($boundaryMatrix0);
	$boundaryMatrix1 = new SparseMatrix<Integer>($boundaryMatrix1);

	my $arrayBoundaryMatrices = new Array<SparseMatrix<Integer>>($boundaryMatrix0,$boundaryMatrix1);
	my $chainComplex = new topaz::ChainComplex<SparseMatrix<Integer>>($arrayBoundaryMatrices,0);
	print "\n";
	print topaz::homology($chainComplex,0)->[1];
	print "\n";
}

sub iterateAllCon {
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
	
	#I'm not sure why this is needed. I should probably initialize it to be sparse in the first place...
	$boundaryMatrix0 = new SparseMatrix<Integer>($boundaryMatrix0);
	$boundaryMatrix1 = new SparseMatrix<Integer>($boundaryMatrix1);

	my $arrayBoundaryMatrices = new Array<SparseMatrix<Integer>>($boundaryMatrix0,$boundaryMatrix1);
	my $chainComplex = new topaz::ChainComplex<SparseMatrix<Integer>>($arrayBoundaryMatrices,0);
	print "\n";
	print topaz::homology($chainComplex,0)->[1];
	print "\n";
}

sub main {
	# 0 = contraction, 1 = deletion
	
	for (my $n =3;$n<11;$n++){
		for (my $r=2;$r<$n;$r++){
				print "\n";
				print "n = :";
				print $n;
				print "\n r =:";
				print $r;
				print "\n";
				print "Contraction: \n";
				iterateAllCon($n,$r,0);
				print "Deletion: \n";
				iterateAllDel($n,$r,1);
		}
	}
}

main()