#From a matroid, generates the subring of the intersection ring graded by # of factors of corank 1 product and determines whether or not the Chow ring relations (in particular the linear ones) are reflected in the intersection ring.
#
#Implemented by Austin Alderete

use application 'matroid';
use Data::Dump qw(dump);

sub moebius {
	my $lattice = $_[0];
	my $F = $_[1];
	my $G = $_[2];
	my $localMobNum = 0;
	
	#Sanity check
	if($F+$G != $G){
		return 0;
	}
	
	elsif($F == $G){
		return 1;
	}
	
	else{
		#This can be optimized by restricting to flats of rank geq to F and less than G.
		for my $flat (@{$lattice->FACES}){
			if($flat+$F == $flat and $flat+$G == $G and $flat->size < $G->size){
				 $localMobNum-=moebius($lattice,$F,$flat);
			}
		}
	}
	return $localMobNum;
}

sub indicatorVectorofCycle {
	my $cycleProduct = $_[0];
	my @nestedMatroids = $cycleProduct->nested_matroids;
	my $nestedCoef = $cycleProduct->NESTED_COEFFICIENTS;
	
	if(scalar(@nestedMatroids) == 0){
		return zero_vector(1);
	}
	
	my $indicatorVector = new Vector<Int>(zero_vector($nestedMatroids[0]->INDICATOR_VECTOR->dim));
	
	for my $index (0..scalar(@nestedMatroids)-1){
			my $tempIndicatorVector = new Vector<Int>($nestedMatroids[$index]->INDICATOR_VECTOR);
			$indicatorVector += $nestedCoef->[$index]*$tempIndicatorVector;
	}
	return $indicatorVector;
}

sub tropicalCycleConstructor {
	my $m = $_[0];
	my $n = $_[1];
	my $r = $_[2];
	my $flat = $_[3];
	my $arrayFlats = $m->LATTICE_OF_FLATS->FACES;
	
	my $tropicalCycle = ();
	
	for my $superFlat (@{$arrayFlats}){
		if($superFlat->size >= $flat->size){
			if($superFlat+$flat == $superFlat){
				my $corank1Matroid = new Matroid(N_ELEMENTS=>$n,RANK=>$n-1,CIRCUITS=>[$superFlat]);
				my $moebius= moebius($m->LATTICE_OF_FLATS,$flat,$superFlat);;
				if($tropicalCycle == ()){
					$tropicalCycle = new tropical::MatroidRingCycle<Max>($corank1Matroid);
					$tropicalCycle = $moebius*$tropicalCycle;
				}
				else{
					my $tempCycle = new tropical::MatroidRingCycle<Max>($corank1Matroid);
					$tropicalCycle += $moebius*$tempCycle;
				}
			}
		}
	}
	$tropicalCycle->attach("FLAT",$flat);
	return $tropicalCycle;
}


sub relationChecker {
	my $m = $_[0];
	my $n = $m->N_ELEMENTS;
	# my $nSet = new Set<Int>(0..$n-1);
	my $r = $m->RANK;
	
	for my $flatIndex (@{$m->LATTICE_OF_FLATS->nodes_of_rank(1)}){
		my $flatRank1 = $m->LATTICE_OF_FLATS->FACES->[$flatIndex];
		if($flatRank1->size == 1){
			print "F = ";
			print $flatRank1;
			print "  :  M(F) = 0 \n";
			print "Reason: F is a singleton, so the corank 1 matroid has a loop. \n";
		}
		else {
			print "F = ";
			print $flatRank1;
			my $tempMat = new Matroid(N_ELEMENTS=>$n,RANK=>$n-1,CIRCUITS=>[$flatRank1]);
			if(is_zero($tempMat->INDICATOR_VECTOR)){
				print " : M(F) = 0 \n";
				print "Reason: the corank 1 matroid has trivial indicator vector. \n";
			}
			else {
				print " : M(F) =/= 0 \n";
				print "Reason: the corank 1 matroid M(F) has nontrivial indicator vector. \n";
			}
		}
	}
	
	#We construct the tropical cycles \sum_{F \subseteq G} \mu(F,G) M(G) along with the data of a flat. We will then iterate over the list twice, taking incomparable sets.
	
	my $arrayTropicalCycles = new Array<tropical::MatroidRingCycle<Max>>(0);

	for my $flat (@{$m->LATTICE_OF_FLATS->FACES}){
		my $tropicalCycle = tropicalCycleConstructor($m,$n,$r,$flat);
		push @{$arrayTropicalCycles}, $tropicalCycle;
	}
	
	for my $cycle1 (@{$arrayTropicalCycles}){
		my $flat1 = $cycle1->get_attachment("FLAT");
		if($flat1->size <2){
			next;
		}
		for my $cycle2 (@{$arrayTropicalCycles}){
			my $flat2 = $cycle2->get_attachment("FLAT");
			if($flat2->size <2){
				next;
			}
			my $unionFlats = $flat1+$flat2;
			if($unionFlats->size > max($flat1->size,$flat2->size)){
				my $cycleProduct = $cycle1*$cycle2;
				my $indicatorVector = indicatorVectorofCycle($cycleProduct);
				if(is_zero($indicatorVector)){
					print "F = ";
					print $flat1;
					print " and G = ";
					print $flat2;
					print "  Quadratic relation is zero. \n";
				}
				else {
					print "F = ";
					print $flat1;
					print " and G = ";
					print $flat2;
					print "    Quadratic relation is nonzero. \n";
				}
			}
		}
	}
	
	return true;
}


sub main {
	my $m = shift @ARGV;
	my $n = $m->N_ELEMENTS;
	my $r = $m->RANK;
	
	return relationChecker($m,$n,$r);
}

main();