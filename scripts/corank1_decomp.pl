#Corank 1 decomposition as laid out in Prop 5.7 of Simon Hampe's Intersection Ring of Matroids
#Implemented by Austin Alderete

use application 'matroid';
use Data::Dump qw(dump);

sub corankDecomposition {
	my $nestedPres = $_[0];
	my $lengthPres = scalar(@{$nestedPres});
	my $nestedPresVec = new Vector<Set>($lengthPres);
	for my $i (0.. $lengthPres-1){
		$nestedPresVec->[$i] = $nestedPres->[$i];
	}
	
	#This generates a presentation whose sets determine the coloops of the various corank 1 matroid
	my $existSingletonDif = true;
	COLOOPER : while($existSingletonDif){
		for my $i (0..$lengthPres-2){
			my $setDif = $nestedPresVec->[$i+1]-$nestedPresVec->[$i];
			if ($setDif->size == 1){
				$nestedPresVec->[$i] = $nestedPresVec->[$i]+$setDif;
				next COLOOPER;
			}
		}
		$existSingletonDif = false;
	}
	
	# We rename for convienience.
	my $coloopVec = $nestedPresVec;
	
	#Combinatorics calculation
	#Alternatively just keep the groundset or new Set(0..n-1)
	#There should be a much more efficient way to do this, but this was first to come to mind
	my $nSet = $nestedPres->[-1];
	my $kSubsets = all_subsets_of_k($nSet,$nSet->size-1);
	
	my $corank1List = ();
	
	for my $i (0..$coloopVec->dim-1){
		my $bases = new Set<Set>();
		if($coloopVec->[$i] == $nSet){
			push(@{$corank1List},uniform_matroid($nSet->size-1,$nSet->size));
		}
		else {
			my $coloopSet = $coloopVec->[$i];
			#The corank 1 matroid is U_|G|,G oplus U_|G^c|-1,G^c but the labels matter.
			# We look at G complement and take all subsets of right size then add back G
			my $coloopComplement = $nSet-$coloopSet;
			my $complementSubsets = all_subsets_of_k($coloopComplement,$coloopComplement->size-1);
			for my $candidateSet (@{$complementSubsets}){
				my $base = $candidateSet+$coloopSet;
				$bases+=$base;
			}
			my $corank1Matroid = new Matroid("N_ELEMENTS"=>$nSet->size,"RANK"=>$nSet->size-1,"BASES"=>$bases,"N_LOOPS"=>0);
			push(@{$corank1List},$corank1Matroid);
		}	
	}


	return $corank1List;
}

sub nestedDecomposition {
	my $m = $_[0];
	my $tropicalCycle = new tropical::MatroidRingCycle<Max>($m);
	my $arrayNestedPresentations = $tropicalCycle->NESTED_PRESENTATIONS;
	my $arrayNestedCoefficients = $tropicalCycle->NESTED_COEFFICIENTS;
	return ($arrayNestedPresentations, $arrayNestedCoefficients);
}

sub main {
	
	my $m = shift @ARGV;
	my ($arrayNestedPresentations,$arrayNestedCoefficients) = nestedDecomposition($m);
	
	#Set up to use push(@{$array},$matroid)
	my $arrayCorankPresentation = ();
	my $arrayCorankCoefficients = ();
	
	for my $i (0..$arrayNestedPresentations->size-1){
		my $corank1List = corankDecomposition($arrayNestedPresentations->[$i]);
		
		for my $corank1Matroid (@{$corank1List}){
			push(@{$arrayCorankPresentation},$corank1Matroid);
			push(@{$arrayCorankCoefficients},$arrayNestedCoefficients->[$i]);
		}
	}
	$arrayCorankPresentation = new Array<Matroid>($arrayCorankPresentation);
	$arrayCorankCoefficients = new Array<Int>($arrayCorankCoefficients);
	return ($arrayCorankPresentation,$arrayCorankCoefficients);
}

main();