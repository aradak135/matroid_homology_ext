#Tutte's algorithm for testing graphicness
#Implemented by Austin Alderete

use application 'matroid';
use Data::Dump qw(dump);

sub isYeven {
	my $m = $_[0];
	my $deleteCols = $_[1];
	my $bridgesY = $_[2];
	my $bridgesYVector = new Vector<Set<Int>>($bridgesY);
	
	my $groundSet = new Set<Int>(0..$m->N_ELEMENTS-1);
	
	#This is a 2 coloring problem.
	#Split into two partition: A = -1, B = 1, unallocated =0;
	my $colors = new Vector<Int>(zero_vector($bridgesYVector->dim));
	my $numColored = 0;
	
	while($numColored < $bridgesYVector->dim){
		my $nextUncolored = 0;
		for my $i (0..$bridgesYVector->dim-1){
			if($colors->[$i] == 0){
				$nextUncolored = $i;
				last;
			}
		}
		my $currentBridge = $bridgesYVector->[$nextUncolored];
		#Color it 1
		$colors->[$nextUncolored] = 1;
		$numColored+=1;
		my $adjacentSetsIndex = new Set<Int>();
		$adjacentSetsIndex += $nextUncolored;
		
		while($adjacentSetsIndex->size>0){
			my $adjIndex = $adjacentSetsIndex->front;
			$adjacentSetsIndex -= $adjIndex;
			for my $j (0..$bridgesYVector->dim-1){
				if($bridgesYVector->[$adjIndex]*$bridgesYVector->[$j] != new Set<Int>() and $adjIndex != $j){
					if($colors->[$adjIndex] == $colors->[$j]){
						return false;
					}
					elsif($colors->[$j] == 0){
						$colors->[$j] = -1*$colors->[$adjIndex];
						$numColored += 1;
						$adjacentSetsIndex += $bridgesYVector->[$j];
					}
				}
			}
		}
	}
	return true;
}

sub bridges {
	my $m = $_[0];
	my $repMatrix = $_[1];
	my $colIndex = $_[2];
	my $rowIndex = $_[3];

	my $deleteCols = new Set();
	for my $j (0..$repMatrix->cols-1){
		if($repMatrix->elem($rowIndex,$j)){
			$deleteCols += $j;
		}
	}
	
	# print "\nFor the row ";
	# print $rowIndex;
	# print " and the column ";
	# print $colIndex;
	# print " the matrix contains 1s in the columns ";
	# print $deleteCols;
	# print ".\n";
	
	my $groundSet = new Set(0..$m->N_ELEMENTS-1);
	my $groundLabels= new Vector<Int>($groundSet);
	my $deletionLabels = new Vector<Int>($groundSet-$deleteCols);
	
	my $mY = deletion($m,$deleteCols);
	
	if($mY->N_CONNECTED_COMPONENTS == 1){
		return false;
		}
	
	my $componentsY =  new Set<Set>($mY->CONNECTED_COMPONENTS);
	
	my $bridgesY = new Set<Set>();
	for my $component (@{$componentsY}){
		my $bridgeY = new Set<Int>();
		for my $index (@{$component}){
			$bridgeY += $deletionLabels->[$index];
		}
		$bridgesY+=$bridgeY;
	}
	
	if(isYeven($m,$deleteCols,$bridgesY)==false){
		return false;
	}
	
	for my $bridgeY (@{$bridgesY}){
		my $bridgeMatroid = contraction($m,$groundSet-($deleteCols+$bridgeY));
		#print transpose($bridgeMatroid->BINARY_VECTORS);
		#print "\n";
		
		if(binaryGraphic($bridgeMatroid) == false){
			return false;
		}
	}
	
	return true;
}

sub binaryGraphic {
	#Takes as input a connected binary matroid
	my $m = $_[0];
	
	my $repMatrix = transpose($m->BINARY_VECTORS);
	
	# print "Binary vectors of current component: \n";
	# print $repMatrix;
	# print "\n";
	
	my $mod2Check = 1;
	my $colIndex;
	for my $i (0..$repMatrix->cols-1){
		my $col = $repMatrix->col($i);
		my $colSum = 0;
		for (@{$col}){
			$colSum += $_;
		}
		if($colSum>2){
			$mod2Check = 0;
			$colIndex = $i;
			last;
		}
	}
	if ($mod2Check){
	#	print "\nTrue by virtue of mod2 check. \n";
		return true;
	}
	
	my $rowIndices = new Set<Int>();
	for my $i (0..$repMatrix->rows-1){
		if($repMatrix->col($colIndex)->[$i] == 1){
			$rowIndices += $i;
		}
	}
	
	# print "\nColumn ";
	# print $colIndex;
	# print " has 1s in Rows ";
	# print $rowIndices;
	# print ". \n";
	
	my $componentCounter = 0;
	for my $rowIndex (@{$rowIndices}){
		if(bridges($m,$repMatrix,$colIndex,$rowIndex)){
			return true;
		}
		$componentCounter +=1;
		if($componentCounter ==3){
			return false;
		}
	}
}

sub connectedComponentBreakdown {
	# Takes as input a matroid
	my $m = $_[0];
	if ($m->CONNECTED){
		#print "\nMatroid is connected. \n";
		return binaryGraphic($m);
	}
	my $groundSet = new Set(0..$m->N_ELEMENTS-1);
	
	for my $component (@{$m->CONNECTED_COMPONENTS}){
		#print "\nMatroid has ";
		#print $m->N_CONNECTED_COMPONENTS;
		#print ". \n";
		my $restMatroid = deletion($m,$groundSet-$component);
		if(binaryGraphic($restMatroid)==false){
			return false;
		}
	}
	return true;
}

sub main {
	
	my $m = shift @ARGV;
	
	#my $m = deletion(fano_matroid(),2);
		
	#TESTING
	#my $repMatrix = new Matrix<Int>([[1,0,0,0,0,0,0,1,1,0,0,0,0,0,1],[0,1,0,0,0,0,0,1,0,1,0,0,0,1,0],[0,0,1,0,0,0,0,1,0,1,1,0,0,1,0],[0,0,0,1,0,0,0,1,1,1,1,1,0,1,0],[0,0,0,0,1,0,0,1,1,1,1,1,1,0,0],[0,0,0,0,0,1,0,1,1,1,0,1,1,0,0],[0,0,0,0,0,0,1,1,1,0,0,1,1,0,1]]);
	#$m = new Matroid(BINARY_VECTORS=>transpose($repMatrix));
	
	#Regularity here implies that each bridge of Y in M partitions Y.
	if ($m->REGULAR){
		return connectedComponentBreakdown($m);
	}
	else {
		return false;
	}
}

main();