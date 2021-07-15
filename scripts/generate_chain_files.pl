use application 'matroid';

sub usage {
	return "Usage:\n$0 --script \"SCRIPT\" FILEPATH R N\n\tFILEPATH or -\n\t0 <= R <= N\n";
}

sub createChainFile {
	my $filePath = $_[0];
	my $r = $_[1];
	my $n = $_[2];

	my $chainMatrices = generate_chain_matrices($r, $n);
	my $chainMatricesCount = scalar(@$chainMatrices);

	# Determine whether to write to file or stdout
	my $isWriteToFile = ($filePath ne '-');
	my $FH;
	if ($isWriteToFile) {
		open($FH, '>', $filePath) or die "Bad filepath!";
	} else {
		open($FH, '>&', \*STDOUT) or die "Bad stdout!";
	}

	# Serialize chain manually since default json/data dumpers don't like Polymake outputs
	my $lineNo = 0;
	for my $chain (@$chainMatrices) {
		if ($lineNo == 0) {
			print $FH "[";
		} else {
			print $FH ",\n ";
		}
		# print the chain e.g. [[0,0],[1,1]]
		print $FH "[";
		my $linkNo = 0;
		for my $link (@$chain) {
			if ($linkNo > 0) {
				print $FH ",";
			}
			# print the link contents e.g. [0,0]
			print $FH "[";
			my $elemNo = 0;
			for my $elem (@$link) {
				if ($elemNo > 0) {
					print $FH ",";
				}
				print $FH $elem;
				$elemNo = $elemNo + 1;
			}
			print $FH "]";
			$linkNo = $linkNo + 1;
		}
		print $FH "]";
		$lineNo = $lineNo + 1;
	}
	if ($lineNo == 0) {
		print $FH "[";
	}
	print $FH "]\n";

	if ($isWriteToFile) {
		print "Wrote $chainMatricesCount chains of r=$r, n=$n to $filePath!\n";
		close (FH);
	}
}

# Check arguments
if (
	# three arguments (filePath, n, r)
	@ARGV != 3 ||
	# filePath is - or writable
	$ARGV[0] != "-" && !-w $ARGV[0] ||
	# r is int
	$ARGV[1] !~ /^\d+$/ ||
	# n is int
	$ARGV[2] !~ /^\d+$/ ||
	# r must be at least 0 and at most n
	int($ARGV[1]) < 0 || int($ARGV[1]) > int($ARGV[2])
) {
	die usage();
}

# Run
createChainFile($ARGV[0], int($ARGV[1]), int($ARGV[2]));
