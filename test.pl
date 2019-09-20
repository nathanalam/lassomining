# reads the fasta files and extacts protein descriptions and sequences
sub readFASTA {
	my $name = shift;
	my $cleanspace = shift;
	$cleanspace = 0 if not defined $cleanspace;
	my @prot;
	my @seq;
	my @seql;
	# print "\n opening $name \n";
	if (open (FASTA, "<$name")) { 
		my $count = -1;
		while (my $myline = <FASTA>) {
			chomp $myline;
			if ($myline =~ m/^>/) {
				$count ++;
				$prot[$count] = $myline;
				$seq[$count] = "";
			} else {
				if ($cleanspace) {$myline =~ s/\s//gi;} 
				if ($count < 0) {print "error: $myline \n";}
				$seq[$count] .= $myline;
			}
		}
		close FASTA;
		foreach (0..$#seq) {
			$seql[$_] = length($seq[$_]);
		}
		return \@prot, \@seq, \@seql;
	}
	else {return '0', '0', '0';}
}

readFASTA(shift)