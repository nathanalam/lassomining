#!/usr/bin/perl -w

# In this version each time MAST needs to be run on a putative maturation enzyme in a particular genome, a lookup in the database is performed to check if this protein has been analyzed by MAST before and if it has then use the values from before
# In this version traseq will be used instead of getorf for finding precursors
# getorf is still used for finding neighbors
# Like v4 except that multiple proteins of the same sequence don't cause wrong locations of maturation enzymes to be reported
# Another change is that the pattern needs to be adjusted to [MVL] instead of ^ in the beginning
# This version fixes the problem of having a useless %AME hash and also of erasing sequences from %AME_scores
# Precursor pattern is output into the log file
# ORF searching behavior is changed from stop-to-stop to [MVL]-to-stop
# Clusters of precursors are saved in clusters.txt
# Warning: transeq doesn't label the -1, -2, -3 frames sequentially. Sometimes it is -1, -3, -2, sometimes some other
# combination. This means the precursor locations on the reverse strand are off by one sometimes
# Take note that rank_hits expects 4 motifs for the B enzyme and 3 motifs for the C enzyme.  Adjust accordingly.
# single quotation marks 	

use Storable;
use Getopt::Long;
use List::Util;
use Cwd;
use Time::HiRes qw( time );
use POSIX qw( floor ceil);
use List::Util qw(shuffle);

GetOptions (
		"ddir=s"  => \my $ddir,
		"nstart=i"  => \my $nstart,   # if ddir, directory number from which to start. First directory = 0
		"nstop=i"   => \my $nstop,    # directory at which to stop
		"motifsdir=s" => \our $motifsdir,
		"outdir=s" => \our $out_dir,
		"list=s" => \our $organism_list);
# nstart and nstop are useful for handling a large number of genomes piecewise

#our $pattern = '[MVL].{5,43}T.G.{6,10}[DE].{5,16}$';
#our $pattern = '[M].{5,15}T.G.{6,10}[DE].{2,12}$';
#our $pattern = '[MVL].{5,45}T.[GASC].{5,10}[DE].{5,20}$'; # . dot means anything, $  The end of the line or string, [MVL] either or 
our $pattern = '^M.{15,45}T..{6,8}[DE].{5,30}$'; # change leader peptide length to 15
our $bins = 100;# system('PATH="$PATH:/tigress/cyso/emboss/bin"');

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

# looks at each precursor sequence and tries to match to the pattern set above
sub pattern_match {
	my @prot = @{$_[0]};
	my @seq = @{$_[1]};
	my @match_prot = ();
	my @match_seq  = ();
	our $pattern;
	foreach (0..$#prot) {
		if ($seq[$_] =~ m/$pattern/) {
			push @match_prot, $prot[$_];
			push @match_seq, $seq[$_];
		}
	}
	return \@match_prot, \@match_seq;
}

sub rank_hits {
	our $num_mast;
	my $currentdir = getcwd;
	our $start_dir;
	our $out_dir;
	our %AME_scores;
	chdir "$start_dir/$out_dir";
	my @neighbors = @{$_[0]};
	my @neighbor_seqs = @{$_[1]};
	my ($Bseq, $Cseq,) = qw(none none);
	my ($Bloc, $Cloc) = ("none\tnone\tNA", "none\tnone\tNA");
	my ($final_B_rank, $final_C_rank) = (0, 0);
	my ($start, $stop);
	foreach (0..$#neighbors) {
		my $seq = $neighbor_seqs[$_];	
		my ($rankB, $rankC) = (0, 0);
		if (defined($AME_scores{$seq})) {                     #check if a neighbor has been analyzed
			$rankB = $AME_scores{$seq}[0];
			$rankC = $AME_scores{$seq}[2];			
			if ($rankB > $final_B_rank) { 
					$final_B_rank = $rankB;
					$Bseq = $seq;
					$Bloc = $neighbors[$_];
			}
			if ($rankC > $final_C_rank) { 
					$final_C_rank = $rankC;
					$Cseq = $seq;
					$Cloc = $neighbors[$_];
			}
		} else {                    #otherwise perform the full analysis
				$AME_scores{$seq} = [0,'none',0,'none']; # rankB, Bloc, rankC, Cloc
				open TEMPFASTA, (">temp_fasta.fasta") or die "can't open fasta";
				print TEMPFASTA "> ".$neighbors[$_]."\n";
				print TEMPFASTA $neighbor_seqs[$_]."\n";
				close TEMPFASTA;
				my @MASTB = `mast "$motifsdir/bmotifsoops.txt" temp_fasta.fasta -ev 10 -remcorr -hit_list -nostatus`;
				$num_mast ++;				
				my @MASTC = `mast "$motifsdir/cmotifsoops.txt" temp_fasta.fasta -ev 10 -remcorr -hit_list -nostatus`;
				$num_mast ++;				
				`rm temp_fasta.fasta`;	
				if ($#MASTB > 2) {
					my @bmotifs_present;
					for (2..$#MASTB-1) {
						my @which_bmotif = split('\s', $MASTB[$_]);
						push @bmotifs_present, $which_bmotif[1];
					}
					for my $mnum (1, 2, 3, 4) {
						if (grep ($_ =~ m/$mnum/, @bmotifs_present) ) {$rankB += 1}
					}
					$AME_scores{$seq}[0] = $rankB; # add rankB
					$AME_scores{$seq}[1] = $neighbors[$_]; # add Bloc				
					if ($rankB > $final_B_rank) { 
						$final_B_rank = $rankB;
						$Bseq = $seq;
						$Bloc = $neighbors[$_];
					}
				}
				if ($#MASTC > 2) {
					my @cmotifs_present;
					for (2..$#MASTC-1) {
						my @which_cmotif = split('\s', $MASTC[$_]);
						push @cmotifs_present, $which_cmotif[1];
					}
					for my $mnum (1, 2, 3) {
						if (grep ($_ =~ m/$mnum/, @cmotifs_present) ) {$rankC += 1}
					}
					$AME_scores{$seq}[2] = $rankC; # add rankC
					$AME_scores{$seq}[3] = $neighbors[$_]; # add Cloc				
					if ($rankC > $final_C_rank) { 
						$final_C_rank = $rankC;
						$Cseq = $seq;
						$Cloc = $neighbors[$_];
					}
				}
		}
	}	
	my $totrank = $final_B_rank + $final_C_rank;
	chdir $currentdir;
	return $totrank, $final_B_rank, $Bseq, $Cseq, $Bloc, $Cloc;
}

sub rank_hits2 {
	our $num_mast;
	my $currentdir = getcwd;
	our $start_dir;
	our $out_dir;
	our %AME_scores_2;
	chdir "$start_dir/$out_dir";
	my @neighbors = @{$_[0]};
	my @neighbor_seqs = @{$_[1]};
	my ($Dseq, $Eseq,) = qw(none none);
	my ($Dloc, $Eloc) = ("none\tnone\tNA", "none\tnone\tNA");
	my ($final_D_rank, $final_E_rank) = (0, 0);
	my ($start, $stop);
	foreach (0..$#neighbors) {
		my $seq = $neighbor_seqs[$_];	
		my ($rankD, $rankE) = (0, 0);
		if (defined($AME_scores_2{$seq})) {                     #check if a neighbor has been analyzed
			$rankD = $AME_scores_2{$seq}[0];
			$rankE = $AME_scores_2{$seq}[2];
			if ($rankD > $final_D_rank) { 
					$final_D_rank = $rankD;
					$Dseq = $seq;
					$Dloc = $neighbors[$_];
			}
			if ($rankE > $final_E_rank) { 
					$final_E_rank = $rankE;
					$Eseq = $seq;
					$Eloc = $neighbors[$_];
			}
		} else {                    #otherwise perform the full analysis
				$AME_scores_2{$seq} = [0,'none',0,'none']; # rankB, Bloc, rankC, Cloc
				open TEMPFASTA, (">temp_fasta.fasta") or die "can't open fasta";
				print TEMPFASTA "> ".$neighbors[$_]."\n";
				print TEMPFASTA $neighbor_seqs[$_]."\n";
				close TEMPFASTA;
				my @MASTD = `mast "$motifsdir/dmotifsoops.txt" temp_fasta.fasta -ev 10 -remcorr -hit_list -nostatus`;
				my @MASTE = `mast "$motifsdir/emotifsoops.txt" temp_fasta.fasta -ev 10 -remcorr -hit_list -nostatus`;
				$num_mast = $num_mast + 2;				
				`rm temp_fasta.fasta`;	
				if ($#MASTD > 2) {
					my @dmotifs_present;
					for (2..$#MASTD-1) {
						my @which_dmotif = split('\s', $MASTD[$_]);
						push @dmotifs_present, $which_dmotif[1];
					}
					for my $mnum (1, 2, 3) {
						if (grep ($_ =~ m/$mnum/, @dmotifs_present) ) {$rankD += 1}
					}
					$AME_scores_2{$seq}[0] = $rankD; # add rankB
					$AME_scores_2{$seq}[1] = $neighbors[$_]; # add Bloc				
					if ($rankD > $final_D_rank) { 
						$final_D_rank = $rankD;
						$Dseq = $seq;
						$Dloc = $neighbors[$_];
					}
				}
				if ($#MASTE > 2) {
					my @emotifs_present;
					for (2..$#MASTE-1) {
						my @which_emotif = split('\s', $MASTE[$_]);
						push @emotifs_present, $which_emotif[1];
					}
					for my $mnum (1, 2, 3) {
						if (grep ($_ =~ m/$mnum/, @emotifs_present) ) {$rankE += 1}
					}
					$AME_scores_2{$seq}[2] = $rankE; # add rankB
					$AME_scores_2{$seq}[3] = $neighbors[$_]; # add Bloc				
					if ($rankE > $final_E_rank) { 
						$final_E_rank = $rankE;
						$Eseq = $seq;
						$Eloc = $neighbors[$_];
					}
				}
		}
	}
	if ($final_D_rank < 2) { #need to change this probably. 
	my ($Dseq, $Dloc) = ("none", "none\tnone\tNA");}
	if ($final_E_rank < 2) {
	my ($Eseq, $Eloc) = ("none", "none\tnone\tNA");}
	chdir $currentdir;
	return $Dseq, $Eseq, $Dloc, $Eloc;
}

sub rank_hits3 {
	our $num_mast;
	my $currentdir = getcwd;
	our $start_dir;
	our $out_dir;
	our %AME_scores_3;
	chdir "$start_dir/$out_dir";
	my @neighbors = @{$_[0]};
	my @neighbor_seqs = @{$_[1]};
	my ($B1seq, $B1loc) = ("none", "NA\tNA\tNA");
	my $final_B1_rank = 0;
	my ($start, $stop);
	foreach (0..$#neighbors) {
		my $seq = $neighbor_seqs[$_];	
		my $rankB1 = 0;
		if (defined($AME_scores_3{$seq})) {                     #check if a neighbor has been analyzed
			$rankB1 = $AME_scores_3{$seq}[0];
			$rankB2 = $AME_scores_3{$seq}[2];	
			if ($rankB1 > $final_B1_rank) { 
					$final_B1_rank = $rankB1;
					$B1seq = $seq;
					$B1loc = $neighbors[$_];
			}
			if ($rankB2 > $final_B2_rank) { 
					$final_B2_rank = $rankB2;
					$B2seq = $seq;
					$B2loc = $neighbors[$_];
			}
			}
		else {                    #otherwise perform the full analysis
				$AME_scores_3{$seq} = [0,'none'];
				open TEMPFASTA, (">temp_fasta.fasta") or die "can't open fasta";
				print TEMPFASTA "> ".$neighbors[$_]."\n";
				print TEMPFASTA $neighbor_seqs[$_]."\n";
				close TEMPFASTA;
				my @MASTB1 = `mast "$motifsdir/b1motifs.txt" temp_fasta.fasta -ev 10 -remcorr -hit_list -nostatus`;
				my @MASTB2 = `mast "$motifsdir/b2motifs.txt" temp_fasta.fasta -ev 10 -remcorr -hit_list -nostatus`;
				$num_mast ++;				
				`rm temp_fasta.fasta`;	
				if ($#MASTB1 > 2) {
					my @b1motifs_present;
					for (2..$#MASTB1-1) {
						my @which_b1motif = split('\s', $MASTB1[$_]);
						push @b1motifs_present, $which_b1motif[1];
					}
					for my $mnum (1) {
						if (grep ($_ =~ m/$mnum/, @b1motifs_present) ) {$rankB1 += 1}
					}
					$AME_scores_3{$seq}[0] = $rankB1; # add rankB
					$AME_scores_3{$seq}[1] = $neighbors[$_]; # add Bloc				
					if ($rankB1 > $final_B1_rank) { 
						$final_B1_rank = $rankB1;
						$B1seq = $seq;
						$B1loc = $neighbors[$_];
					}
				}
				if ($#MASTB2 > 2) {
					my @b2motifs_present;
					for (2..$#MASTB2-1) {
						my @which_b2motif = split('\s', $MASTB2[$_]);
						push @b2motifs_present, $which_b2motif[1];
					}
					for my $mnum (1, 2, 3) {
						if (grep ($_ =~ m/$mnum/, @b2motifs_present) ) {$rankB2 += 1}
					}
					$AME_scores_3{$seq}[2] = $rankB2; # add rankB
					$AME_scores_3{$seq}[3] = $neighbors[$_]; # add Bloc				
					if ($rankB2 > $final_B2_rank) { 
						$final_B2_rank = $rankB2;
						$B2seq = $seq;
						$B2loc = $neighbors[$_];
					}
				}
		}
	}
	my $B2rank = $final_B1_rank + $final_B2_rank;
	chdir $currentdir;
	return $B2rank, $B1seq, $B2seq, $B1loc, $B2loc;
}

sub rank_hits_PTM {
	our $num_mast;
	my $currentdir = getcwd;
	our $start_dir;
	our $out_dir;
	our %AME_scores_2;
	chdir "$start_dir/$out_dir";
	my @neighbors = @{$_[0]};
	my @neighbor_seqs = @{$_[1]};
	my $Bloc = $_[2];
	my $Cloc = $_[3];
	my @Blocs = split("\t", $Bloc); 
	my @Clocs = split("\t", $Cloc); 
	my $BC_start = List::Util::min @Blocs, @Clocs; #something like that...
	my $BC_stop = List::Util::max @Blocs, @Clocs;
	my %starts;
	my %stops;
	my @flanking_genes; 
	my $output;
	$output = "";
	
	foreach (0..$#neighbors) {
		my @locs = split("\t", $neighbors[$_]);
		if ($locs[2] eq 'r') {
			$starts{$_} = $locs[1];
			$stops{$_} = $locs[0];
			}
		else {
			$starts{$_} = $locs[0];
			$stops{$_} = $locs[1];
		}
	}
	my @starts_sorted = sort { $starts{$a} <=> $starts{$b} } keys %starts;
	my @stops_sorted = sort { $stops{$a} <=> $stops{$b} } keys %stops;
	
	#what about overlaps, 300 bp  better 
	
	foreach $i (0..$#starts_sorted) {
		if (abs($start{$i} - $BC_stop) <= 100) {
			$BC_stop = $stop{$i};
			push @flanking_genes, $i;
			}
	}

	foreach $i (reverse 0..$#starts_sorted) {
		if (abs($BC_start - $stop{$i}) <= 100) {
			$BC_start = $start{$i};
			push @flanking_genes, $i;
			}
	}
	
	foreach $i (@flanking_genes) {
		my $seq = $neighbor_seqs[$i];	
				open TEMPFASTA, (">temp_fasta.fasta") or die "can't open fasta";
				print TEMPFASTA "> ".$neighbors[$_]."\n"; #something here is messed up
				print TEMPFASTA $neighbor_seqs[$_]."\n";
				close TEMPFASTA;
				`hmmscan --tblout out.txt "$motifsdir/Pfam-A.hmm" temp_fasta.fasta`;	
				`rm temp_fasta.fasta`;			
				open (OUT, "<out.txt") or die "woops!";
				my @out = <OUT>;
				close (OUT);
				my @outs = split(" ",$out[3]);
				my $top_match = join('', @outs[18..$#outs]);
				$neighbors[$i] =~ s/\t/\|/ig; 
				$output .= $neighbors[$i].'|'.$top_match.'||';
				print "testing".$output;
		}
	chdir $currentdir;
	return $output; # null cases
}



sub get_neighbors {
	our $bins;
	my ($neig_binned, $start, $stop) = @_;
	my %neig_binned = %$neig_binned;	
	my (@neighbors, @neighbor_seqs);
	if ($start > $stop) { ($start,$stop) = ($stop,$start) }
	$start -= 5000; $stop += 5000;
	foreach (floor($start/$bins)..ceil($stop/$bins)) {
		foreach (@{$neig_binned{$_}}) {
			my ($t1, $t2, $t3) = split("\t",@{$_}[0]);
			my $seq = @{$_}[1];
			if (($t1 >= $start && $t1 <= $stop) || ($t2 >= $start && $t2 <= $stop)) {
				push @neighbors, @{$_}[0];
				push @neighbor_seqs, $seq;
			}
		}
	}
	return \@neighbors, \@neighbor_seqs;
}

sub split_seq {
	my $seq = $_[0];
	my $frame = $_[1];
	my $seql = $_[2];
	my @letters = split('', $seq);
	my @reads = ();
	my @locs = ();
	my @chains; # This means the chain is empty
	my ($offset, $dir, $start);
	if ($frame == 0) {$dir = 1; $offset = -1; $start = 1} 
	if ($frame == 1) {$dir = 1; $offset = 0; $start = 1}
	if ($frame == 2) {$dir = 1; $offset = 1; $start = 1} 
	if ($frame == 3) {$dir = -1; $offset = 1; $start = $seql}
	if ($frame == 4) {$dir = -1; $offset = 0; $start = $seql}
	if ($frame == 5) {$dir = -1; $offset = -1; $start = $seql}
	foreach my $pos (0..$#letters) {
		# first check if you have a growing chain
		if (@chains == 0) {
			# Since no chain. Can you can start a chain?
			if ($letters[$pos] =~ m/[M]/) {$chains[0] = $letters[$pos];} # start chain
		} else {
			if ($letters[$pos] =~ m/[M]/) {
				$chains[$#chains + 1] = $letters[$pos];
			}
			if ($letters[$pos] eq 'X') {
				foreach $chain (@chains) {
					push @reads, $chain; 
					my $l = length($chain);
					my $end = $start + $offset + $pos*3*$dir; # see cal.ppt for power point on how to calculater this
					my $begin = $end - $l*3*$dir + $dir; 
					push @locs, $begin.":".$end;
				}
				@chains = (); # zero the chain so that there is no chain				
			} else {
				foreach $chain (@chains) {
					$chain = $chain.$letters[$pos];
				}
			}
		}
	}
	
	return \@reads, \@locs;
}

sub read_transeq{
	my $input = $_[0];
	my $pminl = $_[1];
	my $pmaxl = $_[2];
	my $seql = $_[3];
	my ($prot_ref, $seq_ref) = readFASTA($input);
	if ($prot_ref == 0) {return '0', '0', '0';}
	my @prot = @$prot_ref;
	my @seq = @$seq_ref;
	my (@prec, @prec_seq);
	my %groups;
	my $num_prec = 0;
	my $frame = 0; # frame will be used to figure out the location of the read
	foreach (0..$#prot) {
		my ($reads_ref, $locs_ref) = split_seq($seq[$_], $frame, $seql);	
		my @reads = @$reads_ref; my @locs = @$locs_ref;
	
		# Take note that the reads are reads between [MVL] and stop codons (X)
		
		foreach (0..$#reads) {
			my $read = $reads[$_];
			my @letters = split('', $read);
			my $l = $#letters+1;
			if ($l >= $pminl and $l <= $pmaxl) {
				my @group;
				my $chain = "";
				foreach (reverse 0..$l-1) {
					if ($letters[$_] =~ m/[M]/) {
						$chain = $letters[$_].$chain;
						if (length($chain) >= $pminl) {push @group, $chain}
					} else {
						$chain = $letters[$_].$chain;
					}
				}
				$prec[$num_prec] = $locs[$_];
				$prec_seq[$num_prec] = $group[-1];
				$groups{$group[-1]} = \@group;
				$num_prec ++;
			}
		}
	$frame++;
	}

	return \@prec, \@prec_seq, \%groups;

}

sub cluster_prec {
	my @prec = @{$_[0]};
	my @prec_seq = @{$_[1]};
	my (@prec_sort, @prec_seq_sort);
	# Find the midpoint of each precursor
	foreach (0..$#prec) {
		my ($start, $stop) = split(":", $prec[$_]);
		my $midpoint = ($stop+$start)/2;
		$prec_sort[$_] = [$midpoint, $prec[$_]];
		$prec_seq_sort[$_] = [$midpoint, $prec_seq[$_]];		
	}
	# Sort the precursor names by midpoint in ascending order
	@prec_sort = sort {$a->[0] <=> $b->[0]} @prec_sort;
	# Sort the precursor sequences names by midpoint in ascending order	
	@prec_seq_sort = sort {$a->[0] <=> $b->[0]} @prec_seq_sort;
	
	# map the sorted lists back to the original
	@prec = map {@{$prec_sort[$_]}[1]} 0..$#prec_sort;
	@prec_seq = map {@{$prec_seq_sort[$_]}[1]} 0..$#prec_seq_sort;

#   For testing purpose	
#	my @list = ();
#	foreach (@prec_sort) {
#		push @list, @{$_}[0];
#	}
	
	# initialize the location of the first cluster
	my $avg = $prec_sort[0][0];
	# initialize the first midpoint
	my $i = $prec_sort[0][0];
	# add location of first cluster to list
	my @avg = ($avg);
	# initialize the first cluster in the hash
	my $n = 0; 
	my @clusters; $clusters[$n] = [$prec[0]];
	my @clusters_seq; $clusters_seq[$n] = [$prec_seq[0]];
	foreach (1..$#prec_sort) {
		# read midpoint of next precursor
		my $j = $prec_sort[$_][0];
		if ( ($j-$i) < 1500) {
			# if distance between precursors is less than 50 cluster add jth precursor to the previous cluster
			$avg = ($avg+$j)/2;
			# calculate the new cluster location
			$avg[-1] = $avg;
			# update the new location of the cluster and add the sequence to the list
			push ${clusters[$n]}, $prec[$_];
			push ${clusters_seq[$n]}, $prec_seq[$_];				
		} else {
			# if next precursor is too far away, start a new cluster
			# new location is the midpoint of the first precursor in that cluster
			$avg = $j;
			# add the location to the end of the @avg list
			push @avg, $avg;
			# update the new cluster number
			$n ++;
			# add the sequence to that cluster
			$clusters[$n] = [$prec[$_]];
			$clusters_seq[$n] = [$prec_seq[$_]];
			$i = $j;			
		}
	}
	
	return \@clusters, \@clusters_seq, \@avg;
}


# system('PATH="$PATH:/scratch/network/cyso/emboss/bin"');

if (!$motifsdir) {$motifsdir = getcwd};

my $g_tot_count = 0;
my $g_match_count = 0;
my $g_tot_neig = 0;
our ($start_time, $end_time, @genome_dir);

$ddir = "genomes"; # ddir
my $multiple = 1; #run 1-1004 arrays

srand(1);
@array = shuffle((1..145623));
$nstart = $array[$nstart];
$nstop = $nstart;


# my $genome = 0; #145623 / 27
# my $lim = 0 ; #1000
# my $rem = $genome % $lim;
# my $k = ($genome - $rem) / $lim;
# if ($nstart == $lim) {
	# $nstop = $genome;
	# $nstart = $genome - ($rem - 1);}
# elsif ($nstart == 1) {
	# $nstop = $k;
	# $nstart = 0;}
# else {	$nstop = $nstart*$k;
	# $nstart = ($nstart-1)*$k + 1;}

if ($ddir) {
	our $num_mast = 0;
	our $num_getorf = 0;
	our $get_rank_loops = 0;
	$start_time = time();								# directory of genome directories
	`mkdir $out_dir -p`;
	our $start_dir = getcwd;
	my @genome_dir;
	if ($organism_list) {
	    open (LIST, "<$organism_list");
        @genome_dir = <LIST>;
        chomp @genome_dir;
        close LIST;
	} else {
		opendir DIR, "$ddir";
		@genome_dir = readdir DIR;						# all genome directory
		@genome_dir = grep {$_ !~ m/^\.{1,2}/} @genome_dir;
		@genome_dir = sort {$a cmp $b} @genome_dir; 
		close DIR;
	}
	chdir $ddir; 
	my $progress = 0;	# overall progress
	if (defined($nstart) and defined($nstop)) {
	} else {$nstart = 0, $nstop = $#genome_dir}
	open F3, ">>$start_dir/$out_dir/statistics.txt";
    print F3 "Number of match precursors \t Number of total precursors \t Number of neighbors \t Organism\n";
    close F3;
	
	foreach ($nstart..$nstop) {						# go into individual genome directory
		our $maindir = getcwd;						# current working directory
		my $organism = $genome_dir[$_];
		opendir DIR, "$organism";
		chdir $organism;
		if (-e "contigs.tar") { 
		print "unzipped";
		`tar -xvzf contigs.tar .`; #here
		`rm contigs.tar`;}  
		my @fna_files = readdir DIR;					# files in one genome directory
		print "Current organism: $organism\n";	
		my @fna = grep($_ =~ m/fna$/, @fna_files);			# individual chromosome records
		if (!exists $fna[0]) {
			print "\n Genome already done. \n";
			next;
		}
		my $tot_count = 0;
		my $tot_neig = 0;
		my (@g_match_rank, @g_nucl_name, @g_match_seq, # global lists are for each organism, not all the organisms
        @g_match_B, @g_match_C, @g_match_loc, @g_match_orgn, %g_groups);         # so re-zero this list after starting a new organism
        my (@match_counts, @tot_counts, @organisms, @neig_counts);    # and writing the information to disk
		my @temp = readFASTA($fna[0]);
		my @prot = @{$temp[0]};
		my @seql = @{$temp[2]};
		my @seq = @{$temp[1]};
		foreach (0 .. $#seq) {
				$prot[$_] =~ m/>(.+?)\s/;
				open NEWFNA, ">$1".'.fna';
				print NEWFNA $prot[$_];
				print NEWFNA "\n";
				print NEWFNA $seq[$_];
				close NEWFNA;
				}
		unlink $organism.'.fna';
		`tar -cvzf contigs.tar .`;
		opendir DIR, getcwd;
		@fna_files = readdir DIR;
		@fna = grep($_ =~ m/fna$/, @fna_files);
		foreach (@fna) {							# working with 1 single chromosome record
			my $match_count = 0;
			my $fna_file = $_;
			(my $gname) = $_ =~ m/(.*)\.fna/;			# name of the chromosome
			`transeq -sequence "$fna_file" -frame 6 -table 11 -clean -outseq "$start_dir/$out_dir/temp_orfs.txt" 2>/dev/null`;
			# print getcwd;
			# print "current file: $fna_file \n";
			# Warning: transeq switches the order of the -1, -2, and -3 frames unpredictably. There is no way to avoid having an error of +1/-1 in the 
			# precursor locations report for these strands
			my @temp = readFASTA($fna_file, 1); 
			my @prot = @{$temp[0]};
			my @seql = @{$temp[2]};
			my @seq = @{$temp[1]};
			my ($prec_ref, $prec_seq_ref, $groups_ref) = read_transeq("$start_dir/$out_dir/temp_orfs.txt", 20, 200, $seql[0]);	
			if ($prec_ref == 0) {next;}
			# print "400 \n";	
			my @prec = @$prec_ref; my @prec_seq = @$prec_seq_ref; my %groups = %$groups_ref;     # precursor peptides on that chromosome
			`rm "$start_dir/$out_dir/temp_orfs.txt"`;		# delete transeq output										
			my ($match_prot_ref, $match_seq_ref) = pattern_match(\@prec, \@prec_seq);
			my ($avg_ref, $match_clusters_ref); 
			($match_prot_ref, $match_seq_ref, $avg_ref) = cluster_prec($match_prot_ref, $match_seq_ref);
			my @match_prot = @$match_prot_ref; 			# precursor peptides on the chromosomes
			my @match_seq  = @$match_seq_ref;			# that match the pattern
			my @avg = @$avg_ref; # This list has midpoints of the clusters
			foreach (@match_prot) {
				my @to_count = @{$_};
				$match_count += $#to_count+1;
			}
			push @match_counts, $match_count;
			$tot_count = $#prec+1;
			my (@select_match_rank, @select_nucl_name, @select_match_seq, 
			    @select_match_B, @select_match_loc, @select_match_C, @select_match_orgn, %select_groups);
			our %AME_scores; # Hash of the scores of all maturation enzymes	
			our %AME_scores_2;	
			our %AME_scores_3;
			# Structure AME_scores{seq} = [rankB, Bloc, rankC, Cloc]			    
			my (@neig, @neig_seq); #orf_seqs is only temporary (see below)
			my %neig_binned; # bin for each searching
			
			if (@match_prot and $avg[0]) {	# edited				
				my (@totrank, @B, @C, @loc, @organism);	
									
				# Get all the neighbors at once
				print $fna_file."\n";
				`getorf -sequence "$fna_file" -minsize 200 -maxsize 3000 -find 1 -auto -outseq "$start_dir/$out_dir/neighbor_orfs.txt"`; # get all the neighbors				
				$num_getorf ++;
				my ($neig_ref, $neig_seq_ref) = readFASTA("$start_dir/$out_dir/neighbor_orfs.txt");
				@neig = @$neig_ref; @neig_seq = @$neig_seq_ref; # dereference	
				$tot_neig = $#neig+1;				
				`rm "$start_dir/$out_dir/neighbor_orfs.txt"`;		# delete getorf output	
			    					
				foreach (0..$#neig) {
					$neig[$_] =~ m/\[(.*) - (.*)\]/;					
					my $start = $1; my $stop = $2;
					my $direction = 'f';	
					if ($start > $stop) { 
						($start,$stop) = ($stop,$start);
						$direction = 'r';
					}
					push @{$neig_binned{floor($start/$bins)}}, ["$start\t$stop\t$direction", $neig_seq[$_]]; 	
				}
				
				my @nucl_name;
				foreach (0..$#avg) {
					my $start = $avg[$_]-250; my $stop = $avg[$_]+250;
					my $Aloc = '|';
					foreach (@{$match_prot[$_]}) {
						$Aloc .= $_.'|';
					}
					$nucl_name[$_] = "$gname";					
					my ($neig_ref, $neig_seq_ref) = get_neighbors(\%neig_binned, $start, $stop);
					my ($totrank_ref, $Brank_ref, $B_ref, $C_ref, $Bloc, $Cloc) = rank_hits($neig_ref, $neig_seq_ref); # 
					# if ($totrank_ref > 3) { #also test for B2, B1
						# $B1_ref = "none";
						# $B1loc = "none\tnone\tNA";
						# ($Brank2_ref, $B1_ref, $B2_ref, $B1loc, $B2loc) = rank_hits3($neig_ref, $neig_seq_ref);
						# if ($Brank2_ref > $Brank_ref) {
							# $B_ref = $B2_ref;
							# $Bloc = $B2loc;
							# $totrank_ref = ($totrank_ref - $Brank_ref) + $Brank2_ref;
							# }
						# }
					my ($D_ref, $E_ref, $Dloc, $Eloc, $PTM);
					if ($totrank_ref > 5) {
						#$PTM = rank_hits_PTM($neig_ref, $neig_seq_ref, $Bloc, $Cloc);
						($D_ref, $E_ref, $Dloc, $Eloc) = rank_hits2($neig_ref, $neig_seq_ref);
					}
					else {
						$Dloc = "none\tnone\tNA";
						$Eloc = "none\tnone\tNA";
						$B1loc = "none\tnone\tNA";
					}
					push @totrank, $totrank_ref;
					push @B, $B_ref;
					push @C, $C_ref;
					push @loc, $Aloc."\t".$B1loc."\t".$Bloc."\t".$Cloc."\t".$Dloc."\t".$Eloc;#"\t".$PTM;
					push @organism, $organism;
					$get_rank_loops ++;	
				}
				
				my @ranked = sort {$totrank[$b] <=> $totrank[$a]} 0..$#totrank;
				foreach (@ranked) {
					if ($totrank[$_] >= 1) {
						push @select_match_rank, $totrank[$_];
						push @select_nucl_name, $nucl_name[$_];
						push @select_match_seq,  $match_seq[$_];
						push @select_match_B, $B[$_];
						push @select_match_C, $C[$_];
					    push @select_match_loc, $loc[$_];
						push @select_match_orgn, $organism[$_];
						foreach (@{$match_seq[$_]}) {
							$select_groups{$_} = $groups{$_};
						}
						}
				}
			}
			
			push @g_match_rank, @select_match_rank;		# add all the hits from each fna to the
			push @g_nucl_name, @select_nucl_name;		# overall hit list
			push @g_match_seq,  @select_match_seq;
			push @g_match_B, @select_match_B;
			push @g_match_C, @select_match_C;
			push @g_match_loc, @select_match_loc;
			push @g_match_orgn, @select_match_orgn;
			push @neig_counts, $tot_neig;
			push @tot_counts, $tot_count;
			push @organisms, $organism;
			$g_tot_count = $g_tot_count + $tot_count;
			$g_tot_neig = $g_tot_neig + $tot_neig;
			%g_groups = (%g_groups, %select_groups);
			undef %AME_scores;
			undef %AME_scores_2;
			undef %AME_scores_3;
			unlink $fna_file;
		}
			
		foreach (@match_counts) {
			$g_match_count += $_;
		}
		$end_time = time(); # this gets updated after a single genome is done
		chdir $maindir;
		$progress ++;
		print "".($nstart+$progress)." of ".($nstop+1)." genomes done in ";
		printf("%.2f", $end_time - $start_time);
		print " seconds.   Last genome was of ".$organism." \n";
		open LOG, ">>$start_dir/$out_dir/run_log.txt";
		print LOG "".($nstart+$progress)." of ".($nstop+1)." genomes done in ";
		printf LOG ("%.2f", $end_time - $start_time);
		print LOG " seconds.   Last genome was of ".$organism." \n";
		close LOG;
		open F1, ">>$start_dir/$out_dir/precursors.txt";
        open F2, ">>$start_dir/$out_dir/machinery.txt";
#        open F4, ">>$start_dir/$out_dir/clusters.txt";
        open F5, ">>$start_dir/$out_dir/prec_locs.txt";
#        print F1 "Rank\tNucl ID\tNumber of precursors\tBloc\tBdir\tCloc\tCdir\tPrecursor sequence\tOrganism\n";
        foreach (0..$#g_match_rank) {
          my $Adir;
          my ($Aloc, $B1loc1, $B1loc2, $B1dir, $Bloc1, $Bloc2, $Bdir, $Cloc1, $Cloc2, $Cdir, $Dloc1, $Dloc2, $Ddir, $Eloc1, $Eloc2, $Edir) = split ("\t", $g_match_loc[$_]);
		  my @temp_prec_list = @{$g_match_seq[$_]};
 	      print F1 $g_match_rank[$_]."\t".$g_nucl_name[$_].
 	      "\t".($#temp_prec_list+1)."\t".$B1loc1."\t".$B1loc2."\t".$B1dir."\t".$Bloc1."\t".$Bloc2."\t".$Bdir."\t".$Cloc1."\t".$Cloc2."\t".$Cdir."\t".$Dloc1."\t".$Dloc2."\t".$Ddir."\t".$Eloc1."\t".$Eloc2."\t".$Edir.
 	      "\t".$temp_prec_list[0]."\t".$g_match_orgn[$_]."\n";
 	      print F2 $g_match_B[$_]."\t".$g_match_C[$_]."\n";
	      print F5 "$Aloc\n";
 	      foreach (1..$#temp_prec_list) {
 	      		print F1 "> $temp_prec_list[$_]\t\n";
 	      		print F2 "prev\n";
 	      		print F5 "prev\n";
 	      }
#	      foreach (@temp_prec_list) {
#		      my @groups_seqs = @{$g_groups{$_}};	      
#	      	  foreach (reverse 0..$#groups_seqs) {
#	      	  	my $seq = $groups_seqs[$_];
#	      	  	print F4 "$seq\t";
#	          }
#	          print F4 "\n";		                                     	          
#	      }
        }
        open F3, ">>$start_dir/$out_dir/statistics.txt";
        foreach (0..$#match_counts) {
	    	print F3 $match_counts[$_]."\t".$tot_counts[$_]."\t".$neig_counts[$_]."\t".$organisms[$_]."\n"
        }
        close F1;
        close F2;
        close F3;
#        close F4;
        close F5;
	}	
	chdir $start_dir;								# go back to the original directory
	open F3, ">>$start_dir/$out_dir/statistics.txt";
	print F3 "Precursor pattern that was analyzed: $pattern\n";
	print F3 "\nTotal number of matches: ".$g_match_count."\n";
    print F3 "Total number of precursors: ".$g_tot_count."\n";
    print F3 "Total number of neighbors found: ".$g_tot_neig."\n";
    print F3 "Total number of neighbors analyzed (MAST calls): ".$num_mast."\n";    
    print F3 "getorf has been called $num_getorf times\n";
    print F3 "Total number of loop 'get_neighbors' and 'rank hits' was called: $get_rank_loops\n";
    close F3;
}

print "done";
