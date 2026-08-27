#!/usr/bin/env perl

#AUTHORS
#   Rene Warren
#   Lauren Coombe

#NAME
#   ntRootAncestryPredictor.pl

#SYNOPSIS
#   ntRoot : ntedit-powered human super-population-level ancestry predictions using 1000 Genomes Project integrated variant call set
#
#   This version additionally supports RECOMBINATION-BASED tiling.

#DOCUMENTATION
#   Readme distributed with this software @ www.bcgsc.ca
#   http://www.bcgsc.ca/platform/bioinfo/software/ntroot
#   http://www.bcgsc.ca/platform/bioinfo/software/ntedit
#   We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca
#   If you use ntRoot, ntEdit, the ntEdit code or ideas, please cite our work

#LICENSE
#   ntRoot Copyright (c) 2024-now British Columbia Cancer Agency Branch.  All rights reserved.
#   ntRoot and companion code is released under the GNU General Public License v3

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_t $opt_v $opt_r $opt_i $opt_b);

my $dw = 5000000;
my $verbose = 0;
my $tile_resolution = 0;
my $fai = "";
my $recomb_bed = "";                 ### NEW: recombination-tile BED file
my $ambiguous_category = "Unknown";

getopts('f:t:v:r:i:b:');             ### NEW: added b:

sub usage_page {
	print "\nUsage: $0 -f *variants.vcf [-t TILE_SIZE] [-v VERBOSITY] [-r TILE_OUTPUT] [-i FAI] [-b RECOMB_BED]\n";
	print "\t-f\tVariants VCF file\n";
	print "\t-t\tTile size [$dw bp]\n";
	print "\t-v\tVerbose mode - 0 (False) or 1 (True) [0]\n";
	print "\t-r\tOutput ancestry inferences per tile - 0 (False) or 1 (True) [0]\n";
	print "\t-i\tReference FAI file (Only required when -r specified)\n";
	print "\t-b\tRecombination-tile BED file (chrom start end name [cM]).\n";
	print "\t  \tWhen supplied, ntRoot runs in RECOMBINATION mode and writes\n";
	print "\t  \t*_ancestry-predictions-recomb-resolution.tsv (requires -i FAI).\n\n";
}


if (!$opt_f || ($opt_r && !$opt_i) || ($opt_b && !$opt_i)) {   ### NEW: -b
	usage_page();
	exit(1);
}

my $f = $opt_f;
$dw = $opt_t if ($opt_t);
$verbose = $opt_v if ($opt_v);
$tile_resolution = $opt_r if ($opt_r);
$fai = $opt_i if ($opt_i);
$recomb_bed = $opt_b if ($opt_b);    ### NEW


my $chr;
# Read in the FAI file for chromosome lengths 
if ($tile_resolution || $recomb_bed) {
	open(IN,$fai) || die "Can't read $fai --fatal (is the file in your working directory?)\n";
	while(<IN>){
	   chomp;
	   my @a=split(/\t/);
	   (my $nc = $a[0]) =~ s/^chr//i;         ### NEW: store under normalized key too
	   $chr->{$a[0]}=$a[1];
	   $chr->{$nc}=$a[1];
	}
	close IN; ### filehandle hygiene
}

#############################################################
# NEW: read recombination tiles from BED into per-chr arrays
#############################################################
my %rc_start;   # normalized_chr => [start, start, ...]  (ascending)
my %rc_end;     # normalized_chr => [end, end, ...]
my %rc_name;    # normalized_chr => [name, ...]
my %rc_cm;      # normalized_chr => [cM_width or "", ...]

if ($recomb_bed) {
	open(RB, $recomb_bed) || die "Can't read recombination BED $recomb_bed --fatal.\n";
	my %tmp;   # nchr => [ [start,end,name,cm], ... ]
	while(<RB>){
		chomp;
		next if (/^#/ || /^\s*$/ || /^track/ || /^browser/);
		my @a = split(/\t/);
		next if (@a < 3);
		(my $nc = $a[0]) =~ s/^chr//i;          # normalize: chr21 -> 21
		my $st  = $a[1] + 0;
		my $en  = $a[2] + 0;
		my $nm  = defined $a[3] ? $a[3] : "$nc:$st-$en";
		my $cm  = defined $a[4] ? $a[4] : "";
		push @{$tmp{$nc}}, [$st, $en, $nm, $cm];
	}
	close RB;
	# sort each chromosome's tiles by start, populate parallel arrays
	foreach my $nc (keys %tmp){
		my @sorted = sort { $a->[0] <=> $b->[0] } @{$tmp{$nc}};
		foreach my $t (@sorted){
			push @{$rc_start{$nc}}, $t->[0];
			push @{$rc_end{$nc}},   $t->[1];
			push @{$rc_name{$nc}},  $t->[2];
			push @{$rc_cm{$nc}},    $t->[3];
		}
	}
	my $ntiles = 0; $ntiles += scalar(@{$rc_start{$_}}) for keys %rc_start;
	print "Loaded $ntiles recombination tiles across ".scalar(keys %rc_start)." chromosomes from $recomb_bed\n";
}

# Binary search: rightmost tile whose start <= pos (clamps to first/last tile)
sub find_recomb_tile {
	my ($c, $pos) = @_;
	my $starts = $rc_start{$c};
	return -1 unless defined $starts;
	my $hi = $#{$starts};
	return -1 if $hi < 0;
	return 0 if $pos < $starts->[0];        # clamp below to first tile
	my ($lo, $ans) = (0, 0);
	while ($lo <= $hi) {
		my $mid = int(($lo + $hi) / 2);
		if ($starts->[$mid] <= $pos) { $ans = $mid; $lo = $mid + 1; }
		else { $hi = $mid - 1; }
	}
	return $ans;                            # clamps above to last tile automatically
}

###below support for compressed vcf
my $IN;

if ($f =~ /\.(gz|bgz)$/) {### add support for gz/bgz
    my $cmd = "gunzip -c";  # safe default
    $cmd = "unpigz -c" if system("command -v unpigz >/dev/null 2>&1") == 0;

    open($IN, "-|", split(" ", $cmd), $f) or die "can't read compressed $f -- fatal.\n";###safer open
} else {
    open($IN, "<", $f) or die "can't read $f -- fatal.\n";
}



my $xr=0;
my $s;
my $y;      # fixed-window totals:  $y->{chr}{wn}{'ct'}
my $z;      # fixed-window per-pop: $z->{chr}{wn}{pop}{'sum'|'nzct'}
my $yr;     # NEW recomb totals:    $yr->{nchr}{rn}{'ct'}
my $zr;     # NEW recomb per-pop:   $zr->{nchr}{rn}{pop}{'sum'|'nzct'}
my $populations;

print "Inferring ancestry using SNVs...\n";

while(<$IN>){
	chomp;
	my @a=split(/\t/);
	my $max = -1;###never previously initialized
	my $maxpop;

	if(/_AF/){
		my $wn = int($a[1] / $dw);                       # fixed physical tile index
		(my $nchr = $a[0]) =~ s/^chr//i;                 # NEW: normalized chr
		my $rn = $recomb_bed ? find_recomb_tile($nchr, $a[1]) : -1;   # NEW: recomb tile index
		$xr++;

                my @alleles = split(/\^/,$a[7]);

                #21      5097811 .       G       A       .       PASS    AD=11^AC=115;AN=5096;DP=12758;AF=0.02;EAS_AF=0.03;EUR_AF=0.01;AFR_AF=0.03;AMR_AF=0.03;SAS_AF=0.03;VT=SNP;NS=2548        GT      1/1

                foreach my $allele(@alleles){

			my @b=split(/\;/,$allele);

			foreach my $el(@b){
				my @d=split(/\=/,$el);
				if($d[0]=~/(\S+)\_AF/){
					my $pop=$1;
					my @AFalleles = split(/,/,$d[1]);
					foreach my $afallele(@AFalleles){
						if ($afallele !~ /\d+/) {
							next;
						}
						if (! defined $populations->{$pop}) {
							$populations->{$pop} = 1;
						}
						$s->{$d[0]}{'sum'}+=$afallele;

						#chr  winnum   pop  (fixed physical tiles)
						$z->{$a[0]}{$wn}{$pop}{'sum'}+=$afallele;
						$y->{$a[0]}{$wn}{'ct'}++;

						# NEW: recombination tiles (parallel accumulation)
						if ($recomb_bed && $rn >= 0) {
							$zr->{$nchr}{$rn}{$pop}{'sum'} += $afallele;
							$yr->{$nchr}{$rn}{'ct'}++;
						}

						if($afallele){
							$s->{$d[0]}{'ct'}++;
							$z->{$a[0]}{$wn}{$pop}{'nzct'}++;
							$zr->{$nchr}{$rn}{$pop}{'nzct'}++ if ($recomb_bed && $rn >= 0);  # NEW
							if($a[1]>$max){
								$max=$a[1];
								$maxpop=$pop;
							}
						}
					}
				}
			}
		}
   		}
	}

close $IN;

###calculate metric per tile
my $top;
my $total;
my @ordered_populations = sort keys %$populations;

if(! $xr){
	print "\n! There are no cross-referenced SNV in $f; no ancestry predictions can be reported !\n\n\tDid you:\n\t1) Run ntedit with the correct and properly-formatted human genome input\n\t\te.g., chromosome 14 should be: >14\n\t\te.g., -f GRCh38.fa\n\n\t2) Supply the 1000 Genomes Project integrated variant callset vcf to ntedit with -l\n\n";
	exit(1);
}

#####################################################################
# NEW: RECOMBINATION MODE (recomb-only) -- write recomb TSV and exit
#####################################################################
if ($recomb_bed) {
	my $rout = $f . "_ancestry-predictions-recomb-resolution.tsv";
	open(RCO, ">$rout") || die "Can't write to $rout -- fatal.\n";
	print RCO "chrom\tstart\tend\tancestry_prediction";
	foreach my $population (@ordered_populations) {
		print RCO "\t$population-score";
	}
	print RCO "\ttile_id\tcM_width\n";

	my $rtop;      # recomb-weighted global fraction (bp)
	my $rtotal = 0;

	foreach my $el (sort keys %$zr){
		my $tl = $zr->{$el};
		foreach my $rn (sort {$a<=>$b} keys %$tl){
			my $pl = $tl->{$rn};
			my $winmax = -1;
			my $winpop;
			my $window_population_metric;
			my $ct = $yr->{$el}{$rn}{'ct'} || 0;
			my $start = $rc_start{$el}[$rn];
			my $end   = $rc_end{$el}[$rn];
			my $name  = $rc_name{$el}[$rn];
			my $cm    = $rc_cm{$el}[$rn];
			$cm = "NA" if (!defined $cm || $cm eq "");
			print "WARNING: chr$el $name has $ct only total SNVs -- recombination tile may be too small for a confident call\n" if($ct < 100);
			next if $ct == 0;
			foreach my $pp (keys %$pl){
				my $rate   = $pl->{$pp}{'sum'} / $ct;
				my $nz     = $pl->{$pp}{'nzct'} || 0;
				my $metric = ($nz / $ct) * $rate;
				$window_population_metric->{$pp} = $metric;
				if($metric > $winmax){ $winmax = $metric; $winpop = $pp; }
				elsif($metric == $winmax){ $winmax = $metric; $winpop = $ambiguous_category; }
			}
			my $width = ($end - $start + 1);
			$rtop->{$winpop} += $width;
			$rtotal += $width;
			# 1-based inclusive start for output, consistent with physical-tile output
			my $out_start = $start + 1;
			print RCO "$el\t$out_start\t$end\t$winpop";
			foreach my $population (@ordered_populations) {
				printf RCO "\t%.4f", ($window_population_metric->{$population} || 0);
			}
			print RCO "\t$name\t$cm\n";
		}
	}
	close RCO;

	# Print a recombination-weighted GLOBAL ancestry summary to STDOUT (cross-check)
	print "\n# Recombination-weighted global ancestry (bp fraction across recomb tiles):\n";
	if ($rtotal > 0) {
		foreach my $population (sort {$rtop->{$b}<=>$rtop->{$a}} keys %$rtop){
			printf "#   %-8s %.2f%%\n", $population, ($rtop->{$population}/$rtotal*100);
		}
	}
	print "\nRecombination-resolution ancestry predictions available in:\n$rout\n\n";
	exit(0);
}
#####################################################################
# END recombination mode
#####################################################################

if ($tile_resolution) {
	my $best = $f . "_ancestry-predictions-tile-resolution_tile$dw.tsv";
	open(BEST,">$best") || die "Can't write to $best -- fatal.\n";
	print BEST "chrom\tstart\tend\tancestry_prediction";
	foreach my $population (@ordered_populations) {
		print BEST "\t$population-score";
	}
	print BEST "\n";
}

foreach my $el(sort {$a<=>$b} keys %$z){
	my $wnl=$z->{$el};
	foreach my $wnum(sort {$a<=>$b} keys %$wnl){
		my $pl = $wnl->{$wnum};
		my $winmax = -1;###was never initialized
		my $winpop;
		my $window_population_metric;
		print "WARNING: chr$el tile$wnum has $y->{$el}{$wnum}{'ct'} only total SNVs -- you may need to increase the tile size (currently set at $dw)\n" if($y->{$el}{$wnum}{'ct'}<100);
		foreach my $pp(keys %$pl){
			my $ct = $y->{$el}{$wnum}{'ct'} || 0;###guards div by zero
                        next if $ct == 0;                    ###guards div by zero
			my $rate = $pl->{$pp}{'sum'}/$y->{$el}{$wnum}{'ct'};
			my $metric = ($pl->{$pp}{'nzct'}/$y->{$el}{$wnum}{'ct'}) * $rate;
			$window_population_metric->{$pp} = $metric;
			if($metric>$winmax){
				$winmax = $metric;
				$winpop = $pp;
			}elsif($metric==$winmax){ ###do not assign tiles to a population if ambiguous
				$winmax = $metric;
				$winpop = $ambiguous_category;
			}
		}
		$top->{$winpop} += $dw;
		$total += $dw;
		if ($tile_resolution) {
			my $chunk = $wnum * $dw;
			my $start = $chunk + 1;
			my $end = $chunk + $dw;
			$end = $chr->{$el} if($end>$chr->{$el});
			print BEST "$el\t$start\t$end\t$winpop";
			foreach my $population (@ordered_populations) {
				printf BEST "\t%.4f", ($window_population_metric->{$population});
			}
			print BEST "\n";
		}

	}
}

if ($tile_resolution) {
	close BEST;
}


#calculate/incorporate metrics
foreach my $k(keys %$s){
	my $p = $s->{$k}{'sum'}/$xr;
	my $c = $s->{$k}{'sum'}/$s->{$k}{'ct'};
	my $nzr = $s->{$k}{'ct'}/$xr;
	$s->{$k}{'prob'} = $p*$nzr;
	$s->{$k}{'fract'} = $p*$s->{$k}{'ct'};
}

#output predictions

my $out = $f . "_ancestry-predictions_tile$dw.tsv";
open(OUT,">$out") || die "Can't write to $out -- fatal.\n";


my $header_str = "# GAI score: Average SNV allele frequency * rate of SNVs with non-zero allele frequency\n";
$header_str .= "# Populations ranked by LAI fraction\n";
$header_str .= "# AF: Allele Frequency; nz: Non-zero\n";
$header_str .= "GAI Super-population\tLAI fraction (tile:$dw bp)\tGAI score\tTotal SNV count\tNon-zero AF SNV count";

if ($verbose) {
	$header_str = $header_str . "\tSumAF\tAvgAF\tnzAvgAF\tnzSNVrate\tAvgAF * nzAF_SNV_count\n";
} else {
	$header_str = $header_str . "\n";
}

print OUT $header_str;

my $rank=0;
foreach my $population(sort {$top->{$b}<=>$top->{$a}} keys %$top){
	$rank++;
	my $k = $population . "_AF";
	my $percent = $top->{$population}/$total *100;
	if ($population eq $ambiguous_category) {
		printf OUT "$population\t%.2f%%\tN/A\t$xr\tN/A", $percent;
	} else {
		printf OUT "$population\t%.2f%%\t%.4f\t$xr\t$s->{$k}{'ct'}", ($percent, $s->{$k}{'prob'});
	}
	if ($verbose) {
		if ($population eq $ambiguous_category) {
			printf OUT "\tN/A\tN/A\tN/A\tN/A\tN/A\n";
		} else {
			my $p = $s->{$k}{'sum'}/$xr;
			my $c=$s->{$k}{'sum'}/$s->{$k}{'ct'};
			my $nzr=$s->{$k}{'ct'}/$xr;
			printf OUT "\t%.2f\t%.4f\t%.4f\t%.4f\t%.2f\n", ($s->{$k}{'sum'}, $p, $c, $nzr, $s->{$k}{'fract'});
		}
	} else {
		printf OUT "\n";
	}
}

close OUT;### filehandle hygiene

print "Ancestry predictions available in:\n$out\n";

if ($tile_resolution) {
	print $f . "_ancestry-predictions-tile-resolution_tile$dw.tsv\n";
}
print "\n";

exit(0);
