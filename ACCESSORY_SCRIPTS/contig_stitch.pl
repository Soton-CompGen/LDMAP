#!perl;
use strict;
use warnings;

## Takes file list of split LDMAP outputs, merges contigs together, averaging along the overlapping segment. Must use CSV input.
## Calling input "perl contig_stitch.pl {list file} {overlap} {trim}"



## Read input and validate, prepare files.
my $overlap_max = $ARGV[1]; my $trim = $ARGV[2];
if (((defined $overlap_max) and ($overlap_max >= 0)) and ((defined $trim) and ($trim > 0))) {} else {die "\nEnsure all required variables are supplied on calling.\nSource file, overlap and desired trim distance must be supplied, in that order.\n\n";}
if (($trim * 2) > $overlap_max) {
	die "\nTrim length is greater than the length of overlap_max, would create holes in map!\nTrim is removed from both contig ends on joining.\nCheck input variables.\n\n";
}
open LIST, $ARGV[0] or die "\nUnable to open file containing list.\n";
chomp (my @file_list = <LIST>);
close LIST;
	print "\nRunning on: @file_list\n\n";
print "Settings are: overlap max = $overlap_max\; trim from each touching end = $trim.\n\n";

open LOG, '>concat_' . "$ARGV[0]" . '.log';
print LOG "\nRunning on: @file_list\n";
print LOG "Settings are: overlap max = $overlap_max\; trim from each touching end = $trim.\n\n";

open OUT, '>concat_' . "$ARGV[0]" . '.ter.csv';
close OUT;


## Loop through the file list and concatenate.
my $file1_name; my $file2_name; my $iter = 0;
while (@file_list) {
	## Declare files to be opened for merging
	if ($iter == 0) {
		$file1_name = (shift @file_list);
	}
	$file2_name = (shift @file_list);
	
	## Open first 2 files and discard headers
	open FILE1, $file1_name or die "\n cannot open $file1_name.";
	chomp (my @file1 = <FILE1>);
	close FILE1;
	if ($file1[0] !~ /,/) {die "Input data must be in CSV format\n";}
	if ($file1[0] =~ /^N/) {shift @file1;}
	
	open FILE2, $file2_name or die "\n cannot open $file2_name.";
	chomp (my @file2 = <FILE2>);
	close FILE2;
	shift @file2;
	if ($file2[0] =~ /^N/) {shift @file2;}

	## Count number of shared entries in both arrays to calculate the overlap.
	## Hash based will be quicker, may need to implement in the future.
	my $cnt = 0;
	
	## Hash try
	my %file2_temp;
	foreach (@file2) {
		my @split = split /,/, $_;
		$file2_temp{$split[1]} = $_;
	}
	foreach my $line_1 (@file1) {
		my @split_line_1 = split /,/, $line_1;
		if ( exists $file2_temp{$split_line_1[1]}  ) {
			$cnt ++;
		}
	}
	undef %file2_temp; 
		
	my $overlap = $cnt;
	print LOG "$file1_name\t$file2_name\t$overlap entry overlap\n";
	if ($overlap > $overlap_max) {die "Overlap has exceeded defined max of $overlap_max between $file1_name & $file2_name.\n\n";};

	## Define LDU prior to trim point in file 2
	my @temp = split /,/, $file2[($trim-1)]; 
	my $file2_start_tide = $temp[3];
	undef @temp;
	
	## Trim desired length off arrays
	if ($trim > 0) {
		@file1 = reverse @file1;
		splice @file1, 0, $trim;
		@file1 = reverse @file1;
		splice @file2, 0, $trim; ## Removes from start, counting forwards.
	}
	
	## Count remaining overlapping entries after trim
	$cnt = 0;
	foreach (@file2) {
		my @split = split /,/, $_;
		$file2_temp{$split[1]} = $_;
	}
	foreach my $line_1 (@file1) {
		my @split_line_1 = split /,/, $line_1;
		if ( exists $file2_temp{$split_line_1[1]}  ) {
			$cnt ++;
		}
	}
	undef %file2_temp;
	$overlap = $cnt;
	
	@file1 = reverse @file1;
	my @over1 = splice @file1, 0, $overlap;
	@over1 = reverse @over1;
	@file1 = reverse @file1;
	my @over2 = splice @file2, 0, $overlap; ## Removes from start, counting forwards.
	
	## Define value for end of overlap for subtraction later to narmalise with the meaned value
	@temp = split /,/, $over2[-1];
	my $over2_end_tide = $temp[3];
	undef @temp;
	
	## Define LDUvalue for starting meaned contig
	@temp = split /,/, $file1[-1];
	my $tide1 = $temp[3];
	undef @temp;
	
	
	## Average values within overlap
	## Unpaired values will be added to array too, without alteration. Warning will be printed.
	my @meaned;
	foreach my $line1 (@over1) {
		my @loc1 = split /,/, $line1;
		my $index = 0; my $found = 0;
		foreach my $line2 (@over2) {
			my @loc2 = split /,/, $line2;
			if ($loc1[1] eq $loc2[1]) {
				$found++;
				$loc1[3] = (($loc1[3] + ($loc2[3] + $tide1 - $file2_start_tide)) / 2);
				$line1 = join ',', @loc1,'aved';
				push @meaned, $line1;
				splice @over2, $index, 1;
				last;
			}
			$index++;
		}
		if ($found == 0) { 
			print STDERR "No pair found for $line1\n";
			push @meaned, $line1;
		}
	}
	foreach my $line2 (@over2) {
		#print STDERR "No pair found for $line2\n";
		my @loc2 = split /,/, $line2;
		$loc2[3] = ($loc2[3] + $tide1 - $file2_start_tide);
		$line2 = join ',', @loc2;
		push @meaned, $line2;
	}
	@temp = split /,/, $meaned[-1];
	my $tide2 = $temp[3];
	undef @temp;
	
	## Print values to final file
	open OUT, '>concat_' . "$ARGV[0]" . '.ter.csv';
	print OUT "N,Locus,KBmap,LDUmap,aved?\n";
	
		foreach (@file1) {
			print OUT "$_\n";
		}
	foreach (@meaned) {
		print OUT "$_\n";
	}
	
	
	foreach (@file2) {
		my @loc = split /,/, $_;
		$loc[3] += ($tide2 - $over2_end_tide);
		my $str = join ',', @loc;
		print OUT "$str\n";
	}
	close OUT;
	
	## Define variables for next loop
	$file1_name = 'concat_' . "$ARGV[0]" . '.ter.csv';	
	$iter ++;
}

## Rejoice in final success!
print "File merging complete!\n";

__END__
Takes file list of split LDMAP outputs, merges contigs together, averaging along the overlapping segment. Must use CSV input.
Input file list must be in correct order, errors in this should be fairly easy to spot on a discontinuous line on the resultant plot
Input must have File, overlap of existing contigs, and desired trim distance.
Trimming can be used as the quality of mapping will be greatly reduced at the end of segments, should also work to assemble the non-overlapping 'contig-free' approach, see WLau thesis.
Trim is removed from overlap region!! This does require some pre-emptive mental arithmetic.

Awkward to merge the two contigs as the overlap region that was defined for the VCF splitting is not necessarily the same as that which remains in the LDMAP output.
LDMAP software does further QC check on loci, which can reduce the number of loci from those entered.
