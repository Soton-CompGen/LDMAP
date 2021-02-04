use warnings; use strict;

## Script to convert LDmap output  to CSV
## RJPengelly, 2014, r.j.pengelly@soton.ac.uk

open FILE, $ARGV[0];
open OUT, ">$ARGV[0]" . '.csv';

print OUT 'N,Locus,KBmap,LDUmap';
while (<FILE>) {
	if ($_ =~ /#/) {next;}
	$_ =~ s/^\s+//; 
	$_ =~ s/\s+$//;
	$_ =~ s/\s+/\,/g; 
	print OUT "$_\n";
}
