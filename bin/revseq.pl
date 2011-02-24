#!/usr/bin/perl -w

my $fastafile = $ARGV[0];

open FASTA, "<$fastafile";

my $currentseq = "";
my %seq;
$first=1;

while(my $r = <FASTA>) {
    chomp $r;
    
    if($r=~/\>(.+)/) {
        if(!$first) {
            $seq{$currenthead} = $currentseq;
        } 
        $currenthead=$1;
        $currentseq ="";
        $first=0;
    } else {
        $currentseq .= $r;
    }           
}

$seq{$currenthead} = $currentseq;

foreach $name (keys %seq) {
    my $revseq = iupac_revcomp_polyphred($seq{$name});
    my @seq = split / */,$revseq; 

    print ">$name\n";

    while ( my @seqpart = splice @seq, 0, 80 ) {
	if (@seqpart == 0 ) {
	    last;
	    # shouldn't happen..
	}
	print join("", @seqpart), "\n";
    }
}

sub iupac_revcomp {
    $_ = shift;

    tr/ATUGCRYMKBDVHatugcrymkbdvh/TAACGYRKMACTGtaacgyrkmactg/;
    return scalar(reverse($_));
}

# polyphred iupac codes leaves the revcomped iupac code in place already. just leave them in, and revcomp the standard nucleotides.
sub iupac_revcomp_polyphred {
    $_ = shift;

    tr/ATUGCRYMKBDVHatugcrymkbdvh/TAACGRYMKBDVHtaacgrymkbdvh/;
    return scalar(reverse($_));
}
