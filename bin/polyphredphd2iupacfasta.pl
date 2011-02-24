#!/usr/bin/perl -w

# for contig level 
#my $poly_out;
#my $read_fasta;
#my $read_fasta_out;
#my $contig_fasta;
#my $contig_fasta_out;
#my $ace;

my $WARNING = 1;

#otw ok with phd file alone
my $phd_file =$ARGV[0];

#my $read_fasta_out;
#($read_fasta_out) = ($phd_file =~ m/(.+?)\.phd\.\d+/);

# read sequence

my @seq;
my $in_dna=0;
my $in_tag=0;
my $in_tag_comment=0;
my $type ="";
my $tag_pos=-1;
my $tag_qual = 0;
my $name = $phd_file;

open PHD, "<$phd_file";

while (my $r = <PHD>) {
    chomp $r;

    if ($r=~m/BEGIN_SEQUENCE\s+(.+)/) {
	$name=$1;
    }

    if ( $in_dna == 1 ) {
	if ($r =~ m/^END_DNA/) {
	    $in_dna = 0;
	} else {
	    
	    my @r= split(/\s+/,$r);
	    my $base = $r[0];
	    if ($r[1] > 25) {
		$base = uc $base;
	    }
	    push @seq, $base;
	}
    }

    if ($in_tag == 1) {
	if ($in_tag_comment == 1) {
	    if ($r=~m/^END_COMMENT/) {
		$in_tag_comment = 0;
	    } elsif ($r =~ m/^(\d+)/) {
		$tag_qual = $1;
	    }
	}

	if ($r=~m/^END_TAG/) {
	    $in_tag=0;

	    if($tag_pos != -1 and $type ne "") {
		if ($tag_qual == 99) {
		    $type = uc $type;
		}
		$seq[$tag_pos-1]=$type; 		# tag pos likely 1-based (looks like it from dataNeeded tags at least.. =)
	    }

	} elsif ($r=~m/^TYPE:\s+heterozygote(\w+)/) {
	    $type = $1;
	    if ($type eq "AG" or $type eq "GA") {
		$type = "r";
	    } elsif ($type eq "CT" or $type eq "TC") {
		$type = "y";
	    } elsif ($type eq "CG" or $type eq "GC") {
		$type = "s";
	    } elsif ($type eq "AT" or $type eq "TA") {
		$type = "w";
	    } elsif ($type eq "GT" or $type eq "TG") {
		$type = "k";
	    } elsif ($type eq "AC" or $type eq "CA") {
		$type = "m";
	    }elsif ($type eq "Indel") {
		$WARNING && print STDERR "WARNING: Indel in $name. Not marked on iupacfasta.\n";
		$type = "";
	    } else {

		$WARNING && print STDERR "WARNING: unknown heterozygote type $type detected for $name. \n";
		$type ="x";
	    }
	} elsif ($r=~m/^UNPADDED_READ_POS:\s+(\d+)/) {
	    $tag_pos =$1;
	} elsif ($r=~m/^BEGIN_COMMENT/) {
	    $in_tag_comment =1
	}
    }

    if ($r=~m/^BEGIN_DNA/) {
	$in_dna =1;
    }

    if ($r=~m/^BEGIN_TAG/) {
	$in_tag =1;
	$tag_pos = -1;
	$type = "";
    }

}

print ">$name\n";

while ( my @seqpart = splice @seq, 0, 80 ) {
    if (@seqpart == 0 ) {
	last;
	# shouldn't happen..
    }
    print join("", @seqpart), "\n"
} 
