#!/usr/bin/perl -w

=head1 NAME 

checkpoly.pl - Scan a polyphred polyout file and construct a report based on the results

=head1 SYNOPSIS

 We assume you have already run

 ~/src/acrscreen/bin/run_screen.sh

 Then simply    

 ~/src/acrscreen/bin/checkpoly.pl < edit_dir/hc_mptl_screen.ref.fasta.consedified.polyout

 or

 ~/src/acrscreen/bin/checkpoly.pl --polyout edit_dir/hc_mptl_screen.ref.fasta.consedified.polyout

=head1 OPTIONS AND ARGUMENTS

=over 8

=item C<--polyout>

Polyphred output file to analyse. Can also be given on STDIN to the program.

=item C<--reads>

Basecalled reads in a multifasta file. Same as went into the assembly that generated the polyout. Unused at this point.

=item C<--reference>

Reference fasta file.

=item C<--outfile>

Analysis output file.

=back

=head1 LICENSE AND COPYRIGHT

Copyright held by Daniel Nilsson, 2010.

The package is released under the Perl Artistic License.

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@izb.unibe.ch, daniel.nilsson@ki.se.

=cut

use Getopt::Long;
use Pod::Usage;

my $polyout_file = "";
my $read_fasta_file = "";
my $reference_fasta_file = "";
my $output_file_name="";

my $DEBUG = 0;
my $help = 0;

GetOptions( 'reads=s' => \$read_fasta_file,
	    'polyout=s' => \$polyout_file,
	    'debug' => \$DEBUG,
	    'reference=s', \$reference_fasta_file,
	    'outfile=s' => \$output_file_name,
	    'help|?' => \$help ) || pod2usage(2);

$help && pod2usage(1);

my $polyoutfh = *STDIN;
if($polyout_file ne "") {
    local *POLYIN;
    open POLYIN, "<$polyout_file";
    $polyoutfh = *POLYIN;
} else {
    $DEBUG && print "Reading polyout from stdin.\n";
}

my $outputfh = *STDOUT;
if($output_file_name ne "") {
    local *OUTPUT;
    open OUTPUT, ">$output_file_name";
    $outputfh = *OUTPUT;
} else {
    $DEBUG && print "Writing to stdout.\n";
}

use Bio::Seq;
use Bio::SeqIO;

my $refseq;
if($reference_fasta_file ne "") {
    
    my $in  = Bio::SeqIO->new(-file => "$reference_fasta_file" ,
			   -format => 'Fasta');

    $refseq = $in->next_seq();
}

my %reads;
if($read_fasta_file ne "") {
    
    my $in  = Bio::SeqIO->new(-file => "$read_fasta_file" ,
			   -format => 'Fasta');

    while ( my $read = $in->next_seq() ) {
	$reads{$read->display_id} = \$read;
    }
}

#global parser state
my $in_section=0;
my $section = "";
my @section_stack=();

# poly result vars
my @mod;
my @premod;

# poly analysis input vars
my $presumed_coding_start = 348;
my $presumed_coding_end = 2037;

# polyindel result vars
my $indelpolys = 0;

# genotype deviation result hash {worm}{refpos}
my %variant;

while(my $r=<$polyoutfh>) {
    chomp $r;
    if($r=~m/^BEGIN_(\w+)/) {
	# save last session, if we had one

	if ( $section ne "" ) {
	    push @section_stack, $section;
	}
	$section = $1;
	$in_section = 1;
	# call section open

	if ($section eq "POLY") {
	    @mod = (0,0,0);
	    @premod = (0,0,0);
	    print $outputfh "Displaying entries for non-silent polymorphisms.\n";
	} elsif ($section eq "GENOTYPE") {
	    $synonymous=0;
	    $nonsynonymous=0;
	    $heterozygote_positions=0;
	}
    }
    elsif($r=~m/^END_(\w+)/) {
	my $this_section = $1;
       
	if(scalar(@section_stack) == 0) {
            # No previous section
	    $section = "";
	    $in_section = 0;
	} else {
	    $section = pop @section_stack;
	}

	# call section close
	if ($this_section eq "POLY") {
	    print $outputfh "The phase of polymorphisms:\n";
	    print $outputfh join("\t", ("F1", "F2", "F3")),"\n";
	    print $outputfh join("\t", @premod),"\toutside presumed coding region (position $presumed_coding_start-$presumed_coding_end)\n";
	    print $outputfh join("\t", @mod),"\tin coding region\n";
	} 
	elsif ($this_section eq "INDELPOLY") 
       {
	   if ($indelpolys > 0) {
	       print $outputfh "Found $indelpolys indel polymorphisms.\n";
	   }
	}
	elsif ($this_section eq "GENOTYPE")
	{
	    print $outputfh "Found $heterozygote_positions instances of heterozygosity (positions in reads).\n";
	    print $outputfh "Found $nonsynonymous nonsynonymous and $synonymous synonymous differences from reference on polymorphic positions.\n";

	    print $outputfh "\tWorm\tnS\tnNS\n";
	    foreach my $worm (keys %variant) {
		my $wormsynpos=0;
		my $wormnonsynpos=0;
		
		foreach my $refpos (keys %{$variant{$worm}} ) {
		    if ( $variant{$worm}{$refpos} == 2) {
			$wormsynpos++;
		    } elsif ($variant{$worm}{$refpos} == 3) {
			$wormnonsynpos++;
		    }
		}
		print $outputfh "WORM\t$worm\t$wormsynpos\t$wormnonsynpos"; # farm
	    }
	}
    } elsif($in_section == 1) {
	
	# launch poly section..
	if ( $section eq "POLY" ) {
	    my @col;
	    (@col) = split(/\s+/, $r);
	    my $refpos = $col[1];

	    my $mod3 = ($refpos - $presumed_coding_start) % 3;
	    if ($refpos < $presumed_coding_start or $refpos > $presumed_coding_end) {
		$premod[$mod3]++;
	    } else {
		$mod[$mod3]++;

		my $var1=$col[3];
		my $var2=$col[4];
		my $aa1;
		my $aa2;
		my $printme = 0;
	       
		if ($mod3 == 0) {
		    my @ds = split(/ */,$col[5]);
		    			
		    my $codon1 = $var1 . $ds[0] . $ds[1];
		    my $codon2 = $var2 . $ds[0] . $ds[1];
		    
		    $aa1 = nt2aa($codon1);
		    $aa2 = nt2aa($codon2);

		    $printme = 1;
		}
		if ($mod3 == 1) {
		    my @us = split(/ */,$col[2]);
		    my @ds = split(/ */,$col[5]);
		    
		    my $codon1 =  $us[scalar(@us)-1] . $var1 . $ds[0];
		    my $codon2 =  $us[scalar(@us)-1] . $var2 . $ds[0];
		    
		    $aa1 = nt2aa($codon1);
		    $aa2 = nt2aa($codon2);

		    $printme = 1;
		}
		if ($mod3 == 2) {
		    my @us = split(/ */,$col[2]);

		    my $codon1 =  $us[scalar(@us)-2] . $us[scalar(@us)-1] . $var1 ;
		    my $codon2 =  $us[scalar(@us)-2] . $us[scalar(@us)-1] . $var2 ;
		    
		    $aa1 = nt2aa($codon1);
		    $aa2 = nt2aa($codon2);
		}

		if ($aa1 ne $aa2) {
		    $printme = 1;
		    print $outputfh join(" ",@col), " - phase $mod3 - $aa1 $aa2\n";
		}
	    }

	} elsif ( $section eq "CONTIG" ) {
	    if ($r=~m/^Contig/) {
		print $outputfh "Checking contig $r.\n";
	    } 
	    if ($r=~m/^REFERENCE\s+(.+)/) {
		print $outputfh "Contig has reference $1\n";
	    }
	} elsif ( $section eq "INDELPOLY" ) {
	    
	    $indelpolys++;
	}
	elsif ($section eq "GENOTYPE")
	{
	    # keep tally of syn/non-syn subs
	    # how treat heterozygotes?

	    # print all non-syn polys

	    # any division over motif already? or do that manually later?
	
	    my @col;
	    (@col) = split(/\s+/, $r);

	    my $refpos = $col[1];
	    my $mod3 = ($refpos - $presumed_coding_start) % 3;
		
	    my $refcodonpos = $refpos - $mod3;
	    my $refseqbases = $refseq->seq();
	    my $refcodon = substr($refseqbases, $refcodonpos-1, 3);
	    my $refbase = substr($refseqbases, $refpos-1,1);
	    my $refaa = nt2aa($refcodon);
	   	    
	    my $readname = $col[3];
	    my $readpos = $col[2];
	   
	    my ($worm) = ($readname =~ /^(\w{3}\d{4}\-\d{2})/);
	    
	    $var1 = $col[4];
	    $var2 = $col[5];
	    
	    my @refcodon = split(/ */, $refcodon);

	    if ( $var1 eq $refbase and $var2 eq $refbase) {
		
	    } else {
		# polymorphic..
	    }

	    if ($refpos < $presumed_coding_start or $refpos > $presumed_coding_end) {
		# ignored..
	    } else {
		# same as reference.. count for frequencies?
		my $aa1;
		my $aa2;

		my $codon1;
		my $codon2;

		# what about heterozygous positions next to each other?
		# ignore for now, as these are rare. 
		# as corollary, the surrounding read positions typically are not confidently called as polymorphic,
		# and should be treated as unreliable. We thus call the aa changes using the reference codon, plus
		# the read base from the actually called reference position. We thus no longer require reading reads.
		# And also don't have to solve the issue of reversed read coordinates, which seems to be a little ill
		# defined at the moment?
		
#		my $readseq = ${$reads{$readname}};
#		my $readseqbases = $readseq->seq();

		# note: possible to reach "outside" seq on early/late coords..
		if ($mod3 == 0) {
		    
		    # substr assumes 0 based coords, so -1 for the 1-based readpos
		    # then, we wish to start on the base after, so +1, and take 2 bases
#		    my $ds = substr( $readseqbases, $readpos +1 -1, 2 );
		    		    
		    $codon1 = $var1 . $refcodon[1] . $refcodon[2];
		    $codon2 = $var2 . $refcodon[1] . $refcodon[2];
		    
		    $aa1 = nt2aa($codon1);
		    $aa2 = nt2aa($codon2);
		}
		elsif ($mod3 == 1) {
		    
		    # substr assumes 0 based coords, so -1 for the 1-based readpos
		    # then, we wish to start on the base after, so +1, and take 1 base
#		    my $ds = substr( $readseqbases, $readpos +1 -1, 1 );
		    # then, we wish to start on the base before and, so -1, and take 1 base
#		    my $us=substr( $readseqbases, $readpos - 1 -1, 1 );

		    $codon1 = $refcodon[0] . $var1 . $refcodon[2];;
		    $codon2 = $refcodon[0] . $var2 . $refcodon[2];;
		    
		    $aa1 = nt2aa($codon1);
		    $aa2 = nt2aa($codon2);
		} 
		elsif ($mod3 == 2) {
		    # substr assumes 0 based coords, so -1 for the 1-based readpos
		    # then, we wish to start on the base two steps before and, so -2, and take 2 bases
#		    my $us=substr( $readseqbases, $readpos - 2 -1, 2 );

		    $codon1 = $refcodon[0] . $refcodon[1] . $var1;
		    $codon2 = $refcodon[0] . $refcodon[1] . $var2;

		    $aa1 = nt2aa($codon1);
		    $aa2 = nt2aa($codon2);
		}

		my $heterozygote =0;
		if( $var1 ne $var2 ) {
		    $heterozygote_positions++;
		    $heterozygote=1;
		    $DEBUG && print $outputfh "position $refpos - $readname heterozygous $var1 $var2 at readpos $readpos.\n";
		}

		if( $refaa eq $aa1 ) {
		    $DEBUG && print $outputfh "DEBUG: reference aa $refaa equals variant position, individual worm aa $aa1 for $readname, var 1.\n";
		    if (uc($codon1) ne uc($refcodon)) {
			$synonymous++;
			# log a synonymous variant
			if(defined($variant{$worm}{$refpos}) && $variant{$worm}{$refpos} >= 1 ) {
			    # this variant has previously been logged for this worm but from another read. ignore.
			} else {
			    $variant{$worm}{$refpos}=2;
			}
		    } else {
			# do nothing - no change in this copy
		    }
		} else {
		    $nonsynonymous++;
		    $DEBUG && print $outputfh "DEBUG: refpos $refpos, mod3 $mod3, refcodonpos $refcodonpos, refcodon $refcodon, refbase $refbase, refaa $refaa, length refseq obj ".($refseq->length()).", len refseqbases ".(length($refseqbases)).".\n";
		    $DEBUG && print $outputfh "Found codon 1 $codon1 ($aa1) and 2 $codon2 ($aa2) (P$mod3).\n";
		    print $outputfh "refpos $refpos, reference_aa $refaa($refcodon), worm_aa1 $aa1($codon1),  worm $worm, read $readname var1, readpos $readpos, frame F".($mod3+1).".\n";
		    # log a non-synonymous variant
		    $variant{$worm}{$refpos} = 3;
		    # allow for multiple variants?
		} 

		if ($refaa eq $aa2) {
		    $DEBUG && print $outputfh "DEBUG: reference aa $refaa equals variant position, individual worm aa $aa2 for $readname, var 2.\n";
		    if (uc($codon2) ne uc($refcodon)) {
			$synonymous++;

			if ( defined($variant{$worm}{$refpos}) && $variant{$worm}{$refpos} >= 2 )  {
			    # if variant is homozygous and already marked, no use marking again, right?
			} else {
			    # log a synonymous variant
			    $variant{$worm}{$refpos} = 2;
			}
		    } else {
                        # else do nothing - no change in this copy
		    }
		} else {
		    $nonsynonymous++;

		    # log a non-synonymous variant
		    $variant{$worm}{$refpos}=3;
		    # allow for multiple variants?

		    $DEBUG && print $outputfh "DEBUG: refpos $refpos, mod3 $mod3, refcodonpos $refcodonpos, refcodon $refcodon, refbase $refbase, refaa $refaa, length refseq obj ".($refseq->length()).", len refseqbases ".(length($refseqbases)).".\n";
		    $DEBUG && print $outputfh "Found codon 1 $codon1 ($aa1) and 2 $codon2 ($aa2) (P$mod3).\n";
		    print $outputfh "refpos $refpos, reference_aa $refaa($refcodon), worm_aa2 $aa2($codon2),  worm $worm, read $readname var2, readpos $readpos, frame F".($mod3+1).".\n";
		}
	    }
	}
    }
}

sub nt2aa {
    my $tri=$_[0];
    my $aa="";
    # Make nucleotides lowercase, and have all AAs uppercase.
    $tri=~tr/ATCGNX/atcgnx/;
    my $trilenthird = length($tri) / 3;
    
    # DEBUG: Strict??
    # print $DEBUG && "DEBUG: Length ", length($tri), "and the third is ", length($tri)/3 ,".\n";

    # ok, how much faster is a split -> case instead?    
    
    for(my $k = 0; $k < $trilenthird; $k++) {
	# Make trinucleotide codons... 
	$_=substr($tri,3*$k,3);

	# Genetic standard code according to NCBI
	#   AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	# Starts = ---M---------------M---------------M----------------------------
	# Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
	# Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
	# Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

	s/tt[tc]{1}/F/g;
	s/tt[ag]{1}/L/g;
	s/tc[atgcnx]{1}/S/g;
	s/ta[tc]{1}/Y/g;
	s/ta[ag]{1}/+/g;
	s/tg[tc]{1}/C/g;
	s/tga/+/g;
	s/tgg/W/g;
	s/ct[tcagnx]{1}/L/g;
	s/cc[tcagnx]{1}/P/g;
	s/ca[tc]{1}/H/g;
	s/ca[ag]{1}/Q/g;
	s/cg[tcagnx]{1}/R/g;
	s/at[tca]{1}/I/g;
	s/atg/M/g;
	s/ac[tcagnx]{1}/T/g;
	s/aa[tc]{1}/N/g;
	s/aa[ag]{1}/K/g;
	s/ag[tc]{1}/S/g;
	s/ag[ag]{1}/R/g;
	s/gt[tcagnx]{1}/V/g;
	s/gc[tcagnx]{1}/A/g;
	s/ga[tc]{1}/D/g;
	s/ga[ag]{1}/E/g;
	s/gg[atgcnx]{1}/G/g;
	s/[atgc]*[nx]+[atgc]*[nx]*/U/g; # only substrings of three at the time are processed
	
	# Append the aa the codon encoded to the virtual protein.
	$aa.=$_;
    }
    # Remove the extra nucleotides...

    $aa=~s/[atgcnx]//g;
    return $aa;
}
