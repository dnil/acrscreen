#!/usr/bin/perl -w

my $DEBUG = 0;

my $searchstring_given = 0;
my $searchstringfile = "";
my $fastafile = "";
my $match_whole_words = 0;

while (my $arg = shift @ARGV) {

    if ($arg =~ /^-/) {	
	if ($arg eq '-f') {
	    my $next_arg = shift @ARGV; 
	    if($next_arg eq "") {
		print "-f requires an argument, but non given. Bailing out.\n";
		exit 1;
	    } else {
		$searchstringfile = $next_arg;
		$searchstring_given = 1;
	    }
	}
	
	if ($arg eq '-w') {
	    $match_whole_words = 1;
	}

	if($arg eq '-v') { 
	    # TO BE IMPLEMENTED!
	}

    } else {
	if (! $searchstring_given ) {
	    $searchstring = $arg;
	    $searchstring_given = 1;
	} else {
	    $fastafile = $arg;
	}	
    }
}

if($fastafile eq "") {
    print "No fasta file to search was given. Nothing to do!\n";
    exit 1;
}

my @searchstrings =();

if($searchstringfile ne "") {
    open SEARCHSTRING, $searchstringfile or die  "Could not open $searchstringfile.\n";
    while(my $string = <SEARCHSTRING>) {
	chomp $string;
	push @searchstrings, $string;
    }
    close SEARCHSTRING;
    $DEBUG && print STDERR "Found ".scalar(@searchstrings)." search patterns.\n";
}

# assume considerably larger file than pattern list (else -v)

open FASTAFILE, $fastafile or die "Could not open $fastafile.\n";

my $printentry = 0;

while( my $row = <FASTAFILE> ) {

    chomp $row;
    if ( $row =~ m/^\>/ ) {
	$DEBUG && print STDERR "Checking $row.\n";
	if ( @searchstrings > 0 ) {
	    $printentry = 0;
	    for ( my $i = 0; $i < @searchstrings; $i++ ) {
		$searchstring = $searchstrings[$i];

		if( $match_whole_words == 1) {
		    if($row =~ m/\b$searchstring\b/) {
			print $row."\n";
			$printentry = 1;
			splice(@searchstrings, $i, 1);
			last;
		    }
		} elsif ( $match_whole_words == 0 ) {
		    if($row =~ m/$searchstring/) {
			print $row."\n";
			$printentry = 1;
			splice(@searchstrings, $i, 1);
			last;
		    }
		}
	    }
	} else {
	    if( $match_whole_words == 1) {
		if( $row =~ m/\b$searchstring\b/ ) {
		    print $row."\n";
		    $printentry = 1;
		} else {
		    $printentry = 0;
		}
	    } else {
		if( $row =~ m/$searchstring/ ) {
		    print $row."\n";
		    $printentry = 1;
		} else {
		    $printentry = 0;
		}
	    }
	}
    } elsif ( $printentry ) {
	print $row."\n";
    }
}
