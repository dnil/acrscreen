#!/usr/bin/perl -w
#
# upload-process-shotgun.pl
#
# Adjust the phred-made seq/qual fasta files to conseds liking
#
# Usage: upload-process-shotgun.pl seq_fasta_file_name qual_fasta_file_name phd_dir_name
#
# DN 001218
#

my $DEBUG = 0;

sub getFastaSeq;
sub getFastaQual;
sub getPhdChromTime;

# Check input arguments

if(@ARGV < 3) {
    die("Usage: upload-process-shotgun.pl seq_fasta_filename qual_fasta_filename phd_dir_name\n");
}

my $seq_fasta_filename=$ARGV[0];
my $qual_fasta_filename=$ARGV[1];
my $phd_dir_name=$ARGV[2];

# read in seq_fasta_file and qual_fasta_file

@_=getFastaSeq($ARGV[0]);
$nSeq=shift(@_);
@seq=@_[0..($nSeq-1)];
@name=@_[($nSeq)..(2*$nSeq-1)];

@_=getFastaQual($ARGV[1]);
$nQual=shift(@_);
@qual=@_;

if($nQual != $nSeq) {
    die "ERROR! Number of sequence entries does not equal number of quality entries.\n";
} 

# Read the phd dir for actual phd file names

print "Processing phd files...\n";

opendir(PHD_DIR, $phd_dir_name) || die "Could not open phd directory $phd_dir_name.\n";

my @phd_files = grep /.phd/, readdir PHD_DIR;

closedir(PHD_DIR);

my @found_seq_matching_phd;

foreach $phd_filename (@phd_files) {

    # Parse the phd-file for time and chrom entries
    my ($phd_chrom,$phd_time)=getPhdChromTime("$phd_dir_name/$phd_filename");    

    my ($phd_file_basename) = $phd_filename=~/(.+)\.phd.+/;
    
    my $phd_direction = "";
    ($phd_direction) = ($phd_filename=~/\-([FR]{1})\d+/);

    if ($phd_filename=~/(\.c1)/) {
	# reference sequence for polyscan/polyphred
	$phd_direction = "F";
    }
    
    
    if (! defined $phd_direction) {
	$phd_direction = "";
	
    }

    if ($phd_direction eq "F") {
	$phd_direction = "fwd";
    } elsif ($phd_direction eq "R") {
	$phd_direction = "rev";
    } 
    
    # print "DEBUG: processing phd file $phd_filename, basename $phd_file_basename\n"; 

    # Find any fasta name entries matching the phd-file basename 

    for($i=0; $i<$nSeq; $i++) {	
	if( $name[$i] =~ /$phd_file_basename\s+/ ) {
	    $DEBUG && print "DEBUG: Found match seq $name[$i] to phd_file_basename $phd_file_basename\n";
	    if(defined($found_seq_matching_phd[$i]) && $found_seq_matching_phd[$i] != 0 ) {
		print "Warning: found match nr $found_seq_matching_phd[$i] for seq $name[$i]; this time with phd file $phd_file_basename.\n";
	    }
	    if($name[$i]!~/CHROMAT_FILE/ ) {
		$name[$i].=" CHROMAT_FILE: $phd_chrom";		
	    }
	    if($name[$i]!~/PHD_FILE/ ) {
		$name[$i].=" PHD_FILE: $phd_filename";		
	    }
	    if($name[$i]!~/TIME/ ) {
		$name[$i].=" TIME: $phd_time";		
	    }
	    if($name[$i]!~/DIRECTION/ && $phd_direction ne "") {
		$name[$i].=" DIRECTION: $phd_direction";		
	    }

	    $found_seq_matching_phd[$i]++;
	}
    }    
}

print "Writing modified fasta files...\n";  

# Write the name-modified fasta files

open(SEQ_FILE,">$seq_fasta_filename.consedified");

for($i=0; $i<$nSeq; $i++) {
    print SEQ_FILE ">".$name[$i]."\n";
    print SEQ_FILE $seq[$i]."\n";       
}

close(SEQ_FILE);

open(QUAL_FILE,">$qual_fasta_filename.consedified"); 

for($i=0; $i<$nSeq; $i++) {   
    print QUAL_FILE ">".$name[$i]."\n";
    $q=join(' ',split(/,+/,$qual[$i]));
    print QUAL_FILE $q."\n";   
}

close(QUAL_FILE);

####### END MAIN #######

sub getFastaSeq {
  my $fastaFileName=shift;
  open(FASTAFILE,"<$fastaFileName") || die "Sequence fasta input file $fastaFileName open failed.\n";

  my @fasta;
  my @name;

  # First, get the sequences
  print "Reading sequences...\n";
  my $nFASTA=0;
  while(<FASTAFILE>) {
    chop;
    if(/^\>/) {
      # On fasta description line
      $nFASTA++;
      ($name[$nFASTA-1])=/^\>(.+)/;    
      $fasta[$nFASTA-1]="";
    } else {
      # Well, either the input is broken, or this is sequence data. Let us assume a sentient user.. :)
      # Get all genomic sequence chars into that $fasta string..
      s/[^atgcnxATGCNX]//g;
      $fasta[$nFASTA-1].=$_;
    }
  }

  # Done processing fasta sequence file
  close(FASTAFILE);

  return ($nFASTA, @fasta,@name); 
}

sub getFastaQual {
  my $nQual=0;

  my $fastaQualFileName=shift;
  open(QUALFILE, "<$fastaQualFileName") || die "Quality fasta input file $fastaQualFileName open failed\n";

  print "Reading quality values...\n";
  while(<QUALFILE>) {
    chop;
    if(/^\>/) {
      # On description line.. Name field should be pretty much equal to the fasta file, so ignore it.
      $nQual++;
      $qual[$nQual-1]="";
    } else {
      $qual[$nQual-1].=join(',',split(/\s+/,$_));
      $qual[$nQual-1].=",";
    }  
  }

  # Done processing quality file
  close(QUALFILE);
  
  return($nQual,@qual);
}

sub getPhdChromTime { 
    my $phd_filename=shift;
    
    open(PHD_FILE, "<$phd_filename") || die "Failed opening PHD file $phd_filename for input.\n";

    $chromat_file_found=0;
    $time_found=0;

    while(<PHD_FILE>) {
	chop;
	if(/CHROMAT_FILE:/) {
	    ($chromat_file)= m/CHROMAT_FILE:\s+(.+)\s*/;
	    $chromat_file_found=1;
	} elsif(/TIME:/) {
	    ($time) = m/TIME:\s+(.+)/;
	    $time_found=1;
	}
    } 

    close(PHD_FILE);

    if( ( ! $chromat_file_found ) ) {
	print "WARNING: no chromatogram file name entry found in $phd_filename.\n";
    }

    if( ( ! $time_found ) ) {
	print "WARNING: no time entry found in $phd_filename.\n";
    }

    return ($chromat_file, $time);

}








