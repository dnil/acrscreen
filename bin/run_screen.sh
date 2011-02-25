#!/bin/bash

# (c)2010 Daniel Nilsson, daniel.nilsson@izb.unibe.ch, daniel.k.nilsson@gmail.com
# Copyright 2010, held by the author.
# Released under the Perl Artistic License.
#
# Documentation available in this file. Use e.g. perldoc to view, or pod2man to create a manual.

### BEGIN USER SPECIFIC CONFIG

PHRAP=/home/daniel/src/phrap/phrap
PHRED=/home/daniel/src/phred.071220b/phred
if [ -z "$PHRED_PARAMETER_FILE" ]
then
    PHRED_PARAMETER_FILE=${PHRED%%phred}phredpar.dat
fi
export PHRED_PARAMETER_FILE

# phd2fasta from consed package misc
PHD2FASTA=/home/daniel/bin/phd2fasta
POLYPHRED=/home/daniel/src/polyphred-6.18-binary-x86_64-unknown-linux-gnu/bin/polyphred
SUDOPHRED=/home/daniel/src/polyphred-6.18-binary-x86_64-unknown-linux-gnu/bin/sudophred

CONSED=/home/daniel/src/consed/consed_linux64bit_static

CONSEDIFY=/home/daniel/sandbox/acrscreen/bin/consedify.pl
#/malariamotifs/bin/consedify.pl

REVSEQ=revseq
REVSEQPL=/home/daniel/sandbox/acrscreen/bin/revseq.pl
POLYPHREDPHD2IUPACFASTA=/home/daniel/sandbox/acrscreen/bin/polyphredphd2iupacfasta.pl

CHECKPOLY=/home/daniel/sandbox/acrscreen/bin/checkpoly.pl
PIPELINEFUNK=/home/daniel/sandbox/acrscreen/bin/pipelinefunk.sh

REFERENCE=/home/daniel/sandbox/acrscreen/eiire.fasta
REFSTRING=.c1

WORMFARMNAME=/home/daniel/hc/WormnumberFarmname.txt

### END CONFIG

### BEGIN DOCUMENTATION

: << 'END_OF_DOCS'

=head1 NAME

run_screen.sh

=head1 SYNOPSIS

 mkdir myproject
 cp *ab1 myproject/
 cd myproject
 ~/src/acrscreen/bin/run_screen.sh

=head1 DESCRIPTION

Process a set of chromatogram files with the phred-phrap-polyphred tools to produce some information on polymorphisms.

Also produces a set of IUPAC-aware fasta files from phred phd files.

=over 8 

=item 1

install and check environment variable paths for all required binaries 

=item 2 

generate and check environment variable path for reference sequence

=item 3

create project dir with project name

=item 4

copy all chromatograms (C<*ab1>) into this project dir

=item 5

C<cd> into the project dir

=item 6

 run_screen.sh

=item 7

 consed

=back

=head1 COPYRIGHT AND LICENSE

Copyright Daniel Nilsson, 2010. 

The package is released under the Perl Artistic License.

=head1 BUGS

Surely!

=head1 APPENDIX

The rest of the documentation details individual functions and quirks.

=cut

END_OF_DOCS

### END DOCUMENTATION

### BEGIN SCRIPT 

# will need extended globbing for file name truncations later
shopt -s extglob

# based on bash pipeline functions

. $PIPELINEFUNK

BASE=`basename $PWD`

if [ "$BASE" == "edit_dir" ]
then
    echo Cowardly refusing to run for project "edit_dir" - are you sure you are running this from the right directory? 
    exit
fi

for dir in chromat_dir edit_dir phd_dir poly_dir 
do
    if [ ! -d "$dir" ] 
    then
	mkdir $dir
    fi
done

#cd chromat_dir 
# for file in *ab1 ; do mv $file `echo $file | perl -ne 'chomp; @r = split(/[\_\-]+/); print join ("-",$r[0],$r[1]."_".$r[2],$r[3]."_".$r[4])."\n";'` ;done
	
if [  "" != "`ls -1 | grep ab1`"  ]
then
    mv *ab1 chromat_dir/
fi

#cd chromat_dir

# add direction info in read name..
#for chrom in *ab1 
#do 
#    if [[ "$chrom" =~ \-R[0-9]*\_ ]]
#    then
#	mv $chrom ${chrom%%.ab1}.g.ab1
#    elif [[ "$chrom" =~ \-F[0-9]*\_ ]]
#    then
#	mv $chrom ${chrom%%.ab1}.b.ab1
#    fi
#done

cd edit_dir

if needsUpdate ${BASE}.fasta ../chromat_dir $PHRED
then
    $PHRED -id ../chromat_dir -pd ../phd_dir -dd ../poly_dir -sa ${BASE}.fasta -qa ${BASE}.fasta.qual -trim_alt "" -trim_fasta -trim_phd -trim_cutoff 0.10
    registerFile $PWD/${BASE}.fasta result
    registerFile $PWD/${BASE}.fasta.qual result
fi

# no, turn that sequence around first..
if needsUpdate ${BASE}.ref.fasta.consedified ${BASE}.fasta $REFERENCE $SUDOPHRED $PHD2FASTA $CONSEDIFY
then
    $SUDOPHRED $REFERENCE -r $REFSTRING -q 30

    $PHD2FASTA -id ../phd_dir -os ${BASE}.ref.fasta -oq ${BASE}.ref.fasta.qual

    registerFile $PWD/$BASE.ref.fasta temp
    
    $CONSEDIFY ${BASE}.ref.fasta ${BASE}.ref.fasta.qual ../phd_dir 

    registerFile $PWD/$BASE.ref.fasta.consedified temp
    registerFile $PWD/$BASE.ref.fasta.qual.consedified temp

    ln -s $BASE.ref.fasta.qual.consedified ${BASE}.ref.fasta.consedified.qual
    registerFile $PWD/${BASE}.ref.fasta.consedified.qual temp
fi

acefile=${BASE}.ref.fasta.consedified.ace
if needsUpdate $acefile ${BASE}.ref.fasta.consedified ${BASE}.ref.fasta.qual.consedified $PHRAP 
then
    phrapout=${BASE}.ref.fasta.consedified.phrap_out
    $PHRAP -new_ace -forcelevel 10 ${BASE}.ref.fasta.consedified > $phrapout
    registerFile $PWD/$phrapout temp
    
    registerFile $PWD/$acefile result
fi

: << 'DOC_QUIRK'

=head2 Reference and consensus sequence sense not always identical

VERY annoyingly, I could not find a way to get phrap to respect the reference seq sense unconditionally, so sometimes it will be revved, and indel detection is then not available in polyphred. To be considered a low priority bug (as it currently works with this set of reads/reference ;).

=cut

DOC_QUIRK

refname=`grep \> $REFERENCE |sed -e 's/>//;'`
reforient=`grep "^AF ${refname}" ${BASE}.ref.fasta.consedified.ace  |cut -f 3 -d" "`
if [[ "${reforient}" == "C" ]] 
then
    revread="F"
    polyindel=""
else
    revread="R"
    polyindel="-indel 30"
fi

polyout=${BASE}.ref.fasta.consedified.polyout

if needsUpdate $polyout ${acefile} $POLYPHRED
then 
    # redo can use a clear..
    $POLYPHRED -clear

    # depending on naming style
    # $POLYPHRED -s /- -o ${BASE}.ref.fasta.consedified.polyout -ace ${BASE}.ref.fasta.consedified.ace -ref -indel 30

    $POLYPHRED -s 1 10 -o ${polyout} -ace ${acefile} -ref $REFSTRING $polyindel
    registerFile $PWD/$polyout result
fi 

polyoutcheck=${polyout}.check
if needsUpdate ${polyoutcheck} ${polyout} ${CHECKPOLY}
then 
    ${CHECKPOLY} --polyout $polyout --outfile ${polyout}.check --reference ${REFERENCE}
    registerFile $PWD/${polyoutcheck} result
fi

iupacfasta=${BASE}.all.iupac.fasta
if needsUpdate ${iupacfasta} ../phd_dir $POLYPHREDPHD2IUPACFASTA $REVSEQPL 
then

    if [ -e "${BASE}.all.iupac.fasta" ]
    then
        # remove old versions
	rm ${BASE}.all.iupac.fasta
	rm *.isolate.iupac.fasta
    fi

    for phdfile in ../phd_dir/*phd*;
    do 
	phdfilebase=${phdfile##*/}
	iupacfasta=${phdfilebase%%.phd*}.polyphred.iupac.fasta
	if needsUpdate $iupacfasta $phdfile $POLYPHREDPHD2IUPACFASTA $REVSEQPL
	then
	    $POLYPHREDPHD2IUPACFASTA $phdfile > $iupacfasta
	    
	    if [[ "$phdfile" =~ \-[0-9]?${revread}[0-9]*\_ ]]
	    then
		cp $iupacfasta ${iupacfasta}.rev
		
                #revseq ${iupacfasta}.rev $iupacfasta # is also ok...
		$REVSEQPL ${iupacfasta}.rev > $iupacfasta
		rm ${iupacfasta}.rev
	    fi

            # revcomp reference, if needed
	    if [[ "$phdfile" =~ $refname ]]
	    then
		if [[ "${reforient}" == "C" ]] 
		then
		    cp $iupacfasta ${iupacfasta}.rev
		    
                    #revseq ${iupacfasta}.rev $iupacfasta # is also ok...
		    $REVSEQPL ${iupacfasta}.rev > $iupacfasta
		    rm ${iupacfasta}.rev
		fi
	    fi
	fi
	registerFile ${PWD}/$iupacfasta temp

	cat ${iupacfasta} >> ${BASE}.all.iupac.fasta
	
        # requires shopt -s extglob 
	isolate=${phdfilebase%%-?([0-9])[FR][0-9]*}
	cat ${iupacfasta} >> $isolate.isolate.iupac.fasta

	registerFile ${PWD}/$isolate.isolate.iupac.fasta temp
    done

    registerFile $PWD/${BASE}.all.iupac.fasta result

    for isolatefasta in *isolate.iupac.fasta
    do
	cat $refname$REFSTRING.polyphred.iupac.fasta >> $isolatefasta
    done
fi

# log stats

log=$PWD/${BASE}.run_screen.log

if needsUpdate $log $polyoutcheck $wormfarmname
then
    cd ../chromat_dir

    # number of worms
    number_of_worms=`ls -1 |cut -d- -f1,2 |sort |uniq -c |grep -v EII|wc -l`

    # number of reads
    number_of_reads=`ls -1 *ab1 |wc -l`

    rundate=`date`
    echo "run_screen.sh - $rundate - $BASE" >> $log
    echo "Found $number_of_reads reads (electropherograms) from $number_of_worms worms." >> $log

    # list worms with max number of reads per worm
    echo >> $log
    echo "The worms with the highest number of reads per worm" >> $log
    echo "Reads   Worm" >> $log
    echo "=======|====" >> $log
    ls -1 |cut -d- -f1,2 |sort |uniq -c |grep -v EII|sort -k1nr |head -10 >> $log

    # list worms with min number of reads per worm
    echo >> $log
    echo "The worms with the lowest number of reads per worm" >> $log
    echo "Reads   Worm" >> $log
    echo "=======|====" >> $log
    ls -1 |cut -d- -f1,2 |sort |uniq -c |grep -v EII|sort -k1n |head -10 >> $log
    
    # number of reads per primer
    echo >> $log
    echo "Number of reads per primer" >> $log
    echo "Reads   Primer" >> $log
    echo "=======|======" >> $log


    for file in `ls -1`
    do
	if [[ "$file" =~ \-[0-9]?[FR][0-9]* ]]
	then
	    echo $file
	fi
    done | cut -f3 -d- |cut -d_ -f 1 |sort |uniq -c >> $log
   
    # number of reads per farm - req farm/worm mapping
    echo "Per worm nonsynonymous and synonymous variant positions" >> $log
    grep WORM $polyoutcheck | sort -k2,2 |awk 'BEGIN { n=0; s=0; ns=0; } { s=s+$3; ns=ns+$4; n=n+1} END { print "NS ",ns,"(",ns/n,") S ",s,"(",s/n,") N ",n; }' >> $log

    
    wormfarm=`basename $WORMFARMNAME`
    wormfarmnametab=${WORMFARMNAME%%.txt}.tab
    registerFile $PWD/$wormfarmnametab temp

    perl -ne 'm/(\w{3}\d{4}\-\d{2})\-[\w\d]+_[\w\d]+_([\w\d-]+)/; print "$1\t$2\n";' < $WORMFARMNAME |sort |uniq > $wormfarmnametab
    grep WORM $polyoutcheck |cut -f2,3,4 |sort> worm_s_ns
    join $wormfarmnametab worm_s_ns |cut -f 2,3,4 -d\ |sort | perl -e 'while (<>) { chomp; my @t=split /\s+/; $farm{$t[0]}{"ns"}+=$t[2]; $farm{$t[0]}{"s"}+=$t[1]; $farm{$t[0]}{"worms"}+=1; } map {print $_."\t".$farm{$_}{"s"}." (".sprintf("%.3f",$farm{$_}{"s"}/$farm{$_}{"worms"}).")\t".$farm{$_}{"ns"}." (".sprintf("%.3f",$farm{$_}{"ns"}/$farm{$_}{"worms"}).")\n";} keys %farm; print "Found ", scalar keys %farm," farms.\n";' >> $log
    rm worm_s_ns
    
fi

    # muscle -in $isolatefasta -html -out ${isolatefasta%%.fasta}.muscle.html
    # or split iupac, revcomp rev reads, merger, sixpack..select

### END SCRIPT

### LAUNCH CONSED

#$CONSED &
