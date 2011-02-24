#!/usr/bin/bash 

: << 'DOCUMENTATION'

=head1 NAME 

pipelinefunk.sh

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com

=head1 COPYRIGHT AND LICENSE

Copyright Daniel Nilsson, 2010. 

The package is released under the Perl Artistic License.

=head1 APPENDIX

The rest of the documentation details individual functions and quirks.

=cut

DOCUMENTATION

: << 'FUNCTION_DOC'

=head2 needsUpdate(target, prereq [, prereq]*)

Return true (needsupdate=yes) if target does not yet exist, is older than its prereqs or forceupdate=yes is in effect. set forceupdate=yes to run all available analyses, even if the file modification times advise against it 

=cut

FUNCTION_DOC

if [ -z "$forceupdate" ]
then
    forceupdate=no
fi

function needsUpdate()
{
    needsupdate="no"
    
    if [ "$forceupdate" = "yes" ] 
    then
	needsupdate="yes"
    fi

    target=$1;
    
    for prereq in ${@:2}
    do
	if [ $target -ot $prereq ]
	then
	    needsupdate="yes"
	fi
   done
    
    [ "$needsupdate" = "yes" ]
}

: << 'NUTB_FUNCTION_DOC'

=head2 needsUpdateTimeBased(target, timetoobsolete)

Return true (needsupdate=yes) if target does not yet exist, is older than timetoobsolete (in seconds) or forceupdate=yes is in effect. set forceupdate=yes to run all available analyses, even if the file modification times advise against it. 

    # sample synopsis
   
    seconds_in_two_days=$(( 60 * 60 * 24 * 2))
    update_pathways=no

    org_kegg_list=${org}.list

    if needsUpdateTimeBased ${org}.list $seconds_in_two_days
    then
	wget -c ftp://ftp.genome.jp/pub/kegg/pathway/organisms/${org}/${org}.list
	update_pathways=yes
	updates=yes
    fi

=cut

NUTB_FUNCTION_DOC

function needsUpdateTimeBased()
{
    local file=$1
    local timetoobsolete=$2
        
    local filestamp
    local nowstamp=`date +%s`

    local needsupdate="no"

    if [ "$forceupdate" = "yes" ] 
    then
	needsupdate=yes
    fi

    if [ ! -w $file ]
    then
	needsupdate=yes
    else
	# stat is used for timestamp retrieval, and works slightly differently on OSX
	os=`uname`
	if [ $os = "Darwin" ]
	then
	# osx
	    filestamp=`stat -f %c $file`
	elif [ $os = "Linux" ] 
	then
	# linux 
	    filestamp=`stat --format %Z $file`
	fi

	if [ $(( $nowstamp - $filestamp )) -gt $timetoobsolete ] 
	then
	    needsupdate=yes
	fi
    fi

    [ "$needsupdate" = "yes" ]
}

: <<'REGISTER_FUNCTION_DOC'

=head2 registerFile(file, category) 

 USAGE: registerFile /absolute/path/to/file category
 
Register file with the cleanup system. Basically log it to a file, named after the category. Category is currently limited to C<result|temp>.

=cut

REGISTER_FUNCTION_DOC

# on load, =)
# save current PWD for use with later regs.
pipelineregisterdir=$PWD

function registerFile()
{
    local savefile=$1
    local category=$2

    # test file..

    # check that it's not already on the list?
    
    echo $savefile >> ${pipelineregisterdir}/.pipeline.register.$category
    
}

: <<'CLEAN_FUNCTION_DOC'

=head2 cleanCategory(category)

 USAGE: cleanCategory category
 
Delete files registered with the cleanup system. Will recursively delete directories if registered. Category is currently limited to C<result|temp>.

=cut

CLEAN_FUNCTION_DOC

function cleanCategory()
{
    local category=$1
       
    for file in `cat ${pipelineregisterdir}/.pipeline.register.${category}`
    do
	if [ -d $file ] 
	then
	    rm -rf $file
	else
	    rm $file
	fi
    done

    rm ${pipelineregisterdir}/.pipeline.register.${category}
}

: << 'DOC_DIRECTIVE'

=head1 SYNOPSIS

 [DIRECTIVE='clean|shinyclen']
 . pipelinefunk.sh

If a directive of C<clean> or C<shinyclean> is set already when sourcing this, then clean accordingly. A directive can be given in C<$1> or C<$DIRECTIVE>.

=cut

DOC_DIRECTIVE

if [ -z "$DIRECTIVE" ]
then
    DIRECTIVE=$1
fi

if [ "$DIRECTIVE" == "clean" ]
then
    cleanCategory temp
fi

if [ "$DIRECTIVE" == "shinyclean" ]
then
    cleanCategory temp
    cleanCategory result
fi

if [ "$DIRECTIVE" == "onlyshinyclean" ]
then
    cleanCategory temp
    cleanCategory result
    exit
fi
