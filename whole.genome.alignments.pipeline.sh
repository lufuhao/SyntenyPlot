#!/bin/bash
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
if [ ! -z $(uname -m) ]; then
	machtype=$(uname -m)
elif [ ! -z "$MACHTYPE" ]; then
	machtype=$MACHTYPE
else
	echo "Warnings: unknown MACHTYPE" >&2
fi

#export NUM_THREADS=`grep -c '^processor' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`;
ProgramName=${0##*/}
echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"
RunPath=$PWD
echo "RunDir: $RunPath"



################# help message ######################################
help() {
cat<<HELP

$0 --- Brief Introduction

Version: 20161207

Requirements:
	fasta_splitter_by_numseq.pl
	LAST: lastdb, lastal, maf-convert
	UCSC: faSize, axtChain, chainMergeSort, chainPreNet, pslSwap
	
Descriptions:
	Detect synteny regions using LAST/ChainNet pipeline
	Output mGSV format: http://cas-bioinfo.cas.unt.edu/mgsv/

Output
	$o.axt
	$o

Options:
  -h    Print this help message
  -i    Fasta file to compare with each other except itself
  -g    Annotation GFF3 file only for genes
  -m    Repeat mask GFF3 file
  -r    Fasta file for reference sequences
  -db   Absolute path to lastdb indexed by lastdb using -r reference.fa
  -q    Fasta file for query sequences
  -o    Output.filtered.synteny.file, default: ./out.axt.filter
  -t    Number of threads
  -d    delete temporary files
  -s    Maximum gap to link for netToAxt, default: 1500
  -x    Lastal score metrix preset, default: HOXD70
  -psl  PSL file before doing pslSwap

  *Note specify (-i) or (-r and -q)

Example:
  $0 -i all.sequences.fa -g gene.gff3 -m repeatmask.gff3 -o ./my.synteny
  $0 -r reference.sequences.fa -q query.sequences.fa -g gene.gff3 -m repeatmask.gff3

Author:
  Fu-Hao Lu
  Post-Doctoral Scientist in Micheal Bevan laboratory
  Cell and Developmental Department, John Innes Centre
  Norwich NR4 7UH, United Kingdom
  E-mail: Fu-Hao.Lu@jic.ac.uk
HELP
exit 0
}
[ -z "$1" ] && help
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################
echo -e "\n######################\nProgram initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
threads=1
fastaall=''
referencefasta=''
queryfasta=''
gfffile=''
repeatgff=''
minisize=100
finaloutput="$RunPath/out.axt.filter"
cleantemporary=0
lastdbindex=''
minspace=25
maxgap=1500
pslfile=''
#lastal score metrix http://last.cbrc.jp/doc/last-matrices.html
#HOXD70 is medium. HoxD55 is far. human-chimp.v2 is close.
#Availble:
#	AT77		for weakly-similar AT-rich DNA (~27% substitutions and ~77% A+T) (MC Frith, NAR 2011 39(4):e23)
#	ATMAP		for strongly-similar AT-rich DNA (~4% substitutions and ~76% A+T)
#				for sequences with more than 4% substitution errors, if the excess error rate is explained by quality scores
#	BISF		for aligning bisulfite-converted DNA forward strands to a closely-related genome (MC Frith, R Mori, K Asai, NAR 2012 40(13):e100)
#	BISR		for aligning bisulfite-converted DNA reverse strands to a closely-related genome (MC Frith, R Mori, K Asai, NAR 2012 40(13):e100)
#	BL62 or BLOSUM62	protein scoring scheme is quite good at finding long-and-weak similarities, and not terrible at short-and-strong similarities (S Henikoff & JG Henikoff, PNAS 1992 89(22):10915-9)
#	BL80 or BLOSUM80	protein scoring scheme is good at finding somewhat short-and-strong similarities. (S Henikoff & JG Henikoff, PNAS 1992 89(22):10915-9)
#	HOXD70		often used for weak DNA similarities (F Chiaromonte, VB Yap, W Miller, PSB 2002:115-126)
#	MIQS		finding remote protein homologs (K Yamada & K Tomii, Bioinformatics 2014 30(3):317-25)
#	PAM30		good for finding short-and-strong similarities (MO Dayhoff et al. 1978)
scoremetrix='HOXD70'
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) fastaall=$2;shift 2;;
    -g) gfffile=$2;shift 2;;
    -m) repeatgff=$2;shift 2;;
    -r) referencefasta=$2;shift 2;;
    -db) lastdbindex=$2;shift 2;;
    -q) queryfasta=$2;shift 2;;
    -o) finaloutput=$2;shift 2;;
    -t) threads=$2;shift 2;;
    -d) cleantemporary=1;shift 1;;
    -s) maxgap=$2;shift 2;;
    -x) scoremetrix=$2;shift 2;;
    -psl) pslfile=$2;shift 2;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done


#################### Subfuctions ####################################
###Detect command existence
CmdExists () {
  if command -v $1 >/dev/null 2>&1; then
    echo 0;
  else
#    echo "I require $1 but it's not installed.  Aborting." >&2
    echo 1
  fi
#  local cmd=$1
#  if command -v $cmd >/dev/null 2>&1;then
#    echo >&2 $cmd "  :  "`command -v $cmd`
#    exit 0
#  else
#    echo >&2 "Error: require $cmd but it's not installed.  Exiting..."
#    exit 1
#  fi
}

###Usage: array=(`split delimiter string`)
split () {
	local separator=$1
	local mystring=$2
	echo $mystring | sed -e "s/$separator/\n/g"
}

#Usage: string=$(join delimiter array)
join () {
        local separator=$1
        shift 1
        local -a array=(`echo $@`)
        local returnstr=$(printf "$separator%s" "${array[@]}")
        returnstr=${returnstr:1}
        echo $returnstr
}

abs2rel () { perl -MFile::Spec -e 'print(File::Spec->abs2rel($ARGV[1], $ARGV[0]), "\n")' "$@"; }
#echo $(cd "$(dirname "$param")"; pwd)/$(basename "$param")
#################### Command test ###################################
echo -e "\n\n\n"
echo -e "\n\n\n" >&2
echo "### Command checking..."
echo "### command checking..." >&2
if [ $(CmdExists 'fasta_splitter_by_numseq.pl') -eq 1 ]; then
	echo "Error: script 'fasta_splitter_by_numseq.pl' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'lastdb') -ne 0 ]; then
	echo "Error: CMD 'lastdb' in program 'LAST' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'lastal') -ne 0 ]; then
	echo "Error: CMD 'lastal' in program 'LAST' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'maf-convert') -ne 0 ]; then
	echo "Error: CMD 'maf-convert' in program 'LAST' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'pslSwap') -ne 0 ]; then
	echo "Error: CMD 'pslSwap' in program UCSC KENT is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'faSize') -ne 0 ]; then
	echo "Error: CMD 'faSize' in program UCSC KENT is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'axtChain') -ne 0 ]; then
	echo "Error: CMD 'axtChain' in program UCSC KENT is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'chainMergeSort') -ne 0 ]; then
	echo "Error: CMD 'chainMergeSort' in program UCSC KENT is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'chainPreNet') -ne 0 ]; then
	echo "Error: CMD 'chainPreNet' in program UCSC KENT is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'pslSplitOnTarget') -ne 0 ]; then
	echo "Error: CMD 'pslSplitOnTarget' in program UCSC KENT is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'pslSplitOnTarget') -ne 0 ]; then
	echo "Error: CMD 'pslSplitOnTarget' in program UCSC KENT is required but not found.  Aborting..." >&2 
	exit 127
fi


#################### Defaults #######################################
splitdirref="$RunPath/1.refsplit"
splitdirquery="$RunPath/2.querysplit"
lastrundir="$RunPath/3.last"
lastdb_index='lastdb'
maf2psldir="$RunPath/4.maf2psl"
chaindir="$RunPath/5.chain"
netdir="$RunPath/6.net"

#################### Input and Output ###############################
if [ ! -z "$lastdbindex" ]; then
	lastdbindex=$(cd "$(dirname "$lastdbindex")"; pwd)/$(basename "$lastdbindex")
fi
echo -e "\n\n\n"
echo -e "\n\n\n" >&2
echo "### INPUT and OUTPUT checking..."
echo "### INPUT and OUTPUT checking..." >&2
if [ -z "$fastaall" ]; then
	if [ -z "$referencefasta" ] || [ ! -s "$referencefasta" ]; then
		echo "Error: invalid reference fasta" >&2;
		exit 1;
	fi
	if [ -z "$queryfasta" ] || [ ! -s "$queryfasta" ]; then
		echo "Error: invalid query fasta" >&2;
		exit 1;
	fi
	referencefasta=$(cd "$(dirname "$referencefasta")"; pwd)/$(basename "$referencefasta")
	echo "### Reference Fasta: $referencefasta"
	queryfasta=$(cd "$(dirname "$queryfasta")"; pwd)/$(basename "$queryfasta")
	echo "### Query Fasta: $queryfasta"
else
	fastaall=$(cd "$(dirname "$fastaall")"; pwd)/$(basename "$fastaall")
	echo "### InputFasta: $fastaall"
	referencefasta=$fastaall
	queryfasta=$fastaall
fi
if [[ "$threads" =~ [0-9]+ ]]; then
	echo "### Threads: $threads"
else
	echo "Error: invalid threads number (INT): $threads" >&2
	exit 1
fi

if [ -z "$gfffile" ] || [ ! -s "$gfffile" ]; then
	echo "### Gene annotation GFF files NOT detected, skip..."
else
	gfffile=$(cd "$(dirname "$gfffile")"; pwd)/$(basename "$gfffile")
fi

if [ -z "$repeatgff" ] || [ ! -s "$repeatgff" ]; then
	echo "### Repeatmask GFF files NOT detected, skip..."
else
	repeatgff=$(cd "$(dirname "$repeatgff")"; pwd)/$(basename "$repeatgff")
fi



#################### Cleaning ###############################
CleanTemporary () {
	echo "Info: Cleanning Temporary files"
	if [ -d $splitdirref ]; then
		rm -rf $splitdirref
	fi
	if [ -d $splitdirquery ]; then
		rm -rf $splitdirquery
	fi
	if [ -d $lastrundir ]; then
		rm -rf $lastrundir
	fi
	if [ -d $maf2psldir ]; then
		rm -rf $maf2psldir
	fi
	if [ -d "$maf2psldir.swap" ]; then
		rm -rf "$maf2psldir.swap"
	fi
	if [ -d $chaindir ]; then
		rm -rf $chaindir
	fi
	if [ -d $netdir ]; then
		rm -rf $netdir
	fi
}

CleanTemporary



#################### Main ###########################################
echo -e "\n\n\n"
echo -e "\n\n\n" >&2
echo "### Program starts..."
echo "### Program starts..." >&2



### split reference dir
echo -e "\n"
echo -e "\n" >&2
if [ -z "$pslfile" ] || [ ! -s "$pslfile" ]; then
	echo "### Step1: split reference ..."
	echo "### Step1: split reference ..." >&2
	if [ -d "$splitdirref" ]; then
		rm -rf $splitdirref
	fi
	mkdir -p $splitdirref
	cd $splitdirref
	fasta_splitter_by_numseq.pl $referencefasta 1 'sequence'
	if [ $? -ne 0 ]; then
		echo "Error: split reference failed: $referencefasta" >&2
		exit 1
	fi
else
	echo "### Step1: split reference skiped as PSL file specified: -psl $pslfile"
	echo "### Step1: split reference skiped as PSL file specified: -psl $pslfile" >&2
fi



### split query dir
echo -e "\n"
echo -e "\n" >&2
if [ -z "$pslfile" ] || [ ! -s "$pslfile" ]; then
	echo "### Step2: split query ..."
	echo "### Step2: split query ..." >&2
	if [ -d "$splitdirquery" ]; then
		rm -rf $splitdirquery
	fi
	mkdir -p $splitdirquery
	cd $splitdirquery
	fasta_splitter_by_numseq.pl $queryfasta 1 'sequence'
	if [ $? -ne 0 ]; then
		echo "Error: split query failed: $queryfasta" >&2
		exit 1
	fi
else
	echo "### Step2: split query skiped as PSL file specified: -psl $pslfile"
	echo "### Step2: split query skiped as PSL file specified: -psl $pslfile" >&2
fi



### last
echo -e "\n"
echo -e "\n" >&2
if [ -z "$pslfile" ] || [ ! -s "$pslfile" ]; then
	echo "### Step3: running LAST to detect syntonic region ..."
	echo "### Step3: running LAST to detect syntonic region ..." >&2
	if [ -d $lastrundir ]; then
		rm -rf $lastrundir
	fi
	mkdir -p $lastrundir
	cd $lastrundir
	for refseq in `ls $splitdirref/sequence.*.fa`; do
		if [ ! -z "$lastdbindex" ]; then
			echo "Info: using existing lastdb index: $lastdbindex"
		else
			rm -rf ${lastdb_index}* > /dev/null 2>&1
			lastdb -c -R10 -Q 0 -P $threads $lastdb_index $refseq
			if [ $? -ne 0 ]; then
				echo "Error: index LAST DB failed: $refseq" >&2
				exit 1
			fi
			lastdbindex="$lastrundir/$lastdb_index"
			echo "Info: using generated lastdb index: $lastdbindex"
		fi

		for queryseq in `ls $splitdirquery/sequence.*.fa`; do
			seqname01=${refseq##*/}
			seqname02=${queryseq##*/}
			seqbase01=${seqname01%.*}
			seqbase02=${seqname02%.*}
			if [ "$seqname01" == "$seqname02" ] && [ ! -z "$fastaall" ]; then
				rm $splitdirquery/$seqname01
				continue
				###Compare to ifsef
			fi
			if [ -s ref.$seqbase01.query.$seqbase02.maf ]; then
				echo "Error: last output exists: ref.$seqbase01.query.$seqbase02.maf" >&2
				exit 1
			fi
			lastal -f MAF -Q 0 -p $scoremetrix -M -P $threads $lastdbindex $queryseq > ref.$seqbase01.query.$seqbase02.maf
			if [ $? -ne 0 ] || [ ! -s ref.$seqbase01.query.$seqbase02.maf ]; then
				echo "Error: lastal failed: ref $refseq, query $queryseq" >&2
				exit 1
			fi
		done
	done
else
	echo "### Step3: running LAST skiped as PSL file specified: -psl $pslfile"
	echo "### Step3: running LAST skiped as PSL file specified: -psl $pslfile" >&2
fi



### maf2psl
echo -e "\n"
echo -e "\n" >&2
if [ -z "$pslfile" ] || [ ! -s "$pslfile" ]; then
	echo "### Step4: convert MAF to PSL ..."
	echo "### Step4: convert MAF to PSL ..." >&2
	if [ -d $maf2psldir ]; then
		rm -rf $maf2psldir
	fi
	mkdir -p $maf2psldir
	cd $maf2psldir

	for maffile in `ls $lastrundir/*.maf`; do
		mafname=${maffile##*/}
		mafbase=${mafname%.*}
		maf-convert psl $maffile > $mafbase.psl
		if [ $? -ne 0 ] || [ ! -s $mafbase.psl ]; then
			echo "Error: maf2psl: $maffile -> $mafbase.psl" >&2
			exit 1
		fi
	done
	cat *.psl > $finaloutput.psl
	if [ $? -ne 0 ] || [ ! -s "$finaloutput.psl" ]; then
		echo "Error: collect all psl: $finaloutput.psl" >&2
		exit 1
	fi
	pslfile=$(echo $(cd $(dirname "$finaloutput.psl"); pwd)/$(basename "$finaloutput.psl"))
else
	echo "### Step4: convert MAF to PSL skiped as PSL file specified: -psl $pslfile"
	echo "### Step4: convert MAF to PSL skiped as PSL file specified: -psl $pslfile" >&2
fi



echo -e "\n"
echo -e "\n" >&2
echo "### Step5: PSL swap ..."
echo "### Step5: PSL swap ..." >&2
pslSwap $pslfile all.swap.psl
if [ $? -ne 0 ] || [ ! -s "all.swap.psl" ]; then
	echo "Error: pslSwap" >&2
	exit 1
fi

if [ -d "$maf2psldir.swap" ]; then
	rm -rf "$maf2psldir.swap"
fi
mkdir -p "$maf2psldir.swap"
pslSplitOnTarget all.swap.psl "$maf2psldir.swap"/ -lump
if [ $? -ne 0 ]; then
	echo "Error: pslSplitOnTarget" >&2
	exit 1
fi



### chain
echo -e "\n"
echo -e "\n" >&2
echo "### Step6: chaining ..."
echo "### Step6: chaining ..." >&2
if [ -d $chaindir ]; then
	rm -rf $chaindir
fi
mkdir -p $chaindir
cd $chaindir
reference2bit=$chaindir/reference.2bit
faToTwoBit $referencefasta $reference2bit
if [ $? -ne 0 ] || [ ! -s $reference2bit ]; then
	echo "Error: faToTwoBit reference: $referencefasta" >&2
	exit 1
fi
query2bit=$chaindir/query.2bit
faToTwoBit $queryfasta $query2bit
if [ $? -ne 0 ] || [ ! -s $query2bit ]; then
	echo "Error: faToTwoBit query: $queryfasta" >&2
	exit 1
fi

faSize -detailed $referencefasta > reference.size
if [ $? -ne 0 ] || [ ! -s $chaindir/reference.size ]; then
	echo "Error: faSize reference: $referencefasta" >&2
	exit 1
fi
faSize -detailed $queryfasta > query.size
if [ $? -ne 0 ] || [ ! -s $chaindir/query.size ]; then
	echo "Error: faSize query: $queryfasta" >&2
	exit 1
fi
for pslfile in `ls "$maf2psldir.swap"/*.psl`; do
	pslname=${pslfile##*/}
	pslbase=${pslname%.*}
#	axtChain -psl -faQ -faT -linearGap=loose $pslfile $queryfasta $referencefasta $pslbase.chain
	axtChain -psl -linearGap=loose $pslfile $query2bit $reference2bit $pslbase.chain
	if [ $? -ne 0 ] || [ ! -s $pslbase.chain ]; then
		echo "Error: axtChain: $pslfile" >&2
		exit 1
	fi
done
echo "### Step6: chainMergeSort ..."
echo "### Step6: chainMergeSort ..." >&2
chainMergeSort $chaindir/*.chain > all.chain
if [ $? -ne 0 ] || [ ! -s $chaindir/all.chain ]; then
	echo "Error: chainMergeSort: $chaindir/all.chain" >&2
	exit 1
fi
echo "### Step6: chainPreNet ..."
echo "### Step6: chainPreNet ..." >&2
chainPreNet $chaindir/all.chain $chaindir/query.size $chaindir/reference.size  all.pre.chain
if [ $? -ne 0 ] || [ ! -s $chaindir/all.pre.chain ]; then
	echo "Error: chainPreNet: $chaindir/all.chain" >&2
	exit 1
fi



### Netting
echo -e "\n"
echo -e "\n" >&2
echo "### Step7: neting ..."
echo "### Step7: neting ..." >&2
if [ -d $netdir ]; then
	rm -rf $netdir
fi
mkdir -p $netdir
cd $netdir
chainNet -minSpace=1 $chaindir/all.pre.chain $chaindir/query.size $chaindir/reference.size stdout /dev/null | netSyntenic stdin noClass.net
if [ $? -ne 0 ] || [ ! -s $netdir/noClass.net ]; then
	echo "Error:chainNet: $netdir/noClass.net" >&2
	exit 1
fi
#netClass -noAr noClass.net ci2 cioSav2 cioSav2.net

echo "### Step7: netToAxt ..."
echo "### Step7: netToAxt ..." >&2
#netToAxt -maxGap=$maxgap noClass.net $chaindir/all.pre.chain $query2bit $reference2bit $finaloutput.axt
netToAxt -maxGap=$maxgap noClass.net $chaindir/all.pre.chain $query2bit $reference2bit stdout | axtSort stdin $finaloutput.axt
if [ $? -ne 0 ] || [ ! -s $finaloutput.axt ]; then
	echo "Error: netToAxt: $finaloutput.axt" >&2
	exit 1
fi



### filter
echo -e "\n"
echo -e "\n" >&2
echo "### Step8: Preparing Final synteny file ..."
echo "### Step8: Preparing Final synteny file ..." >&2
export MINISIZE=$minisize
perl -e 'print "#org1\torg1_start\torg1_end\torg2\torg2_start\torg2_end\tscore\tevalue\n"' > $finaloutput
#perl -ne 'BEGIN{$minisize=$ENV{"MINISIZE"};}chomp; next unless (/^\d+/); @arr=split(/\s+/); print $arr[1], "\t", $arr[2], "\t", $arr[3], "\t", $arr[4], "\t", $arr[5], "\t", $arr[6], "\t", $arr[8], "\t", 0, "\n" if(($arr[3]-$arr[2])>=$minisize and ($arr[6]-$arr[5])>=$minisize);' $finaloutput.axt >> $finaloutput
perl -ne 'BEGIN{$minisize=$ENV{"MINISIZE"};}chomp; next unless (/^\d+/); @arr=split(/\s+/); next unless (($arr[3]-$arr[2])>=$minisize or ($arr[6]-$arr[5])>=$minisize); print $arr[1], "\t", $arr[2], "\t", $arr[3], "\t", $arr[4], "\t", $arr[5], "\t", $arr[6], "\t", $arr[8], "\t", $arr[7], "\n";' $finaloutput.axt >> $finaloutput


if [ $? -ne 0 ] || [ ! -s $finaloutput ]; then
	echo "Error: final axt: $finaloutput" >&2
	exit 1
fi


### Annotation
cd $RunPath
if [ ! -z "$gfffile" ] && [ -s $gfffile ]; then
	echo -e "\n\n\n"
	echo -e "\n\n\n" >&2
	echo "### Step9: Preparing gene annotation ..."
	echo "### Step9: Preparing gene annotation ..." >&2
	echo "### Detect Gene annotation GFF files: $gfffile"
	perl -e 'print "#org_id\tstart\tend\tstrand\tfeature_name\tfeature_value\ttrack_name\ttrack_shape\ttrack_color\n";' > $finaloutput.annot.gff3
	perl -ne 'chomp;next if (/^#/);@arr=split(/\t/);$arr[8]=~s/^.*ID=//;$arr[8]=~s/;.*$//; print "$arr[0]\t$arr[3]\t$arr[4]\t$arr[6]\t$arr[8]\t.\tgene\tarrow\tbrown\n";' $gfffile >> $finaloutput.annot.gff3
else
	echo "### Gene annotation GFF files NOT detected, skip..."
fi


### repeatmask
if [ ! -z "$repeatgff" ] && [ -s "$repeatgff" ]; then
	echo "### Detect repeatmask GFF files: $repeatgff"
	#perl -MFuhaoPerl5Lib::MiscKit=MergeRanges -MData::Dumper=Dumper -ne 'BEGIN{%hash=();} chomp; next if (/^#/); @arr=split(/\t/);print "Tewst: @arr\n"; $hash{$arr[0]}{$arr[3]}=$arr[4]; END {print Dumer $hash{$seq};foreach $seq (sort (keys %hash)){($test, $id)=MergeRanges($hash{$seq}); print STDERR "Test: $seq\t$test\n"; foreach $j (sort {$a<=>$b} keys %{$id}) {print $seq, "\t", $j, "\t", ${$id}{$j}, "\n";}}}' $repeatgff
	perl -MFuhaoPerl5Lib::MiscKit=MergeRanges -ne 'BEGIN{%hash=();} chomp; next if (/^#/); @arr=split(/\t/);$hash{$arr[0]}{$arr[3]}=$arr[4]; END {foreach $seq (sort (keys %hash)){($test, $id)=MergeRanges($hash{$seq}); print STDERR "Test: $seq\t$test\n"; foreach $j (sort {$a<=>$b} keys %{$id}) {print $seq, "\t", $j, "\t", ${$id}{$j}, "\t.\t.\t.\trepeat\tbox\tblack\n";}}}' $repeatgff >> $finaloutput.annot.gff3
else
	echo "### Repeatmask GFF files NOT detected, skip..."
fi



if [ $cleantemporary -eq 1 ]; then
	CleanTemporary
fi



exit 0;
