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

Version: 20170110

Requirements:
    fasta_splitter_by_numseq.pl
    LAST: lastdb, lastal, maf-convert
    UCSC: faSize, axtChain, chainMergeSort, chainPreNet, pslSwap,
          pslSplitOnTarget, faToTwoBit, netToAxt
    
Descriptions:
    Detect synteny regions using LAST/ChainNet pipeline
    Output mGSV format: http://cas-bioinfo.cas.unt.edu/mgsv/

Output
    \$o.axt
    \$o

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
  -s    Maximum gap to link for netToAxt, default: 100
  -x    Lastal score metrix preset, default: HOXD70
  -psl  PSL file before doing pslSwap
  -l    Minimum alignment length to filter, default: 100
  -LG   axtChain -linearGap option, default: medium
        <medium|loose|filename> Specify type of linearGap to use.
        loose is chicken/human linear gap costs.
        medium is mouse/human linear gap costs.
        Or specify a piecewise linearGap tab delimited file.

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
min_align_length=100
finaloutput="$RunPath/out.axt.filter"
cleantemporary=0
lastdbindex=''
minspace=25
maxgap=100
pslfile=''
linearGap='medium'
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
    -l) min_align_length=$2;shift 2;;
    -LG) linearGap=$2;shift 2;;
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
File2Fullpath () {
	local FFinput=$1
	local FFsubinfo="SUF(File2Fullpath)"
	
#	echo "#Test1: $FFinput #" >&2
	if [[ "$FFinput" =~ .*\/.* ]]; then
		if [[ "$FFinput" =~ ^\/ ]]; then
#			echo "#Test2: $FFinput #" >&2
			echo $FFinput
		else
#			echo "#Test3: $FFinput #" >&2
			echo $(cd $(dirname "$FFinput"); pwd)/$(basename "$FFinput")
		fi
	elif [[ "$FFinput" =~ ^[-0-9a-zA-Z\._]+$ ]]; then
#		echo "#Test4: $FFinput #" >&2
		echo $(pwd)/$(basename "$FFinput")
	else
#		echo "#Test5: $FFinput #" >&2
		echo ''
	fi
#	echo "#Test6: $FFinput #" >&2
}
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
if [ $(CmdExists 'faToTwoBit') -ne 0 ]; then
	echo "Error: CMD 'faToTwoBit' in program UCSC KENT is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'netToAxt') -ne 0 ]; then
	echo "Error: CMD 'netToAxt' in program UCSC KENT is required but not found.  Aborting..." >&2 
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
finaloutput=$(File2Fullpath $finaloutput)
if [ ! -z "$pslfile" ]; then
	pslfile=$(File2Fullpath $pslfile)
	if [ ! -s "$pslfile" ]; then
		echo "Error: invalid PSL file" >&2
		exit 100
	else
		echo "Info: PSL file detected: $pslfile"
	fi
fi



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
if [[ $min_align_length =~ ^[0-9]+$ ]]; then
	echo "Parameter -l : $min_align_length" >/dev/null
else
	echo "Error: invalid INT -l parameter" >&2
	exit 10;
fi
if [[ "$linearGap" == "medium" ]]; then
	echo "Info: axtChain -linearGap was set to $linearGap"
elif [[ "$linearGap" == "loose" ]]; then
	echo "Info: axtChain -linearGap was set to $linearGap"
elif [ -s "$linearGap" ]; then
	echo "Info: axtChain -linearGap was set to file $linearGap"
	linearGap=$(File2Fullpath $linearGap)
else
	echo "Error: invalid -LG option for axtChain -linearGap" >&2
	exit 100;
fi




#################### Cleaning ###############################
CleanTemporary () {
	echo "Info: Cleanning Temporary files"
	if [ -d $splitdirref ]; then
		rm -rf $splitdirref > /dev/null 2>&1
	fi
	if [ -d $splitdirquery ]; then
		rm -rf $splitdirquery > /dev/null 2>&1
	fi
	if [ -d $lastrundir ]; then
		rm -rf $lastrundir > /dev/null 2>&1
	fi
	if [ -d $maf2psldir ]; then
		rm -rf $maf2psldir > /dev/null 2>&1
	fi
	if [ -d "$maf2psldir.swap" ]; then
		rm -rf "$maf2psldir.swap" > /dev/null 2>&1
	fi
	if [ -d $chaindir ]; then
		rm -rf $chaindir > /dev/null 2>&1
	fi
	if [ -d $netdir ]; then
		rm -rf $netdir > /dev/null 2>&1
	fi
}



### Clean previous folder and files
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
		rm -rf $splitdirref > /dev/null 2>&1
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
		rm -rf $splitdirquery > /dev/null 2>&1
	fi
	mkdir -p $splitdirquery
	cd $splitdirquery
	fasta_splitter_by_numseq.pl $queryfasta 1 'sequence'
	if [ $? -ne 0 ]; then
		echo "Error: split query failed: $queryfasta" >&2
		echo "CMD used: fasta_splitter_by_numseq.pl $queryfasta 1 'sequence'" >&2
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
		rm -rf $lastrundir > /dev/null 2>&1
	fi
	mkdir -p $lastrundir
	cd $lastrundir
	for refseq in `find $splitdirref -type f -name sequence.*.fa`; do
		finalindex=''
		if [ ! -z "$lastdbindex" ]; then
			echo "Info: using existing lastdb index: $lastdbindex"
			finalindex=$lastdbindex
		else
			rm -rf ${lastdb_index}* > /dev/null 2>&1
			lastdb -c -R10 -Q 0 -P $threads $lastdb_index $refseq
			if [ $? -ne 0 ]; then
				echo "Error: index LAST DB failed: $refseq" >&2
				echo "CMD used: lastdb -c -R10 -Q 0 -P $threads $lastdb_index $refseq" >&2
				exit 1
			fi
			finalindex="$lastrundir/$lastdb_index"
			echo "Info: using generated lastdb index: $finalindex"
		fi
		for queryseq in `find $splitdirquery -type f -name sequence.*.fa`; do
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
			lastal -f MAF -Q 0 -p $scoremetrix -M -P $threads $finalindex $queryseq > ref.$seqbase01.query.$seqbase02.maf
			if [ $? -ne 0 ]; then
				echo "Error: lastal failed: ref $refseq, query $queryseq" >&2
				echo "CMD used: lastal -f MAF -Q 0 -p $scoremetrix -M -P $threads $finalindex $queryseq > ref.$seqbase01.query.$seqbase02.maf" >&2
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
		rm -rf $maf2psldir > /dev/null 2>&1
	fi
	mkdir -p $maf2psldir
	cd $maf2psldir
	pslfile=$(File2Fullpath "$finaloutput.psl")
	if [ -e "$pslfile" ]; then
		rm -rf $pslfile > /dev/null 2>&1
	fi
	touch "$pslfile"

	for maffile in `find $lastrundir -type f -name *.maf`; do
		mafname=${maffile##*/}
		mafbase=${mafname%.*}
		maf-convert psl $maffile > $mafbase.psl
		if [ $? -ne 0 ]; then
			echo "Error: maf2psl: $maffile -> $mafbase.psl" >&2
			echo "CMD used: maf-convert psl $maffile > $mafbase.psl" >&2
			exit 100
		fi
		if [ -s "$mafbase.psl" ]; then
			cat "$mafbase.psl" >> "$pslfile"
			rm -f $maffile "$mafbase.psl" > /dev/null 2>&1
		fi
	done
	if [ $? -ne 0 ] || [ ! -s "$pslfile" ]; then
		echo "Error: collect all psl: $pslfile" >&2
		exit 100
	fi
	
else
	echo "### Step4: convert MAF to PSL skiped as PSL file specified: -psl $pslfile"
	echo "### Step4: convert MAF to PSL skiped as PSL file specified: -psl $pslfile" >&2
fi



echo -e "\n"
echo -e "\n" >&2
echo "### Step5: PSL swap ..."
echo "### Step5: PSL swap ..." >&2
if [ ! -d $maf2psldir ]; then
	mkdir -p $maf2psldir
fi
cd $maf2psldir
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
	echo "CMD used: pslSplitOnTarget all.swap.psl $maf2psldir.swap/ -lump" >&2
	exit 1
fi



### chain
echo -e "\n"
echo -e "\n" >&2
echo "### Step6: chaining ..."
echo "### Step6: chaining ..." >&2
if [ -d $chaindir ]; then
	rm -rf $chaindir > /dev/null 2>&1
fi
mkdir -p $chaindir
cd $chaindir
reference2bit=$chaindir/reference.2bit
faToTwoBit $referencefasta $reference2bit
if [ $? -ne 0 ] || [ ! -s $reference2bit ]; then
	echo "Error: faToTwoBit reference: $referencefasta" >&2
	echo "CMD used: faToTwoBit $referencefasta $reference2bit" >&2
	exit 1
fi
query2bit=$chaindir/query.2bit
faToTwoBit $queryfasta $query2bit
if [ $? -ne 0 ] || [ ! -s $query2bit ]; then
	echo "Error: faToTwoBit query: $queryfasta" >&2
	echo "CMD used: faToTwoBit $queryfasta $query2bit" >&2
	exit 1
fi

faSize -detailed $referencefasta > reference.size
if [ $? -ne 0 ] || [ ! -s $chaindir/reference.size ]; then
	echo "Error: faSize reference: $referencefasta" >&2
	echo "CMD used: faSize -detailed $referencefasta > reference.size" >&2
	exit 1
fi
faSize -detailed $queryfasta > query.size
if [ $? -ne 0 ] || [ ! -s $chaindir/query.size ]; then
	echo "Error: faSize query: $queryfasta" >&2
	echo "CMD used: faSize -detailed $queryfasta > query.size" >&2
	exit 1
fi
for pslfile in `find "$maf2psldir.swap" -type f -name *.psl`; do
	pslname=${pslfile##*/}
	pslbase=${pslname%.*}
#	axtChain -psl -faQ -faT -linearGap=$linearGap $pslfile $queryfasta $referencefasta $pslbase.chain
	axtChain -psl -linearGap=$linearGap $pslfile $query2bit $reference2bit $pslbase.chain
	if [ $? -ne 0 ] || [ ! -s $pslbase.chain ]; then
		echo "Error: axtChain: $pslfile" >&2
		echo "CMD used: axtChain -psl -linearGap=loose $pslfile $query2bit $reference2bit $pslbase.chain" >&2
		exit 1
	fi
done
echo "### Step6: chainMergeSort ..."
echo "### Step6: chainMergeSort ..." >&2
chainMergeSort $chaindir/*.chain > all.chain
if [ $? -ne 0 ] || [ ! -s $chaindir/all.chain ]; then
	echo "Error: chainMergeSort: $chaindir/all.chain" >&2
	echo "CMD used: chainMergeSort $chaindir/*.chain > all.chain" >&2
	exit 1
fi
echo "### Step6: chainPreNet ..."
echo "### Step6: chainPreNet ..." >&2
chainPreNet $chaindir/all.chain $chaindir/query.size $chaindir/reference.size  all.pre.chain
if [ $? -ne 0 ] || [ ! -s $chaindir/all.pre.chain ]; then
	echo "Error: chainPreNet: $chaindir/all.chain" >&2
	echo "CMD used: chainPreNet $chaindir/all.chain $chaindir/query.size $chaindir/reference.size  all.pre.chain" >&2
	exit 1
fi



### Netting
echo -e "\n"
echo -e "\n" >&2
echo "### Step7: neting ..."
echo "### Step7: neting ..." >&2
if [ -d $netdir ]; then
	rm -rf $netdir > /dev/null 2>&1
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
#netToAxt -maxGap=$maxgap noClass.net $chaindir/all.pre.chain $query2bit $reference2bit ${finaloutput%.*}.ORIGIN.axt
netToAxt -maxGap=$maxgap noClass.net $chaindir/all.pre.chain $query2bit $reference2bit stdout | axtSort stdin ${finaloutput%.*}.ORIGIN.axt
if [ $? -ne 0 ] || [ ! -s "${finaloutput%.*}.ORIGIN.axt" ]; then
	echo "Error: netToAxt: ${finaloutput%.*}.ORIGIN.axt" >&2
	exit 1
fi



### filter
echo -e "\n"
echo -e "\n" >&2
echo "### Step8: Preparing Final synteny file ..."
echo "### Step8: Preparing Final synteny file ..." >&2
export MINISIZE=$min_align_length
perl -e 'print "#org1\torg1_start\torg1_end\torg2\torg2_start\torg2_end\tscore\tevalue\n"' > $finaloutput
#perl -ne 'BEGIN{$minisize=$ENV{"MINISIZE"};; print STDERR "Info: Min_alignment_length: $minisize\n";} chomp; next unless (/^\d+/); @arr=split(/\s+/); print $arr[1], "\t", $arr[2], "\t", $arr[3], "\t", $arr[4], "\t", $arr[5], "\t", $arr[6], "\t", $arr[8], "\t", 0, "\n" if(($arr[3]-$arr[2])>=$minisize and ($arr[6]-$arr[5])>=$minisize);' ${finaloutput%.*}.ORIGIN.axt >> $finaloutput
perl -ne 'BEGIN{$minisize=$ENV{"MINISIZE"}; print STDERR "Info: Min_alignment_length: $minisize\n";} chomp; next unless (/^\d+/); @arr=split(/\s+/); next unless (($arr[3]-$arr[2])>=$minisize or ($arr[6]-$arr[5])>=$minisize); print $arr[1], "\t", $arr[2], "\t", $arr[3], "\t", $arr[4], "\t", $arr[5], "\t", $arr[6], "\t", $arr[8], "\t", $arr[7], "\n";' ${finaloutput%.*}.ORIGIN.axt >> $finaloutput


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



#1. Generate PSL alignments (e.g., with BLAT or lastz).
#2. Turn those alignments into chains with axtChain.
#3. Merge the short chains using chainMergeSort, chainSplit, and chainSort.
#4. You may wish to filter your chains at this point with chainPreNet, to remove chains that don't have a chance of being part of the final file.
#5. Create a net from the chains using the chainNet program, pass that to netSyntenic to add synteny information, use netChainSubset to create a liftOver file, and finally (optionally) join chain fragments with chainStitchId (this is skipped on the wiki page).
#
#The result is a liftOver chain file of the sort used to begin the reciprocal best pipeline.
#
#Here is an example segment from the scripts that created hg38.oviAri3.over.chain.gz. This corresponds to steps 4 and 5 above.
#
# Make nets ("noClass", i.e. without rmsk/class stats which are added later):
#chainPreNet  hg38.oviAri3.all.chain.gz /hive/data/genomes/hg38/chrom.sizes /hive/data/genomes/oviAri3/chrom.sizes stdout \
#| chainNet  stdin -minSpace=1 /hive/data/genomes/hg38/chrom.sizes /hive/data/genomes/oviAri3/chrom.sizes stdout /dev/null \
#| netSyntenic stdin noClass.net
#
# Make liftOver chains:
#netChainSubset -verbose=0 noClass.net hg38.oviAri3.all.chain.gz stdout \
#| chainStitchId stdin stdout | gzip -c > hg38.oviAri3.over.chain.gz
