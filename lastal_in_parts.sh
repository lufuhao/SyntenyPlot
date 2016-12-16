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

$0 --- Align large genome to reference by spliting into parts and merge back

Version: 20161216

Requirements:
	perl && File::Spec
	samtools
	maf-convert
	lastal

Descriptions:
	Split a long query sequence into shorter length (-l)
	And do the lastal alignment
	For genome synteny purpose
	
	**NOTE: -g option is set to avoid the last part is too small
	For example: seq length 1000199 and you set -l 10000 -v 200
	So the last file is 990001-1000199 to prevent the 199 short seq

Options:
  -h    Print this help message
  -i    Query sequence(s) in fasta
  -l    Step length for each lastal run, default: 10000000
  -g    Overhang length, default: 1/10 of step_length (-l)
  -p    lastal options
  -x    index of lastdb
  -o    Final output in PSL format
  -a    Start position for interrupted jobs, default: 1

Example:
  lastal_in_parts.sh -x $rundir/$lastdb_index -l 10000000 -i $queryseq \
      -o chr3Bbac.last.chr3Bnegene.psl -a 10000001\
      -p "-f MAF -Q 0 -p $scoremetrix -P $threads -M" 

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
echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
fastainput=''
steplen=10000000
lastaloption=''
lastalindex=''
finalout=''
overhang=0;
startpos=1;
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) fastainput=$2;shift 2;;
    -l) steplen=$2;shift 2;;
    -p) lastaloption=$2;shift 2;;
    -x) lastalindex=$2;shift 2;;
    -o) finalout=$2; shift 2;;
    -g) overhang=$2; shift 2;;
    -a) startpos=$2; shift 2;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 100;;
    *) break;;
  esac
done


#################### Subfuctions ####################################
###Detect command existence
CmdExists () {
  if command -v $1 >/dev/null 2>&1; then
    echo 0
  else
    echo 1
  fi
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
RunLastal () {
	local RLquery=$1
	local RLoutput=$2
	local RLsubinfo="SH(RunLastal)"
	
	lastal $lastaloption $lastalindex $RLquery > $RLquery.maf
	if [ $? -ne 0 ]; then
		echo "${RLsubinfo}Error: lastal running failed" >&2
		echo "CMD used: lastal $lastaloption $lastalindex $RLquery > $RLquery.maf" >&2
		return 1;
	elif [ ! -s "$RLquery.maf" ]; then
		echo "${RLsubinfo}Error: lastal output failed" >&2
		echo "CMD used: lastal $lastaloption $lastalindex $RLquery > $RLquery.maf" >&2
		return 1;
	fi
	maf-convert psl $RLquery.maf > $RLoutput
	if [ $? -ne 0 ]; then
		echo "${RLsubinfo}Error: maf2psl running failed" >&2
		echo "CMD used: maf-convert psl $RLquery.maf > $RLoutput" >&2
		return 1;
	elif [ ! -s "$RLoutput" ]; then
		echo "${RLsubinfo}Error: maf2psl output failed" >&2
		echo "CMD used: maf-convert psl $RLquery.maf > $RLoutput" >&2
		return 1;
	fi
	return 0;
}
#abs2rel () { perl -MFile::Spec -e 'print(File::Spec->abs2rel($ARGV[1], $ARGV[0]), "\n")' "$@"; }
#rel2abs () { perl -MFile::Spec -e 'print(File::Spec->rel2abs($ARGV[0]), "\n")' "$@"; }

#################### Command test ###################################
if [ $(CmdExists 'samtools') -ne 0 ]; then
	echo "Error: CMD 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'lastal') -ne 0 ]; then
	echo "Error: CMD 'lastal' in PROGRAM 'LAST' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################
start=1;
if [ $overhang -eq 0 ]; then
	((overhang=steplen/10));
	overhang=$(printf "%.f" $overhang)
fi
if [[ ! $overhang =~ ^[0-9]+$ ]]; then
	echo "Error: invalid overhang, must be INT number" >&2
	exit 100;
fi



#################### Input and Output ###############################
if [ -z "$fastainput" ] || [ ! -s "$fastainput" ]; then
	echo "Error: invalid query sequence" >&2
	exit 100;
fi
if [ -e "$finalout" ]; then
	rm -rf "$finalout" > /dev/null 2>&1
fi
	



#################### Main ###########################################
echo "##### SUMMARY #####"
echo "Fasta input:    $fastainput"
echo "Step: length:   $steplen"
echo "LAST option:    $lastaloption"
echo "LAST index:     $lastalindex"
echo "Output:         $finalout"
echo "Overhang:       $overhang"
echo "##### SUMMARY ends #####"

if [ ! -s "$fastainput.fai" ]; then
	samtools faidx "$fastainput"
	if [ $? -ne 0 ];then
		echo "Error: index query sequence failed" >&2
		echo "CMD used: samtools faidx $fastainput"
		exit 100;
	elif [ ! -s "$fastainput.fai" ]; then
		echo "Error: index query output failed" >&2
		echo "CMD used: samtools faidx $fastainput"
		exit 100;
	fi
fi
declare -a seqnamelen=($(perl -lane 'print "$F[0]::::$F[1]";' "$fastainput.fai"))
for indseq in "${seqnamelen[@]}"; do
	echo "Seq: ${indseq%%::::*} length ${indseq##*::::}";
	i=1;
	rm query*.fa > /dev/null 2>&1
	for ((start=$startpos; start<=${indseq##*::::}; start+=steplen)); do
		let end=$start+$steplen-1
		if [[ $end -gt ${indseq##*::::} ]]; then
			end=${indseq##*::::}
		elif [[ $((${indseq##*::::}-$overhang)) -le $end ]];then
			end=${indseq##*::::}
		fi
		echo "${indseq%%::::*}:$start-$end"
		samtools faidx "$fastainput" ${indseq%%::::*}:$start-$end > query$start-$end.fa
		if [ $? -ne 0 ]; then
			echo "Error: subseq extract failed: ${indseq%%::::*}:$start-$end" >&2
			echo "CMD used: samtools faidx $fastainput ${indseq%%::::*}:$start-$end > query$start-$end.fa" >&2
			exit 100;
		elif [ ! -s "query$start-$end.fa" ]; then
			echo "Error: subseq output failed: ${indseq%%::::*}:$start-$end" >&2
			echo "CMD used: samtools faidx $fastainput ${indseq%%::::*}:$start-$end > query$start-$end.fa" >&2
			exit 100;
		fi
		if `RunLastal query$start-$end.fa query$start-$end.psl`; then
			echo "Info: lastal running succeeds: ${indseq%%::::*}:$start-$end"
		else
			echo "Info: lastal running failed: ${indseq%%::::*}:$start-$end" >&2
			exit 100;
		fi
		export MAFSTARTNUM=$(($start-1))
		perl -ne 'BEGIN {$addnum=$ENV{"MAFSTARTNUM"};print STDERR "Info: Addnumber=$addnum\n"; $linenum=0;} if (/^#/) {print;} chomp; @arr=split(/\t/); die "Error: invalid PSL at line$linenum: $_\n" unless (scalar(@arr)==21); if ($arr[9]=~/:(\d+)-\d+$/) {die "Error: invalid addnumber $addnum, should be $1\n" unless ($addnum == ($1-1)); $arr[9]=~s/:\d+-\d+$//;}else {die "Error: invalid sequence name \n";} $arr[11]+=$addnum; $arr[12]+=$addnum; if ($arr[8] eq "+") {@arr2=split(/,/, $arr[19]); for ($i=0;$i<scalar(@arr2);$i++){$arr2[$i]+=$addnum if ($arr2[$i]=~/^\d+$/);} $arr[19]=join(",", @arr2)}; print join("\t", @arr), "\n";' query$start-$end.psl >> $finalout
		if [[ $end -eq ${indseq##*::::} ]]; then
			break
		fi
	done
	rm query*.fa > /dev/null 2>&1
done


exit 0;
