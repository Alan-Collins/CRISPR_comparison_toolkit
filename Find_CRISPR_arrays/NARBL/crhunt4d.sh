#!/bin/bash
#Usage: crhunt INFILE DRFILE OUTDIR CPU
#infile = fasta of reads
#drfile = fasta of known repeats to BLAST
#outdir = output directory
#cpu = # threads to use

#variables
infile=$1
drfile=$2
outdir=$3
maxcpu=$4
infile_base=$(basename $1)
drfile_base=$(basename $2)

#create output directories
if [ ! -d $outdir ]; then
mkdir $outdir
fi
mkdir $outdir/dbfiles

#split fasta for smaller DBs
numlines=$(wc -l $infile | awk '{print $1}')
lpf=$(expr $numlines / $maxcpu)
lpf_rd=$(printf %.0f $lpf)
if [ $((lpf%2)) -eq 0 ];
then
    #echo "even"
lpf=$((lpf+2))

else
    #echo "odd";
lpf=$((lpf+3))
fi
split -d -l $lpf $infile $outdir/$infile_base.tmpfasta_

#make job completion checkfile
rm -f $outdir/dbfiles/checkif.jobdone
touch $outdir/dbfiles/checkif.jobdone

#make blast DB of each split read file
for file in $outdir/*tmpfasta_* ; do
fb=$(basename $file)
(makeblastdb -max_file_sz 1GB -in $file -dbtype nucl -parse_seqids -out $outdir/dbfiles/$fb.db && echo "$fb done!" >> $outdir/dbfiles/checkif.jobdone &)
done

jd=0
while [ $jd -lt $maxcpu ]; do
sleep 10
jd=$(cat $outdir/dbfiles/checkif.jobdone | wc -l)
#echo $jd
done
rm -f $outdir/*tmpfasta_*
echo "BLAST database construction complete"

#get dblist
rm -f $outdir/*.nal
dblist=$(ls -1 $outdir/dbfiles/*.nhr | tr '\n' ' ' | sed 's/\.nhr//g')
blastdb_aliastool -dblist "$dblist" -dbtype nucl -out $outdir/dbfiles/$infile_base.all -title $infile_base.all

#blast all repeats
blastn -task blastn-short -query $drfile -db $outdir/dbfiles/$infile_base.all -num_threads $maxcpu -max_target_seqs 100000000 -evalue 0.1 -outfmt '6 std sstrand' -out $outdir/allrepeats_vs_"$infile_base"
echo "BLAST complete"

#sort out good hits, split by query, pull associated reads
querylist=$(grep ">" $drfile | sed 's/>//')
for q in $querylist
do
echo $q
awk -v query=^$q '{if ($1 ~ query && $3 >= 90) print $2,$13}' $outdir/allrepeats_vs_"$infile_base" | sort -u > $outdir/$q.list
qlines=$(wc -l < $outdir/$q.list)
if [ $qlines -gt 0 ]
then
blastdbcmd -db $outdir/dbfiles/$infile_base.all -dbtype nucl -line_length 500 -entry_batch $outdir/$q.list -out $outdir/$q.reads
fi
done
echo "Reads of interest acquired."
echo "Cleaning temp files"
rm -fr $outdir/dbfiles $outdir/*.list #$outdir/allrepeats_vs_"$infile_base"
echo "Mission accomplished"
