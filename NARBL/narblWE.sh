#!/bin/bash
#usage: narbl.sh readfile_in repfile_in
#expectations: readfile is named <repeatid>.reads
#requires: fuzznuc, blast+, fastx

#variables
readin=$1
readin_base=$(basename $1)
repin=$2

#kill any duplicated reads (by read name)
mv $readin $readin.orig
cat $readin.orig | paste - - | sort -u -k1,1 | tr '\t' '\n' > $readin

#this will make fuzznuc happy
sed 's/:/_/g' $readin > $readin_base.fix

#acquire repeat
export repname=$(basename $1 .reads)
repeat=$(perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw('$repname')}print if $c' $repin | tail -n +2)
echo "repeat is: $repeat"
#make fuzznuc pattern - repeat w/ 12 bases on either side
pattern=$(echo "n(12)""$repeat""n(12)")

#hunt down pattern in reads
fuzznuc -sequence $readin_base.fix -pattern $pattern -pmismatch 8 -complement yes -outfile $repname.fuzz -rformat_outfile excel

#remove header columns
grep -v "SeqName	Start	End	Score	Strand	Pattern	Mismatch" $repname.fuzz > $repname.fuzz.headless

#eliminate bidirectional matches
#find dupes
#awk '{print $1}' $repname.fuzz.headless | sort | uniq -d > $repname.bid
#save best (fewest mismatches) hit for each.  eliminate everything else.
#rm -f $repname.bid.tmp
#if [ -s $repname.bid ]
#then
#while read line
#do
#grep $line $repname.fuzz.headless | sort -g -k7 | tail -n +2 >> $repname.bid.tmp
#done < $repname.bid
#grep -v -f $repname.bid.tmp $repname.fuzz.headless > $repname.fuzz.clean
#else
cat $repname.fuzz.headless > $repname.fuzz.clean
#fi

#get read lengths
rm -f $repname.rdlen
names=$(awk '{print $1}' $repname.fuzz.clean)
for n in $names
do
grep -A 1 "$n" $readin_base.fix | tail -n +2 | awk '{print length($1)}' >> $repname.rdlen
done

#get reads and positions, flip around positions if on wrong strand
paste $repname.fuzz.clean $repname.rdlen | awk '{if ($5 ~ "+") {print $1, $2, $3} else {print $1,$8-$3+1,$8-$2+1}}' | sort -k 1,1 -k 2,2n > $repname.flankpos
rm $repname.rdlen

#take names of reads that need to be RCed.
#this should really be unnecessary
awk '{if ($5 ~ "-") {print $1}}' $repname.fuzz.clean | sort -u > $repname.rcthese

#split into fasta to RC and fasta to not RC
rm -f rc.tmp.fasta norc.tmp.fasta
namesu=$(awk '{print $1}' $repname.fuzz.clean | sort -u)
for nu in $namesu
do
if [ $(grep -c "$nu" $repname.rcthese) -eq 1 ]
then
grep -A 1 "$nu" $readin_base.fix >> rc.tmp.fasta
else 
grep -A 1 "$nu" $readin_base.fix >> norc.tmp.fasta
fi
done
#wc -l rc.tmp.fasta
#wc -l norc.tmp.fasta
#RC those who need it
if [ -s rc.tmp.fasta ]
then
echo "Danger Will Robinson! RCing underway!"
fastx_reverse_complement -i rc.tmp.fasta -o rced.tmp.fasta
fi

#put humpty dumpty together again
if [ -a rced.tmp.fasta -a -a norc.tmp.fasta ]
then
echo "both"
cat rced.tmp.fasta norc.tmp.fasta > $readin_base.fullrep.fa
rm rced.tmp.fasta norc.tmp.fasta rc.tmp.fasta
elif [ -a rced.tmp.fasta ]
then
echo "all RC"
mv rced.tmp.fasta $readin_base.fullrep.fa
rm rc.tmp.fasta
elif [ -a norc.tmp.fasta ]
then
echo "no RC"
mv norc.tmp.fasta $readin_base.fullrep.fa
fi

#blastdb of full-repeat reads
makeblastdb -dbtype nucl -in $readin_base.fullrep.fa -parse_seqids -out $readin_base.frdb
#get the actual repeat sequence
awk '{print $1,$2+12"-"$3-12}' $repname.flankpos > $repname.repeatpos
blastdbcmd -db $readin_base.frdb -entry_batch $repname.repeatpos -out $repname.repeats.fasta
grep -v ">" $repname.repeats.fasta | sort -u > $repname.realreps

#pull up/down chunks, remove duplicates
awk '{print $1,$2"-"$2+11"\n"$1,$3-11"-"$3}' $repname.flankpos > $repname.chunkpos
blastdbcmd -db $readin_base.frdb -entry_batch $repname.chunkpos -out $repname.chunks
#to track L/R status
awk '{print "L\nR"}' $repname.flankpos > $repname.chunkside
grep -v ">" $repname.chunks | paste - $repname.chunkside | sort | uniq -c > $repname.chunkcounts
awk '{if ($1 !~ "^1$") {print ">"$2"\n"$2}}' $repname.chunkcounts > $repname.chunks.mult
awk '{if ($1 ~ "^1$") {print ">"$2"\n"$2}}' $repname.chunkcounts > $repname.chunks.single
makeblastdb -dbtype nucl -in $repname.chunks.mult -out $repname.mchunkdb
#problem: singletons which overlap chunks.  When you go look for matches in the full read set, both come up - splits the graph.
#eliminate singletons with significant overlap to other chunks
#blast singletons against non-singletons
#blastn-short but longer word size
blastn -word_size 8 -reward 1 -query $repname.chunks.single -db $repname.mchunkdb -evalue 0.1 -outfmt 6 -out $repname.schunk.hits
awk '{print $1}' $repname.schunk.hits | sort -u > $repname.schunk.uhits
grep -v -f $repname.schunk.uhits $repname.chunkcounts | awk '{print $2,$3}'> $repname.chunks.uniq.full
awk '{print $1}' $repname.chunks.uniq.full > $repname.chunks.uniq

#eliminate chunks that are a perfect match to the actual repeats
fuzznuc -auto -sequence $repname.realreps -pattern @$repname.chunks.uniq -pmismatch 0 -complement no -outfile $repname.chunks.toss -rformat_outfile excel 
if [ -s $repname.chunks.toss ]
then
echo "noooooooooo"
grep -v "SeqName	Start	End	Score	Strand	Pattern	Mismatch" $repname.chunks.toss | awk 'BEGIN {FS = "[:\t]"} ; {print $7}' | sort -u > $repname.chunks.toss.uniq
rm -f $repname.chunks.toss.tmp
while read line
do
grep $line $repname.chunks.uniq.full >> $repname.chunks.toss.tmp
done < $repname.chunks.toss.uniq
comm -23 $repname.chunks.uniq.full $repname.chunks.toss.tmp > $repname.chunks.uniq.clean
awk '{print $1}' $repname.chunks.uniq.clean > $repname.chunks.uniq.clean.pat
else
echo "wooooo"
cat $repname.chunks.uniq > $repname.chunks.uniq.clean.pat
cat $repname.chunks.uniq.full > $repname.chunks.uniq.clean
fi

#locate chunks on reads
fuzznuc -sequence $readin_base.fix -pattern @$repname.chunks.uniq.clean.pat -pmismatch 0 -complement no -outfile $repname.chunks.uniq.clean.fuzz -rformat_outfile excel

grep -v "SeqName	Start	End	Score	Strand	Pattern	Mismatch" $repname.chunks.uniq.clean.fuzz | tee $repname.chunks.uniq.clean.fuzz.nh | awk 'BEGIN {FS = "[:\t]"} ; {print $7}' > $repname.chunks.uniq.clean.fuzz.tmp
rm -f $repname.chunks.uniq.clean.flank.tmp
while read line
do
grep $line $repname.chunks.uniq.clean | awk '{print $2}' >> $repname.chunks.uniq.clean.flank.tmp
done < $repname.chunks.uniq.clean.fuzz.tmp
paste $repname.chunks.uniq.clean.fuzz.nh $repname.chunks.uniq.clean.flank.tmp | tee $repname.chunks.final | awk '{print $1}' | sort -u > $repname.int

rm -f $repname.links.all
rm -f $repname.spacerpos.all
rm -f $repname.reppos.all
while read line
do
grep "$line" $repname.chunks.final | awk 'BEGIN {FS = "[:\t]"} ; {print $7,$9,$2,$3}' | sort -u > $repname.chunk.tmp
#chunkarray=( $( cat $repname.chunk.tmp ) )
chunkarray=( $(awk '{print $1}' $repname.chunk.tmp) )
flankarray=( $(awk '{print $2}' $repname.chunk.tmp) )
startarray=( $(awk '{print $3}' $repname.chunk.tmp) )
stoparray=( $(awk '{print $4}' $repname.chunk.tmp) )

for ((num=0;num<${#chunkarray[@]}-1;num++))
do 
	for((numplus=$((num + 1));numplus<=${#chunkarray[@]}-1;numplus++))
	do
		if [ ${startarray[$num]} -lt ${startarray[$numplus]} ]
		then
			#if RL link, get spacer position
			if [ ${flankarray[$num]} == R -a ${flankarray[$numplus]} == L ]
			then
			#echo  "$((${startarray[$num]}-${stoparray[$numplus]}))"
				if [ $((${startarray[$num]}-${stoparray[$numplus]})) -lt -60 ]
				then
				echo "${chunkarray[$num]}	${flankarray[$num]}	${chunkarray[$numplus]}	${flankarray[$numplus]}" >> $repname.links.long
				echo "$line	${startarray[$num]}-${stoparray[$numplus]}" >> $repname.spacerpos.long
				else
				echo "${chunkarray[$num]}	${flankarray[$num]}	${chunkarray[$numplus]}	${flankarray[$numplus]}" >> $repname.links.short
				echo "$line	${startarray[$num]}-${stoparray[$numplus]}" >> $repname.spacerpos.short
				fi
			#if LR link, get repeat position
			elif [ ${flankarray[$num]} == L -a ${flankarray[$numplus]} == R ]
			then
			#echo "$(($((${stoparray[$num]}+1))-$((${startarray[$numplus]}-1))))"
				if [ $(($((${stoparray[$num]}+1))-$((${startarray[$numplus]}-1)))) -lt -60 ]
				then
				echo "${chunkarray[$num]}	${flankarray[$num]}	${chunkarray[$numplus]}	${flankarray[$numplus]}" >> $repname.links.long
				echo "$line	$((${stoparray[$num]}+1))-$((${startarray[$numplus]}-1))" >> $repname.reppos.long
				else
				echo "${chunkarray[$num]}	${flankarray[$num]}	${chunkarray[$numplus]}	${flankarray[$numplus]}" >> $repname.links.short
				echo "$line	$((${stoparray[$num]}+1))-$((${startarray[$numplus]}-1))" >> $repname.reppos.short
				fi
			fi
		elif [ ${startarray[$num]} -gt ${startarray[$numplus]} ]
		then
			if [ ${flankarray[$numplus]} == R -a ${flankarray[$num]} == L ]
			then
			#echo "$((${startarray[$numplus]}-${stoparray[$num]}))"
				if [ $((${startarray[$numplus]}-${stoparray[$num]})) -lt -60 ]
				then
				echo "${chunkarray[$numplus]}	${flankarray[$numplus]}	${chunkarray[$num]}	${flankarray[$num]}" >> $repname.links.long
				echo "$line	${startarray[$numplus]}-${stoparray[$num]}" >> $repname.spacerpos.long
				else
				echo "${chunkarray[$numplus]}	${flankarray[$numplus]}	${chunkarray[$num]}	${flankarray[$num]}" >> $repname.links.short
				echo "$line	${startarray[$numplus]}-${stoparray[$num]}" >> $repname.spacerpos.short
				fi
			elif [ ${flankarray[$numplus]} == L -a ${flankarray[$num]} == R ]
			then
			#echo "$(($((${stoparray[$numplus]}+1))-$((${startarray[$num]}-1))))"
				if [ $(($((${stoparray[$numplus]}+1))-$((${startarray[$num]}-1)))) -lt -60 ]
				then
				echo "${chunkarray[$numplus]}	${flankarray[$numplus]}	${chunkarray[$num]}	${flankarray[$num]}" >> $repname.links.long
				echo "$line	$((${stoparray[$numplus]}+1))-$((${startarray[$num]}-1))" >> $repname.reppos.long
				else
				echo "${chunkarray[$numplus]}	${flankarray[$numplus]}	${chunkarray[$num]}	${flankarray[$num]}" >> $repname.links.short
				echo "$line	$((${stoparray[$numplus]}+1))-$((${startarray[$num]}-1))" >> $repname.reppos.short
				fi
			fi
		else
		echo "Oh no! Something is wrong with ${flankarray[$num]} ${flankarray[$numplus]}"
		fi
		done
	done
done < $repname.int
#go get spacers & repeats
makeblastdb -dbtype nucl -in $readin_base.fix -parse_seqids -out $readin_base.fulldb
blastdbcmd -line_length 500 -db $readin_base.fulldb -entry_batch $repname.spacerpos.short -out $repname.spacers.final.fasta
blastdbcmd -line_length 500 -db $readin_base.fulldb -entry_batch $repname.reppos.short -out $repname.repeats.final.fasta
#final links table
sort $repname.links.short | uniq -c | awk 'BEGIN {OFS = "\t"}{print $1,$2,$4,$3$5}' > $repname.links.count.txt
sort $repname.links.long | uniq -c | awk 'BEGIN {OFS = "\t"}{print $1,$2,$4,$3$5"X"}' >> $repname.links.count.txt
#| grep 'RL\|LR' 

#spacer-pulling commands for convenience
awk -v repname=$repname '{print $0, "grep -E '\''"$2".*"$3"'\'' "repname".spacers.final.fasta"}' $repname.links.count.txt > $repname.spacerpull

