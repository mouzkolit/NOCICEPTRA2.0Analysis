for i in *gz
do
	#trims the selected adapter sequence
	cutadapt -a TGGAATTCTCGGGTGCCAAGG -o "trimmed_${i}" --cores 8  --minimum-length 23 $i
done
