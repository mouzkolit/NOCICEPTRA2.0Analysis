# this is used to trim randomly the first 4 bases from the beginning and the end
# this is necessary as indicated in the NextFlex documentation

for i in *trimmed*
do
	cutadapt -u 4 -u -4 -o "random_${i}" $i
done
