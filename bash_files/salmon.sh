 # add the parameters path of the fastq files and the path of the index 

cd $1
echo $PWD
echo "~/scratch/genome/${2}"

for i in *.gz
do
        echo $i         
        base="${i}_quant"
        salmon quant -i "~/scratch/genome/${2}"  -l A -r $i -o $base -p 16 --numBootstraps 30  --useVBOpt -$
done

