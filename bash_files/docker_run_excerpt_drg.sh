for i in  $PWD/Human_DRG/*.gz
do
	sample="$(basename -- $i)"
	echo $sample

	docker run -v $PWD/Human_DRG:/exceRptInput \
           -v $PWD/ExcerptDRG:/exceRptOutput \
           -v $PWD/hg38:/exceRpt_DB/hg38 \
           -t rkitchen/excerpt \
        INPUT_FILE_PATH=/exceRptInput/$sample \
        N_THREADS=8 \
        ADAPTER_SEQ='guessKnown' 
	MIN_READ_LENGTH=17

done
