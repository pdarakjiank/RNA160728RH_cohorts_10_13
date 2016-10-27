for cnt in `seq 1 24`
do
        echo "Working on sample ${cnt}"
        fastqc RNA160728RH*_S${cnt}_*gz 
        echo "Finished fastqc"
done
