

for F in Phenotypes/*  ; do
	N=$(basename $F)
	sbatch Run_association.sh $N GWAS $F
	
	
done


q()



for F in Phenotypes_2/*  ; do
        N=$(basename $F)
        #echo $F

        sbatch Run_association.sh $N GWAS_low_priority $F

done
