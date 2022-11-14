#BSUB -J greedy_selection_v2
#BSUB -P acc_ipm2
#BSUB -q premium
#BSUB -n 1
#BSUB -R "span[hosts=1] affinity[core(10, same=socket, exclusive=(socket, injob))
#BSUB -R rusage[mem=33000]
#BSUB -W 24:00
#BSUB -o 270922_greedy_selection_v2.out
#BSUB -e 270922_greedy_selection_v2.err
#BSUB -L /bin/bash
ml plink
ml anaconda3

raw=/hpc/users/shemir03/roohy/hprc/1kg/combined_samples_plink
for(( iteration = 1 ; iteration <= 25 ; iteration++))
do

    for(( i = 1 ; i <= 22 ; i++))
    do
        plink --bfile ${raw}/chr${i} --make-bed --out ./temp/sample_chr${i} --keep ./temp/ids.txt --mac 1 --memory 32000
        plink --bfile ${raw}/chr${i} --out ./temp/oos_chr${i} --keep ./temp/oos_ids.txt --maf 0.01 --recode vcf-iid --memory 32000
    done
    python generate_ids_greedy.py
done