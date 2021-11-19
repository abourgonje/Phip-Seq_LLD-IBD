from pathlib import Path
from subprocess import call

Path_results = "/groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Plink_GWAS/Meta_analysis/Results/"

#agilent_32775
#twist_69666


for File in Path(Path_results).glob("*.tbl"):
	#§:File = Path("/groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Plink_GWAS/Meta_analysis/Results/agilent_139646.txt1.tbl")
	Name = File.name.split(".")[0]
	temp_variants= "tmp/Significant_SNPs_{Name}".format(Name=Name)
	temp_subset= "tmp/subset_{Name}".format(Name=Name)
	temp_variants2= "tmp/Clumped_SNPs_{Name}".format(Name= Name)
	VCF= "Genotypes/Clumped_SNPs_{Name}".format(Name= Name)
	Annotation= "Annotation/VEP_{Name}".format(Name = Name)
	if Path(Annotation).exists(): continue

	for Folder in ["tmp", "Genotypes", "Annotation","logs"]:
		if not Path(Folder).exists():
			command = "mkdir "+ Folder
			print(command)
			call(command, shell = True)
		else: pass
	Command ="""#!/bin/sh
#SBATCH  --job-name={N}.job
#SBATCH --output=logs/{N}.out
#SBATCH --error=logs/{N}.err
#SBATCH --time=0:10:00
#SBATCH --mem-per-cpu=6G
#SBATCH --nodes=1

ml VEP/104.2-foss-2018b-Perl-5.28.0
ml PLINK/1.9-beta6-20190617

#vep -i {VCF}.vcf  --database --grch37 --output_file {Annotation} --force_overwrite
#exit

echo "Using awk to extract variants 5e-8"
#1. Identify significant and suggestive variants
awk '{{if ($4<5e-8) print $1}}' {File}  > {temp_variants}
N=$(wc -l {temp_variants}) ; N=$(echo "$N" | awk -F' ' '{{print $1}}')

echo "Number of variants extracted $N"
if [[ $N -lt 1 ]] ; then
	exit
fi

echo "Extracting subset of variants from genotype file"
#2.  Extract
plink --bfile /groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Genetic_Preprocess/2_Relatedness_ethnicity/LLD --extract {temp_variants} --make-bed --out {temp_subset}

#2.5 Get unique IDs
cut -f 2 tmp/subset_{N}.bim | sort | uniq -d > tmp/{N}.dups

echo "Clumping SNPs based on Pvalues (R2: 0.1)"
#3. Get independent SNPs
plink --bfile {temp_subset} --exclude tmp/{N}.dups --clump {File} --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.1 --clump-kb 1000 --out Lead/{N}_leadSNPs_clumps --clump-snp-field MarkerName --clump-field Weight

#4. Extract independent SNPs
echo Extract independent SNPs on {temp_variants2}
if [ ! -f Lead/{N}_leadSNPs_clumps.clumped ]; then
	echo "Clumping failed, quitting"
	exit
fi

awk -F' {{1,}}' '{{print $4}}' Lead/{N}_leadSNPs_clumps.clumped | tail -n+2 > {temp_variants2}
plink --bfile /groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Genetic_Preprocess/2_Relatedness_ethnicity/LLD --extract {temp_variants2} --recode vcf --out {VCF}
echo Annotation of SNPs
vep -i {VCF}.vcf  --database --grch37 --output_file {Annotation} --force_overwrite

rm {temp_variants} {temp_subset} {temp_variants2} 

""".format(N=Name, temp_variants=temp_variants, temp_subset=temp_subset, temp_variants2=temp_variants2, VCF=VCF, Annotation=Annotation, File=File)
	
	Script = "run_"+Name
	with open(Script, "w") as F:
		F.write(Command)
	Execute = "sbatch "+ Script
	#exit(Script)
	call(Execute, shell=True)
	Execute2 = "rm "+ Script
	call(Execute2, shell=True)




