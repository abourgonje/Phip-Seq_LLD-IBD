from pathlib import Path
from subprocess import call

Path_LLD = "/groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Plink_GWAS/GWAS_low_priority"
Path_IBD = "/groups/umcg-gastrocol/tmp01/Shixian/Antibody_GWAS/Summary.statistics.txt"

All_LLD = "Input/LLD_GWAS.tsv"
All_IBD = "Input/IBD_GWAS.tsv"
Redo = False
if Redo == True:
	with open(All_LLD, "w") as F:
		F.write("\t".join(["CHR","SNP", "BP", "A1", "TEST", "NMISS", "OR", "STAT", "P", "probe", "CHR:BP"])+"\n")
	#/groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Plink_GWAS/GWAS_low_priority/twist_28862.assoc.logistic
	for F in Path(Path_LLD).glob("*.assoc.logistic"):
		if Redo == False: break
		Name = F.name.split(".ass")[0]
		AWK_command = """awk -v v={v} '(NR>1) {{print $0,"\\t",v,"\\t",$1":"$3}}' {FILE} """.format(v=Name, FILE=str(F))
		AWK_command += " >> "+All_LLD
		call(AWK_command, shell=True)


Common_probes = []
Not_common = []
IBD_probes = []
LLD_probes = []
for Probe in Path("Split_IBD/").glob("*.txt"):
	IBD_probes.append(Probe.name)
for Probe in Path("Split_LLD/").glob("*.txt"):	
	P = Probe.name
	LLD_probes.append(Probe.name)
	if P in IBD_probes: Common_probes.append(P)
	else: Not_common.append(P)

Not_common.extend(list(set(IBD_probes) - set(Common_probes)) )

print("Length IBD: {N} \nLength LLD: {N2} \nLength not matching: {N3}".format(N= len(IBD_probes), N2 = len(LLD_probes), N3= len(list(Not_common)) ))
print("In IBD but not in LLD: {N}".format(N= len( list(set(IBD_probes) - set(LLD_probes)) ) ) )
print("In LLD but not in IBD: {N}".format(N= len( list(set(LLD_probes) - set(IBD_probes)) ) ) )



#F = """ 
#SEPARATOR WHITESPACE
#MARKER probe:CHR:BP
#ALLELE1 A1
#EFFECT log(OR)
#PVALUE P
#WEIGHT NMISS
#PROCESS Input/LLD_GWAS_SNP2.tsv 
## === THE SECOND INPUT FILE HAS THE SAME FORMAT AND CAN BE PROCESSED IMMEDIATELY ===
#PROCESS Input/IBD_GWAS_SNP2.tsv 

#OUTFILE Results/All .tbl
#ANALYZE
#QUIT
#"""
#SCRIPT = "scripts/All"
#with open(SCRIPT, "w") as O: O.write(F)
#F2 = """#!/bin/sh
##BATCH  --job-name=metal.job
##SBATCH --time=24:00:00
##SBATCH --mem-per-cpu=30G
##SBATCH --nodes=1
#ml Metal/2011-03-25-foss-2018b
#metal SOURCE {S} 
#""".format(S=SCRIPT)
#SCRIPT2 = "scripts/All.sh"
#with open(SCRIPT2, "w") as O: O.write(F2)
#exit(SCRIPT2)
#call("sbatch " + SCRIPT2, shell=True)

#exit()

for P in Common_probes:
	if Path("Results/" + P + "1.tbl.info").exists(): continue
	print(P)
	F = """
SEPARATOR WHITESPACE
MARKER CHR:BP
ALLELE1 A1
EFFECT log(OR)
PVALUE P
WEIGHT NMISS
PROCESS Split_LLD/{N}
# === THE SECOND INPUT FILE HAS THE SAME FORMAT AND CAN BE PROCESSED IMMEDIATELY ===
PROCESS Split_IBD/{N}

OUTFILE Results/{N} .tbl
ANALYZE
QUIT
""".format(N=P)
	SCRIPT = "scripts/{P}".format(P=P)
	with open(SCRIPT, "w") as O: O.write(F)
	F2 = """#!/bin/sh
#SBATCH  --job-name=metal.job
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
ml Metal/2011-03-25-foss-2018b
metal SOURCE {S}
""".format(S=SCRIPT)
	SCRIPT2 = "scripts/{P}.sh".format(P=P.split(".")[0])
	with open(SCRIPT2, "w") as O: O.write(F2)
	#aexit(SCRIPT2)
	call("sbatch " + SCRIPT2, shell=True)

	
