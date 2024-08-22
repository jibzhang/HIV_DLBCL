#Each Group/Core is limited to 216 cores and 4.6TB in defq, and 192 Cores, 2.2TB and 24 V100 GPU in gpu
configfile: 'config_sample_HIVDLBCL.yaml'

import pandas as pd

#sample list
intake=pd.read_table(config["intake"])
samplelist=list(intake.sample_id.drop_duplicates())
dir_out=config["dir_out"]

rule all:
    input:
        expand(dir_out+"/smk_code/{sample}.sh",sample=samplelist)

rule write_code:
    input:
        key_smk="HIV_DLBCL.smk",
    output:
        code=dir_out+"/smk_code/{sample}.sh"
    params:
        cores=4,
        mem="60G",
        log=dir_out+"/smk_log/{sample}_sv_%j.log",
        key_file="/coh_labs/jochan/HIV_DLBCL/{sample}/EXOME/log/done_sv.txt"
    shell:
        '''
        echo "#!/bin/bash" > {output.code}
        echo "#SBATCH --job-name=sv_{wildcards.sample}           # Job name" >> {output.code}
        echo "#SBATCH -n {params.cores}                          # Number of cores" >> {output.code}
        echo "#SBATCH -N 1-1                        # Min - Max Nodes" >> {output.code}
        echo "#SBATCH -p compute                    # gpu queue" >> {output.code}
        echo "#SBATCH --mem={params.mem}                     # Amount of memory in GB" >> {output.code}
        echo "#SBATCH --time=8:00:00               # Time limit hrs:min:sec" >> {output.code}
        echo "#SBATCH --output={params.log}    # Standard output and error log" >> {output.code}
        echo "snakemake -s {input.key_smk} -p {params.key_file}  -j{params.cores}"  >> {output.code}
        #sbatch {output.code}
        #scontrol update jobid=179841 partition=gpu-scavenge
        '''

















