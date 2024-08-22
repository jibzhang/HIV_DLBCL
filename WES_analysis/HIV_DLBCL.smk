configfile: 'config_sample_HIVDLBCL.yaml'
configfile: 'config_software.yaml'

include: "common_rule.smk"
include: "GenerateAnalysisReadyBam.smk"
include: "Somatic_mutation.smk"
include: "copy_number_nopair_analysis.smk"
include: "CNVkit.smk"
include: "structural_variation.smk"

ALL_cosmic = expand(dir_out + "/{sample}/" + squence_type + "/log/{sample}.search_cosmic.bmk", sample=samplelist)
#ALL_dbsnp = expand(dir_out + "/{sample}/" + squence_type + "/log/{sample}.dbsnp_gnomad.bmk", sample=samplelist_tumor)
ALL_mutation = expand(dir_out + "/{sample}/" + squence_type + "/log/{sample}.filter_public.bmk", sample=samplelist_tumor)
ALL_CNV = expand(dir_out + "/{sample}/" + squence_type + "/log/doneCNVkit.txt", sample=samplelist)
ALL_SV = expand(dir_out + "/{sample}/" + squence_type + "/log/done_sv.txt", sample=samplelist)

TASKS = []
TASKS.extend(ALL_cosmic)
#TASKS.extend(ALL_dbsnp)
TASKS.extend(ALL_mutation)
TASKS.extend(ALL_CNV)
TASKS.extend(ALL_SV)

rule all:
    input: TASKS
    output:
        dir_out + "/{sample}/" + squence_type + "/log/done_analysis.txt"
    shell:
        '''
        touch {output}
        '''