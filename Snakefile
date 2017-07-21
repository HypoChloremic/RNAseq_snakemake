#################################################################################
#################################################################################
# Snakemake pipe for RNAseq
# Copyright (c) 2017 Ali Rassolie
# Karolinska Institutet
# V.0.0.8 GIMLET
# Versions:
## snakemake: 3.13.3
## cutadapt: 1.14
## multiqc: 1.0
## fastqc: 0.11.5
## salmon: 0.8.2
## conda: 4.3.21

# We are equating the files of interest in here

# Note that the script is GLOBAL. Further, we had previously discussed
# the incorporation of a csv file with the files of interest.
# We can either make a script that generates that.

# Here we import the scripts package.
from scripts.file_name_producer import FileExpander
file_init = FileExpander()

# Refers to the configfile where we got everything
configfile: "scripts/config.yaml"
PARENT_DIR    = config["PARENT_DIR"]
SALMON_INDEX  = config["SALMON_INDEX"]
TRANSCRIPTOME = config["TRANSCRIPTOME"]
URQT          = config["URQT_PATH"]
MULTIQC_REPORT_PATH = config["MULTIQC_REPORT_PATH"]
TASK_CPUS      = config["TASK_CPUS"]
SYMLINKED_PATH = config["SYMLINKED_DIR"]
# Refers to the scripts files, where we wish to get the paths
FULL_PATHS = file_init.SYMLINKED_PATHS
LOG_PATHS  = file_init.LOG_PATHS
# Files och NO_DATA_PATH är essentiella, år att de tillåter os
# separera filerna vi vill analysera och the folders in which they lie.
FILES        = file_init.FILES
NO_DATA_PATH = file_init.NO_DATA_PATH


# Jeff tells me that there may be different sequences for specific
# fastq files, and that the sequences can be found based upon some index
# from nextera?
ADAPTOR_SEQUENCE = config["ADAPTOR_SEQUENCE"]

#################################################################################
#################################################################################

#################
### The rules ###
#################
  input:
    # For cutadapt rule
    expand("{path}/data/cutadapt_output_1_{files}_fastq.gz", zip, path=NO_DATA_PATH, files=FILES),
    expand("{path}/data/cutadapt_output_2_{files}_fastq.gz", zip, path=NO_DATA_PATH, files=FILES),
    # For URQT rule
    expand("{path}/data/{files}URQT_output_1_fastq.gz", zip, path=NO_DATA_PATH, files=FILES),
    expand("{path}/data/{files}URQT_output_2_fastq.gz", zip, path=NO_DATA_PATH, files=FILES),
    # For fastq rule
    expand("{path}/log/{files}_untrim/", zip, path=NO_DATA_PATH, files = FILES),
    expand("{path}/log/{files}_trim/", zip, path=NO_DATA_PATH, files = FILES),
    # For salmon rule
    expand("{path}/data/{files}_transcripts_quant", zip, path=NO_DATA_PATH, files = FILES)


    first_file = "{path}/data/{files}_1.fastq.gz",
    second_file= "{path}/data/{files}_2.fastq.gz"

  output:
    first_file = "{path}/data/cutadapt_output_1_{files}_fastq.gz",
    second_file= "{path}/data/cutadapt_output_2_{files}_fastq.gz"

  log:
    "{path}/log/{files}_log.txt"

  shell:

rule urqt:
  input:
    first_file = "{path}/data/cutadapt_output_1_{files}_fastq.gz",
    second_file = "{path}/data/cutadapt_output_2_{files}_fastq.gz"

  output:
    urqt_trimmed_1 = "{path}/data/{files}URQT_output_1_fastq.gz",
    urqt_trimmed_2 = "{path}/data/{files}URQT_output_2_fastq.gz"

  log:
    "{path}/log/{files}_URQT_report.txt"

  shell:


rule fastqc_untrim:
  input:
    untrimmed_fastq_1 = "{path}/data/{files}_1.fastq.gz",
  output:
    "{path}/log/{files}_untrim/"

  shell:
    "fastqc --quiet --threads {TASK_CPUS} --outdir {output} {input.untrimmed_fastq_1} {input.untrimmed_fastq_2}"


rule fastqc_trim:
  input:
    trimmed_urqt_1 = "{path}/data/{files}URQT_output_1_fastq.gz",
    trimmed_urqt_2 = "{path}/data/{files}URQT_output_2_fastq.gz"
    "{path}/log/{files}_trim/"

  shell:
    "fastqc --quiet --threads {TASK_CPUS} --outdir {output} {input.trimmed_urqt_1} {input.trimmed_urqt_2}"

rule multiqc:
  shell:
    "multiqc {SYMLINKED_PATH}/ --filename  ./reports"


rule salmon:
  input:
    urqt_trimmed_1 = "{path}/data/{files}URQT_output_1_fastq.gz",
    urqt_trimmed_2 = "{path}/data/{files}URQT_output_2_fastq.gz"

  output:
    transcripts_quant = "{path}/data/{files}_transcripts_quant"

  shell:
    "salmon quant -i {SALMON_INDEX} -l IU -p 8 --useVBOpt --numBootstraps 100 --seqBias --gcBias --posBias -1 {input.urqt_trimmed_1} -2 {input.urqt_trimmed_2} -o {output.transcripts_quant}"