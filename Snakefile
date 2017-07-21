#################################################################################
#################################################################################
# Snakemake pipe for RNAseq
# (c) 2017 Ali Rassolie
# V.0.0.7 GIMLET
# Versions:
## snakemake: 3.13.3
## cutadapt: 1.14
## multiqc: 0.9
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
PARENT_DIR = config["PARENT_DIR"]
SALMON_INDEX    = config["SALMON_INDEX"]
URQT = config["URQT_PATH"]
MULTIQC_REPORT_PATH = config["MULTIQC_REPORT_PATH"]
TASK_CPUS = config["TASK_CPUS"]

# Refers to the scripts files, where we wish to get the paths
FULL_PATHS = file_init.SYMLINKED_PATHS
LOG_PATHS = file_init.LOG_PATHS

# Files och NO_DATA_PATH är essentiella, år att de tillåter os
# separera filerna vi vill analysera och the folders in which they lie.
FILES = file_init.FILES
NO_DATA_PATH = file_init.NO_DATA_PATH


# Jeff tells me that there may be different sequences for specific
# fastq files, and that the sequences can be found based upon some index
# from nextera?
#################################################################################
#################################################################################

#################
### The rules ###
#################

rule all:
  input:
    expand("{path}/data/{files}_1_trimmed_out", zip, path=NO_DATA_PATH, files = FILES),
    expand("{path}/data/{files}_2_trimmed_out", zip, path=NO_DATA_PATH, files = FILES),
    expand("{path}/log/{files}_fastqc_log", zip, path=NO_DATA_PATH, files = FILES)

  input:
    first_file = "{path}/data/{files}_1.fastq.gz",
    second_file= "{path}/data/{files}_2.fastq.gz"

  output:
    first_file = "{path}/data/cutadapt_output_1_{files}",
    second_file= "{path}/data/cutadapt_output_2_{files}"

  log:
    "{path}/log"

  shell:

rule urqt:
  input:
    first_file = "{path}/data/{files}_1.fastq.gz",
    second_file = "{path}/data/{files}_2.fastq.gz"

  output:
    "{path}/data/{files}URQT_output"

  log:
    "{path}/log/{files}_URQT_report.txt"

  shell:


rule fastqc_untrim:
  input:
    "{path}/data/{files}_1.fastq.gz",
    "{path}/data/{files}_2.fastq.gz"

  output:
    "{path}/log/{files}_fastqc_log"

  shell:
    "fastqc {input} --quiet --threads {TASK_CPUS} --outdir {output}"


rule fastqc_trim:
  input:
    "{path}/data/{files}_URQT_output"

  output:
    "{path}/log/{files}_fastqc_log"
    "fastqc {input} --quiet --threads {TASK_CPUS} --outdir {output}"

rule multiqc:
  input:
    expand("{path}/log/", path=NO_DATA_PATH)

  output:
    "{MULTIQC_REPORT_PATH}/"

  shell:
    "multiqc {input} --filename --o {output}"


rule salmon_index:
  shell:
    "salmon index -t {SALMON_INDEX}"


rule salmon:
  input:
    "{path}/data/{files}_clipped_R1.fq.gz",
    "{path}/data/{files}_clipped_R2.fq.gz"

  output:
    "{path}/data/{files}_salmoned"

  shell:
    "salmon quant -i {SALMNON_INDEX} -l IU -p 8 --useVBOpt --numBootstraps 100 --seqBias --gcBias --posBias -1 {input[0]} -2 {input[0]} -o {output}"
