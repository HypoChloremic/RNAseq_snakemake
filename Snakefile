#################################################################################
#################################################################################
# Snakemake pipe for RNAseq
# (c) 2017 Ali Rassolie
# V.0.0.3 GIMLET
# Versions:
# snakemake: 3.13.3
# cutadapt: 1.14
# multiqc: 0.9
# fastqc: 0.11.5
# trimmomatic: 0.35
# conda: 4.3.21

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
INDEX    = config["INDEX"]
#################################################################################
#################################################################################
# Snakemake pipe for RNAseq
# (c) 2017 Ali Rassolie
# V.0.0.4 GIMLET
# Versions:
# snakemake: 3.13.3
# cutadapt: 1.14
# multiqc: 0.9
# fastqc: 0.11.5
# trimmomatic: 0.35
# conda: 4.3.21

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
INDEX    = config["INDEX"]
FULL_PATHS = file_init.FULL_PATHS
# Jeff tells me that there may be different sequences for specific
# fastq files, and that the sequences can be found based upon some index
# from nextera?
ADAPTOR_SEQUENCE = config["ADAPTOR_SEQUENCE"]

#################################################################################
#################################################################################

#################
### The rules ###
#################

rule all:
  input:
    expand("{something}_trimmed_adapt_1", something=FULL_PATHS),
    expand("{something}_clipped", something=FULL_PATHS)

rule cutadapt:
  input:
    "{something}_1.fastq.gz",
    "{something}_2.fastq.gz"

  output:
    "{something}_trimmed_adapt_1",
    "{something}_trimmed_adapt_2"

  shell:
    "cutadapt {ADAPTOR_SEQUENCE} -o  {output[0]} -p  {output[1]} {input[0]} {input[1]} > ${LOG}/cutadapt_report.txt"

rule urqt:

rule fastqc_untrim:

rule fastqc_trim:

rule salmon_index:

rule salmon_run:

rule ns:
  # this implies that we will be creating the entire path with a func.
  # pro of doing this is that it doesnt seem like the expand func can deal with not combining
  # list contents of different lists.

rule multiqc:

rule salmon_index:

rule salmon_run:

rule ns:
  # this implies that we will be creating the entire path with a func.
  # pro of doing this is that it doesnt seem like the expand func can deal with not combining
  # list contents of different lists.
  input:
    "{something}_1.fastq.gz",
    "{something}_2.fastq.gz"

  output:
    "{something}_clipped"

  shell:
    "nesoni clip --length 45 --quality 25 --out-separate yes {output} pairs: {input[0]} {input[1]}"

rule salmon:
  input:
    "{input}_clipped_R1.fq.gz"
    "{input}_clipped_R2.fq.gz"

  output:
    "{something}_salmoned"