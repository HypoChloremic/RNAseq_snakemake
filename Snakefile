# Snakemake pipe for RNAseq
# (c) 2017 Ali Rassolie

# Following will be including script import for the pipe
import file_name_producer

rule all:
  input: expand("{something}_clipped", something=full_paths)

rule ns:
  # this implies that we will be creating the entire path with a func. 
  # pro of doing this is that it doesnt seem like the expand func can deal with not combining
  # list contents of different lists.
  input:
    "{something}_1.fastq.gz",
    "{something}_2.fastq.gz"

  output:
    "{something}_clipped"
  shell: "nesoni clip --length 45 --quality 25 --out-separate yes {output} pairs: {input[0]} {input[1]}"

rule salmon:
  shell: "salmon quant -i {INDEX} -l IU -p 8 --useVBOpt --numBootstraps 100 --seqBias --gcBias --posBias -1 {something}_clipped_R1.fq.gz -2 {something}_clipped_R2.fq.gz -o {something}_salmoned"

