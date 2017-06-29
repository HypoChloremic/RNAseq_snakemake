from subprocess import call
import os, re

configfile: "./configfile.yaml"
# Parent directory
PARENT_DIR = config["PARENT_DIR"]
# Global path for ERCC used in salmon.
ERCC_IN = config["ERCC_PATH"]

print(config[helloworld()])

def isgz(i):
  bool = re.search(r"\w\.gz", i)

  try:
    assert bool
    return True

  except AssertionError:
    return False


# This lambda func will replace remove the suffix of the files that we find
clean = lambda i: [ re.findall(r"(\w+)_\d\.fastq\.gz", z) for z in i ] 

sample = PARENT_DIR
dir_list = os.listdir(sample)

# This is so as to present the folders in the current directory
first_dir_child = [ i for i in dir_list if os.path.isdir(i) ]

# Now we would like to explore every single folder that is in first_dir_child
full_paths = []
filenames = []

for folder in first_dir_child:

  if not folder == ".snakemake":
    temp_path = f"{sample}/{folder}"
    temp_list = os.listdir(temp_path)

    "We have an index 0 here, which may be problematic if there are several folders within the nest"
    temp_path = f"{sample}/{folder}/{temp_list[0]}"
    temp_list2 = os.listdir(f"{temp_path}") 
    act_paths = [ f"{temp_path}/{i}" for i in temp_list2 if isgz(i) ]
    temp_act_paths = clean(act_paths)
    full_paths.append(f"{temp_path}/{temp_act_paths[0][0]}")
  else:
    pass

# We are equating the files of interest in here

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

  shell: "nesoni clip --length 45 --quality 25 --out-separate yes {output} reads: {input[0]} {input[1]}"

rule salmon:
  shell: "salmon quant -i {ERCC_IN} -l IU -p 8 --useVBOpt --numBootstraps 100 --seqBias --gcBias --posBias -1 {something}_clipped_paired -o {something}_salmoned"

