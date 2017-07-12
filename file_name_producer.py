# Snakemake pipe for RNAseq
# (c) 2017 Ali Rassolie

print(
"""
*********************************
*********************************
WARNING!
The following script is global.
Press RETURN to continue the run.
To terminate run, press ctrl+Z
*********************************
*********************************

""")
input()


from subprocess import call
import os, re

class Constants:
  pass


configfile = "config.yaml"



with open(configfile, "r") as file:
  text_data = file.read().replace('"','').replace("'","").split("\n")
  text_data = [ i.split(": ") for i in text_data ]
  input (text_data)
  PARENT_DIR = txt_convert("PARENT_DIR", text_data)
  INDEX = txt_convert("INDEX", text_data)


  PARENT_DIR = config["PARENT_DIR"]
  INDEX = config["INDEX"]
  

def txt_convert(i, txt):
  pass

def isgz(i):
  bool = re.search(r"\w\.gz", i)
  try:
    assert bool
    return True
  except AssertionError:
    return False



def file_expander(n=3):
  pass



# This lambda func will replace remove the suffix of the files that we find
clean = lambda i: [ re.findall(r"(\w+)_\d\.fastq\.gz", z) for z in i if z]

sample = PARENT_DIR
dir_list = os.listdir(sample)

print(dir_list)

# This is so as to present the folders in the current directory
first_dir_child = [ i for i in dir_list if os.path.isdir(f"{sample}/{i}") ]

# Now we would like to explore every single folder that is in first_dir_child
full_paths = []
filenames = []

print(first_dir_child)

for folder in first_dir_child:
  if not folder == ".snakemake":
    temp_path = f"{sample}/{folder}"
    temp_list = os.listdir(temp_path)

   # "We have an index 0 here, which may be problematic if there are several folders within the nest"
    temp_path = f"{sample}/{folder}/{temp_list[0]}"
    print(temp_path)
    temp_list2 = os.listdir(f"{temp_path}")
    act_paths = [ f"{temp_path}/{i}" for i in temp_list2 if isgz(i) ]
    temp_act_paths = clean(act_paths)
    print(f"temp {temp_act_paths}")

    for index in range(len(temp_act_paths)):
      try:
        full_paths.append(f"{temp_path}/{temp_act_paths[index][0]}")
        break
      except IndexError as e:
        pass
  else:
    pass

 
