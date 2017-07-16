# Snakemake pipe for RNAseq
# (c) 2017 Ali Rassolie
__version__ = "0.0.6 GIMLET"


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

import os, re

class FileExpander:
  def __init__(self):

# This is the config file, necessary effective functioning.
# Please do note that we are doing this relative to the Snakefile folder, and not the scripts
# folder. Some strange python thing.
    configfile = "scripts/config.yaml"
    clean = lambda i: [ re.findall(r"(\w+)_\d\.fastq\.gz", z) for z in i if z]

# Here we open the configfile to access the parent directory; we should maybe use JSON formatting,
# snakemake might not like it alas.
    with open(configfile, "r") as file:
      text_data = file.read().replace('"','').replace("'","").split("\n")
      text_data = { i.split(": ")[0]:i.split(": ")[1] for i in text_data if i}
      PARENT_DIR = text_data["PARENT_DIR"]
      input(f"parent dir {PARENT_DIR}")
# Now we are defining the global variable FULL_PATHS

# Please do make sure that the input PARENT_DIR is a list.
    self.FULL_PATHS = file_expander([PARENT_DIR])
    print("\n".join(self.FULL_PATHS))
    input()

# This lambda func will replace remove the suffix of the files that we find
    self.FULL_PATHS = clean(self.FULL_PATHS)
    print(self.FULL_PATHS)



def isgz(i):
  param = re.search(r"\w\.gz$", i)
  try:
    assert param
    return True
  except AssertionError:
    return False


def file_expander(parent, n=3, delimiter="/"):
  for i in range(n):
    if i == n-1:
      files = []
      for parent_path in parent:

        files.extend([f"{parent_path}{delimiter}{k}" for k in os.listdir(parent_path) if not os.path.isdir(f"{parent_path}{delimiter}{k}") and isgz(f"{parent_path}{delimiter}{k}")])
      return files
    else:
      parent = expand(parent, delimiter)


def expand(f, delimiter="/"):
  l = []
  for each in f:
    files_in_d = [ k for k in os.listdir(each) if os.path.isdir(f"{each}{delimiter}{k}")]

    for i in files_in_d:
      l.append(f"{each}{delimiter}{i}")
  return l