# Snakemake pipe for RNAseq
# (c) 2017 Ali Rassolie
__version__ = "0.1.0 VAGRANT"


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

class FileExpander:
  def __init__(self):

# This is the config file, necessary effective functioning.
# Please do note that we are doing this relative to the Snakefile folder, and not the scripts
# folder. Some strange python thing.
    configfile = "scripts/config.yaml"
    clean = lambda i: [ re.findall(r"(.*)_\d\.fastq\.gz", z)[0] for z in i if re.findall(r"(.*)_\d\.fastq\.gz", z)]

# Here we open the configfile to access the parent directory; we should maybe use JSON formatting,
# snakemake might not like it alas.
    try:
      PARENT_DIR, SYMLINKED_DIR = self.opener(configfile)
    except FileNotFoundError as e:
      PARENT_DIR, SYMLINKED_DIR = self.opener("config.yaml")
# Now we are defining the global variable FULL_PATHS
# The variable is important to produce the file-names
# the rules require

# Please do make sure that the input PARENT_DIR is a list.
    self.FULL_PATHS = file_expander([PARENT_DIR])
    self.SYMLINKED_PATHS = symlink_creator(origin=self.FULL_PATHS, target=SYMLINKED_DIR)
    print(self.SYMLINKED_PATHS)

# This lambda func will replace remove the suffix of the files that we find
    self.FULL_PATHS = clean(self.FULL_PATHS)
    self.SYMLINKED_PATHS = clean(self.SYMLINKED_PATHS)
    print("symlinked")
    print("\n".join(self.SYMLINKED_PATHS))
    input()


  def opener(self, i):
    with open(i, "r") as file:
      text_data = file.read().replace('"','').replace("'","").split("\n")
      text_data = { i.split(": ")[0]:i.split(": ")[1] for i in text_data if i}
      PARENT_DIR = text_data["PARENT_DIR"]
      SYMLINKED_DIR = text_data["SYMLINKED_DIR"]
      return PARENT_DIR, SYMLINKED_DIR