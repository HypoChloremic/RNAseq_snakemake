# Snakemake pipe for RNAseq
# (c) 2017 Ali Rassolie
__version__ = "0.0.9 GIMLET"


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
    except FileNotFoundError as e:
      PARENT_DIR, SYMLINKED_DIR = self.opener("config.yaml")

# Now we are defining the global variable FULL_PATHS
# The variable is important to produce the file-names
# the rules require

# Please do make sure that the input PARENT_DIR is a list.
    self.FULL_PATHS = file_expander([PARENT_DIR])
    print(self.FULL_PATHS)
    input()
    self.SYMLINKED_PATHS = symlink_creator(origin=self.FULL_PATHS, target=SYMLINKED_DIR)

# This lambda func will replace remove the suffix of the files that we find
    self.FULL_PATHS = clean(self.FULL_PATHS)
    print("\n".join(self.FULL_PATHS))

  def opener(self, i):
    with open(i, "r") as file:
      text_data = file.read().replace('"','').replace("'","").split("\n")
      text_data = { i.split(": ")[0]:i.split(": ")[1] for i in text_data if i}
      PARENT_DIR = text_data["PARENT_DIR"]
      SYMLINKED_DIR = text_data["SYMLINKED_DIR"]
      input(f"parent dir {PARENT_DIR}")
      return PARENT_DIR, SYMLINKED_DIR

# Please do note that the origin and the target must be a list; perhaps a dictionary
  # This is the complete list of all the files that the parent_dir contained.
  # symlink_creator is a method that will create symbolic links between the origins and the targets.
  origins = origin
  target = target
  file_names = [ "/".join(i.split("/")[-n:]) for i in origins]
  target_paths = [f"{target}/{i}" for i in file_names]
  for each_target in zip(origin, target_paths):
      call(["ln", "-s", each_target[0], each_target[1])


def isgz(i):
  param = re.search(r"\w\.gz$", i)
  try:
    assert param
    return True
  except AssertionError:
    return False

# In this case we have the following file-tree structure:
# .../Experiment/sample/run/... and we thus provide n=3
# so as to expand the third dir, and return its contents.
def file_expander(parent, n=3, delimiter="/"):
  for i in range(n):
    if i == n-1:
      files = []
      for parent_path in parent:

        files.extend([f"{parent_path}{delimiter}{k}" for k in os.listdir(parent_path) if not os.path.isdir(f"{parent_path}{delimiter}{k}") and isgz(f"{parent_path}{delimiter}{k}")])
      return files
    else:
      parent = expand(parent, delimiter)

# expand() method will take an array containing paths, and return a list where the
# paths have been opened and returned in a list.
def expand(f, delimiter="/"):
  l = []
  for each in f:
    # This list compr will expand the path of spec index in the input path list
    # and it will be appending children that are directories themselves.
    files_in_d = [ k for k in os.listdir(each) if os.path.isdir(f"{each}{delimiter}{k}")]

    for i in files_in_d:
      l.append(f"{each}{delimiter}{i}")
  return l

