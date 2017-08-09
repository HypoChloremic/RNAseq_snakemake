# This will produce matrices from specified salmon quant outputs.
# (cc) 2017 Ali Rassolie
# Karolinska Institutet
 
__version__ = "0.2.0 GIMLET"
__doc__ =""" 
This script has been adapted to be used with the following versions:
 snakemake: 3.13.3
 cutadapt: 1.14
 multiqc: 1.0
 fastqc: 0.11.5
 salmon: 0.8.2
 conda: 4.3.21
"""

import file_name_producer as fp
import re, os
import matplotlib.pyplot as plt
from pandas import DataFrame
from pandas import concat
from collections import Counter, OrderedDict
from time import time, sleep
from math import log10
from multiprocessing import Process
from subprocess import call

class MatrixCreator(fp.FileExpander):

  def __init__(self, directory=None, overwrite=False, folder="temp"):

    _, SYMLINK_DIR, self.ANNOTATION_PATH = self.opener("config.yaml", term="ANNOTATION_PATH")
    files        = file_expander([directory], n=6)
    self.quants  = [ i for i in files if self.isquant(i) ]
    if overwrite is True:
      call(["rm", "-r", f"{folder}"])

  def annotation(self, file_name = "temp/temp_dict.dat", **kwargs):
    # This method is going to look for tmp folders, where 
    # a dictionary file exists.
    print("\nAnnotation")
    try:
      with open(file_name, "r") as file: 
        print(f"Opening {file_name}")
        t   = time()

      # temp_creator have included newlines. This allows us 
      # to use readlines to produce the necessary entries in
      # the list.       
        dat = file.readlines()
        self.an_dict = dict()
        for each in dat:
          self.an_dict[each.split()[0]] = each.split()[1]

      # We create a set for the genes, so as to remove duplicates.
        self.gene_names = set(list(self.an_dict.values())) 
        print(f"Finished reading the gene-transcript dictionary for {time()-t}s")
        
    except FileNotFoundError as e:
      self.create_annotation()

  def create_annotation(self, **kwargs):
    print("Creating gene-transcript dictionary...")
    t = time()
    with open(self.ANNOTATION_PATH) as file:
      dat  = file.read().replace('"', "").replace(";","").split("\n")
      data = [ i.split() for i in dat ]

      self.an_dict   = dict()
      self.temp_gene = None
      self.annotation_index_errors = 0
     
    for each in data:
      bool_ = False
      for i, n in enumerate(each):          
          
        if "ENST" in n:
          self.transcript = n
          self.an_dict[self.transcript] = None
          bool_ = True

        elif "gene_name" in n and bool_ is True:
          self.an_dict[self.transcript] = each[i+1]

        elif "gene_name" in n and bool_ is False:
          try:
            self.an_dict[f"None_{i}"] = each[i+1]

          except IndexError as e:
            self.annotation_index_errors += 1

        else:
          pass

    self.gene_names = set(list(self.an_dict.values()))
    print(f"Finished creating dictionary for: {time()-t}s")
    self.temp_creator(self.an_dict)

  def temp_creator(self, d, file_name= "temp/temp_dict.dat", folder="temp"):

    print(f"Creating directory for {file_name}")
    call(["mkdir", f"{folder}"])

    print(f"Creating {file_name}")
    with open(file_name, "w") as file:
      for key in iter(d):
        file.write(f"{key} {d[key]}\n")

  def run(self, func):
    pass       

  def run_matrix(self):
    self.annotation()
    print("Running both FiltMatrix and UnfiltMatrix")
    run_filt   = Process(target=self.filtered_matrix)
    run_unfilt = Process(target=self.unfiltered_matrix)
    run_filt.start()
    run_unfilt.start()

  def run_prf(self):
    self.annotation()
    run = Process(target=self.filtered_matrix)
    run.start()

  def run_pruf(self):
    self.annotation()
    run = Process(target=self.unfiltered_matrix)
    run.start()

  def unfiltered_matrix(self, tmp="tmp", count="count"):
    tmp_name_unfilt   = tmp
    count_name_unfilt = count
    textgen = self.text_processing(self.quants)

   # Initiating the generator.
    tmp, count, pstmp, pscount = next(textgen)  
#    func_calls = (concat([tmp, new_tmp], axis=1),
#                  concat([count, new_count], axis=1),
#                  concat([pstmp, new_pstmp], axis=1),
#                  concat([pscount, new_pscount], axis=1))
   # Driving the generator. 
    print("Driving unfilt generator...")
    for new_tmp, new_count, new_pstmp, new_pscount in textgen:

   # In this form, concat is used to add new columns
   # to new dataframes. Note that the dataframes are
   # large, causing time consumption. 
      tmp   = concat([tmp, new_tmp], axis=1)
      count = concat([count, new_count], axis=1)
      pstmp = concat([pstmp, new_pstmp], axis=1)
      pscount = concat([pscount, new_pscount], axis=1)
    print("Finished driving the unfilt generators")
   
   # Filtering the matrices, in order to produce
   # the necessary files, where gene isoform values
   # have been added. 
   # The number of rows should equal to the number of genes expressed. 

   # Writing the data to files. Note that these have not
   # been converted to tsv yet. 
    print("Saving the dataframes to csv...")
    tmp.to_csv(f"{tmp_name_unfilt}_unfiltered.csv")
    count.to_csv(f"{count_name_unfilt}_unfiltered.csv")

    pstmp.to_csv(f"{tmp_name_unfilt}_ps.csv")
    pscount.to_csv(f"{count_name_unfilt}_ps.csv")
  

  def filtered_matrix(self, tmp="tmp", count="count"):

    tmp_name_unfilt   = tmp
    count_name_unfilt = count
    textgen = self.text_processing(self.quants)

    tmp, count, pstmp, pscount = next(textgen)
    print("Driving text-generator...")
    t = time()
    for new_tmp, new_count, new_pstmp, new_pscount in textgen:
      tmp   = concat([tmp, new_tmp],     axis=1)
      count = concat([count, new_count], axis=1)
      #pstmp = concat([pstmp, new_pstmp], axis=1)
      #pscount = concat([pscount, new_pscount], axis=1)

    print(f"Finished driving text-generator for {time()-t}s")

   # Filtering the matrices, in order to produce
   # the necessary files, where gene isoform values
   # have been added.
   # The number of rows should equal to the number of genes expressed.
    filtered_tmp   = self.filtered_matrix_producer(dataframe=tmp)
    filtered_count = self.filtered_matrix_producer(dataframe=count)

    filtered_tmp.to_csv(f"{tmp_name_unfilt}_filtered.csv")
    filtered_count.to_csv(f"{count_name_unfilt}_filtered.csv")

  def filtered_matrix_producer(self, dataframe=None, ncolumns=190):
    # The dataframe of interest. 
    # It will contain all the compounded values for each gene. 
    print(f"Producing filtered dataframes\nID: {id(dataframe)}")
    filtered_matrix = DataFrame(data    = [[0.]*len(dataframe.columns)]*len(self.gene_names),
                                index   = self.gene_names, 
                                columns = list(dataframe.columns))

    t = time()

    for i in range(len(dataframe.index.values)):
    # Create the necessary gene-names. 
      if "ERCC" in dataframe.index[i]:
        continue
      
    # Extract the gene-name from the dataframe index. 
      gene_name = dataframe.index[i].split("|")[-1]
      filtered_matrix.loc[gene_name] = filtered_matrix.loc[gene_name].add(dataframe.iloc[i])
    print(f"Finished producing filtered matrix: ")
    return filtered_matrix

  def open_matrix(self, i):
    df = DataFrame.from_csv(i)
    return df

  def text_processing(self, i):
    #self.name_count = Counter()
    tmp_creator   = lambda data, sample_tmp: DataFrame(data, columns=[sample_tmp], index=data["Name"])
    count_creator = lambda data, sample_count: DataFrame(data, columns=[sample_count], index=data["Name"])

    for pos,n in enumerate(i):
      print(pos)
#      if pos == 2:
#        return StopIteration

      with open(n, "r") as file:
      # This readline method call is problematic, for
      # input files with differing headers.
 
        file.readlines(1)
        sample = n.split("/")[-2]
        sample_tmp = f"{sample}_tmp"
        sample_count = f"{sample}_count"

        gene_match    = {"Name": [], sample_tmp:[], sample_count:[]}
        no_gene_match = {"Name": [], sample_tmp:[], sample_count:[]}

        k = file.readlines()
        o = time()

        for i in k:
          i    = i.replace("\n", "").split("\t")
          name = i[0].split("|")[0]
          try:
            if "ERCC" not in name:
              gene = self.an_dict[name]
              name = f"{name}|{gene}"
              gene_match["Name"].append(name)
              gene_match[sample_tmp].append(float(i[-2]))
              gene_match[sample_count].append(float(i[-1]))
              
            else:
              gene_match["Name"].append(name)
              gene_match[sample_tmp].append(i[-2])
              gene_match[sample_count].append(i[-1])
            
          except KeyError as e:
          # Receiving a KeyError here implies that 
          # the ENST transcript code found no match
          # in the transcript-gene dictionary.  
            print("No math")
            sleep(1) 
            no_gene_match["Name"].append(name)
            no_gene_match[sample_tmp].append(i[-2])
            no_gene_match[sample_count].append(i[-1])

      # Here we create the dataframes, where we specify what columns, or keys rather
      # because we in this case are dealing with dictionaries and not dataframe, to use. 
        
        tmp     = tmp_creator(gene_match, sample_tmp)
        count   = count_creator(gene_match, sample_count)
        pstmp   = tmp_creator(no_gene_match, sample_tmp)
        pscount = count_creator(no_gene_match, sample_count)

        yield tmp, count, pstmp, pscount


  # Returns files with quant prefix. 
  # Boolean function. 

  def isquant(self, i):
  # The regex param. It is going to be used to
  # identify the correct data files of interest. 
    param = re.search(r"\w/quant\.sf$", i)

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

        files.extend([f"{parent_path}{delimiter}{k}" for k in os.listdir(parent_path) if not os.path.isdir(f"{parent_path}{delimiter}{k}") ])
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

if __name__ == "__main__":
  process = MatrixCreator(directory="/media/box2/Experiments/Jeff/RNAseq/paired_end/jeff2/output/P1316_Female_A2_Data/")
  #process.run_matrix()
  process.run_pruf()
