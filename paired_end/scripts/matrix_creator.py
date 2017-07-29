# This will produce matrices from specified salmon quant outputs.
# (cc) 2017 Ali Rassolie
# Karolinska Institutet
# 
__version__ = "0.2.0 GIMLET"
__doc__ =""" 
This script has been adapted to be used with o

"""

import file_name_producer as fp
import re, os
import matplotlib.pyplot as plt
from pandas import DataFrame
from pandas import concat
from collections import Counter, OrderedDict
from time import time, sleep
from math import log10


class MatrixCreator(fp.FileExpander):

  def __init__(self, directory=None):

    _, SYMLINK_DIR, self.ANNOTATION_PATH = self.opener("config.yaml", term="ANNOTATION_PATH")
    files        = file_expander([directory], n=6)
    self.quants  = [ i for i in files if self.isquant(i) ]

  def annotation(self, **kwargs):

    with open(self.ANNOTATION_PATH) as file:
      #t = time()

      data = [ i.split() for i in file.readlines() ]
      #data = file.read().split()
      #t2 = time()
      #indices = [ i+(i2/(10**(int(log10(i2))+1))) for i, k in enumerate(data) for i2, n in enumerate(k) if "ENST" in n ]
      #indices = [ f"{i}.{i2}" for i, k in enumerate(data) for i2, n in enumerate(k) if "ENST" in n ]
  #    t3 = time()
 #     print(t2-t, t3-t2)

#      print(data[int(indices[0].split(".")[0])][int(indices[0].split(".")[1])])
      t = time()
      l = OrderedDict()
      self.temp_gene = None
      for each in data:
        for i, n in enumerate(each):
          try:
            if "gene_name" in n:
              l[n[i+1]] = None
              self.temp_gene = n[i+1]

            elif "ENST" in n:
              l[list(l.items())[-1][0]] = n
              


            else:
              continue
          except IndexError as e:
            pass 

      print(time()-t)
      print(l)

  def unfiltered_matrix(self, tmp="tmp", count="count"):
    tmp_name_unfilt = tmp
    count_name_unfilt = count
    textgen = self.text_processing(self.quants)

   # Initiating the generator.
    tmp, count, pstmp, pscount = next(textgen)  

   # Timing a future loop
    o = time()

   # Driving the generator. 
    for new_tmp, new_count, new_pstmp, new_pscount in textgen:
      tmp = concat([tmp, new_tmp], axis=1)
      count = concat([count, new_count], axis=1)
      pstmp = concat([pstmp, new_pstmp], axis=1)
      pscount = concat([pscount, new_pscount], axis=1)
      print(time()-o)
      o = time()


   # Filtering the matrices, in order to produce
   # the necessary files, where gene isoform values
   # have been added. 
   # The number of rows should equal to the number of genes expressed. 

   # Writing the data to files
    tmp.to_csv(f"{tmp_name_unfilt}_unfiltered.csv")
    count.to_csv(f"{count_name_unfilt}_unfiltered.csv")
    pstmp.to_csv(f"{tmp_name_unfilt}_ps.csv")
    pscount.to_csv(f"{count_name_unfilt}_ps.csv")
  


  def filtered_matrix_producer(self, tmp="tmp", count="count"):

    tmp_name_unfilt = tmp
    count_name_unfilt = count
    textgen = self.text_processing(self.quants)

    tmp, count, pstmp, pscount = next(textgen)

    for new_tmp, new_count, new_pstmp, new_pscount in textgen:
      tmp = concat([tmp, new_tmp], axis=1)
      count = concat([count, new_count], axis=1)
      pstmp = concat([pstmp, new_pstmp], axis=1)
      pscount = concat([pscount, new_pscount], axis=1)


   # Filtering the matrices, in order to produce
   # the necessary files, where gene isoform values
   # have been added.
   # The number of rows should equal to the number of genes expressed.
    filtered_tmp   = self.filtered_matrix(dataframe=tmp)
    filtered_count = self.filtered_matrix(dataframe=count)

    filtered_tmp.to_csv(f"{tmp_name_unfilt}_filtered.csv")
    filtered_count.to_csv(f"{count_name_unfilt}_filtered.csv")

  def filtered_matrix(self, dataframe=None, ncolumns=190):
    # The dataframe of interest. 
    # It will contain all the compounded values for each gene. 

    filtered_matrix = DataFrame( data    = [[0.]*len(dataframe.columns)]*len(self.gene_names),
                                 index   = self.gene_names, 
                                 columns = list(dataframe.columns))

    t = time()

    #for i in dataframe.index:
    print(len(dataframe.index.values))
    for i in range(len(dataframe.index.values)):
    # Create the necessary gene-names. 
      if "ERCC" in dataframe.index[i]:
        continue
      
    #  temp_name = 
      #print(dataframe.iloc[i].values)
    # Extract the gene-name from the dataframe index. 
      gene_name = dataframe.index[i].split("|")[-1]

    # Make sure that the gene has not been appended before we move on.       
      #print(filtered_matrix.loc[gene_name].values)
      #print(dataframe.iloc[i].values, dataframe.index[i])
      filtered_matrix.loc[gene_name] = filtered_matrix.loc[gene_name].add(dataframe.iloc[i])
      #print(filtered_matrix.loc[gene_name].values)
      print(time()-t)
      t = time()
      #sleep(0.1)


  def matrix_opener(self, i):
    with open(i, "r") as file:
      print(len(file.readlines()))

  def text_processing(self, i):
    self.name_count = Counter()
    self.annotation()
    tmp_creator   = lambda data, sample_tmp: DataFrame(data, columns=[sample_tmp], index=data["Name"])
    count_creator = lambda data, sample_count: DataFrame(data, columns=[sample_count], index=data["Name"])

    for pos,n in enumerate(i):
      print(pos)
#      if pos == 2:
#        return StopIteration

      with open(n, "r") as file:

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
              gene = self.an_dict[name.split(".")[0]]
              #name = f"{gene}"
              name = f"{name}|{gene}"
              gene_match["Name"].append(name)
              gene_match[sample_tmp].append(float(i[-2]))
              gene_match[sample_count].append(float(i[-1]))
           
              
            else:
              gene_match["Name"].append(name)
              gene_match[sample_tmp].append(i[-2])
              gene_match[sample_count].append(i[-1])
            
          except KeyError as e:
            
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
  # The regex param
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
  process = MatrixCreator(directory=)
  #process.unfiltered_matrix(tmp = "tmp3", count="count3")
  process.annotation()
  process.matrix_opener("tmp3_ps.csv")
