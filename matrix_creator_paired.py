# This will produce matrices from specified salmon quant outputs.
# (cc) 2017 Ali Rassolie
# Karolinska Institutet
 
__version__ = "0.2.5 GIMLET"
__doc__ =""" 
This script has been adapted to be used with the following versions:
 python:    3.6.1
 snakemake: 3.13.3
 cutadapt:  1.14
 multiqc:   1.0
 fastqc:    0.11.5
 salmon:    0.8.2
 conda:     4.3.21
"""

# We will need to import several methods from fp, 
# among others the opener. 
import paired_end.jeff2.scripts.file_name_producer as fp
import paired_end.jeff2.scripts.ERROR_DEFS as ERROR_DEFS
import re, os
import matplotlib.pyplot as plt
from pandas import DataFrame
from pandas import concat
from collections import Counter, OrderedDict
from time import time, sleep
# log10 import has been depreciated. An implementation in 
# cython needs to be tested before reintroduction. 
from math import log10
from multiprocessing import Process, Queue
from subprocess import call


######################
### The Main Class ###
######################

class MatrixCreator(fp.FileExpander):

  def __init__(self, directory=None, 
               overwrite=False, 
               folder="temp", 
               config="matrix_config.yaml", 
               delimiter="/", 
               output_folder="matrix_output"):

    # The output folder name
    self.output_folder = output_folder

    # Importing error defs
    self.ERDEFS = ERROR_DEFS.Definitions()
   
    # Finding out what the name of the sample is, by splitting the directory
    self.file_name = directory.split(delimiter)[-1] 

    # Using the opener to get config paths.
    _, SYMLINK_DIR, self.ANNOTATION_PATH = self.opener(config, term="ANNOTATION_PATH")
    files        = file_expander([directory], n=6)
    self.quants  = [ i for i in files if self.isquant(i) ]

   
    # Will remove folder, so as to create a new annotation file that is written. 
    if overwrite is True:
      call(["rm", "-r", f"{folder}"])

  def annotation(self, file_name = None, **kwargs):
    if file_name:
      pass
    else:
      file_name = f"temp/temp_dict_{self.file_name}.dat"
    
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
        dat          = file.readlines()
        self.an_dict = dict()

        for each in dat:
          self.an_dict[each.split()[0]] = each.split()[1]

      # We create a set for the genes, so as to remove duplicates.
        self.gene_names = set(list(self.an_dict.values())) 
        print(f"The amount of rows in the gene dictionary is: {len(self.gene_names)}")
        print(f"Finished reading the gene-transcript dictionary for {time()-t}s")
        
    except FileNotFoundError as e:
      self.create_annotation()

  def create_annotation(self, **kwargs):
    # This methods will create the gene-transcript dictionary.

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
    print(f"The amount of rows in the gene dictionary is: {len(self.gene_names)}")
    print(f"Finished creating dictionary for: {time()-t}s")
    self.temp_creator(self.an_dict)

  def temp_creator(self, d, file_name=None, output_folder="matrix_output", folder="temp"):
    # temp_creator is a method that will create new directories.  
    if file_name:
      pass
    else:
      file_name = f"temp/temp_dict_{self.file_name}.dat"
    
    print(f"{self.ERDEFS.HEADER}Creating directory for {file_name}{self.ERDEFS.END}")
    call(["mkdir", f"{folder}", f"{output_folder}"])
    
    print(f"Creating {file_name}")
    with open(file_name, "w") as file:
      for key in iter(d):
        file.write(f"{key} {d[key]}\n")

  def run(self, func):
    pass       

  def run_matrix(self):
    self.annotation()
    print("{self.ERDEFS.HEADER}Running both FiltMatrix and UnfiltMatrix{self.ERDEFS.END}")
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
    textgen = self.text_processing(self.quants, parent="UNFILT")

   # Initiating the generator.
    tmp, count, pstmp, pscount = next(textgen)  

   # Driving the generator. 
    print(f"{self.ERDEFS.OKGREEN} (UNFILT) Driving unfilt generator...{self.ERDEFS.END}")
    for new_tmp, new_count, new_pstmp, new_pscount in textgen:

   # In this form, concat is used to add new columns
   # to new dataframes. Note that the dataframes are
   # large, causing time consumption. 
      tmp   = concat([tmp, new_tmp], axis=1)
      count = concat([count, new_count], axis=1)
      pstmp = concat([pstmp, new_pstmp], axis=1)
      pscount = concat([pscount, new_pscount], axis=1)
    print(f"{self.ERDEFS.OKGREEN} (UNFILT) Finished driving the unfilt generators{self.ERDEFS.END}")
   
   # Filtering the matrices, in order to produce
   # the necessary files, where gene isoform values
   # have been added. 
   # The number of rows should equal to the number of genes expressed. 

   # Writing the data to files. Note that these have not
   # been converted to tsv yet. 
    print(f"{self.ERDEFS.OKGREEN}(UNFILT) Saving the dataframes to csv...{self.ERDEFS.END}")

   # Here we are saving the dataframes to csv, by using a method
   # inherent to dataframes. The reason for not using a personal method
   # is to gain in performance, as time for optimisation is unavailable. 
    tmp.to_csv(f"{self.output_folder}/{tmp_name_unfilt}_unfiltered_{self.file_name}.csv")
    count.to_csv(f"{self.output_folder}/{count_name_unfilt}_unfiltered_{self.file_name}.csv")

    pstmp.to_csv(f"{self.output_folder}/{tmp_name_unfilt}_ps_{self.file_name}.csv")
    pscount.to_csv(f"{self.output_folder}/{count_name_unfilt}_ps_{self.file_name}.csv")
  
    print(f"{self.ERDEFS.OKGREEN}(UNFILT) Finished saving the dataframes to csv{self.ERDEFS.END}")

  def filtered_matrix(self, tmp="tmp", count="count"):

    tmp_name_unfilt   = tmp
    count_name_unfilt = count
    textgen = self.text_processing(self.quants, parent="FILT")

    tmp, count, pstmp, pscount = next(textgen)
    print(f"{self.ERDEFS.BLUE}{self.ERDEFS.BOLD}(FILT) Driving text-generator...{self.ERDEFS.END}")
    t = time()
    for new_tmp, new_count, new_pstmp, new_pscount in textgen:
      tmp   = concat([tmp, new_tmp],     axis=1)
      count = concat([count, new_count], axis=1)

    print(f"{self.ERDEFS.BLUE}{self.ERDEFS.BOLD}(FILT) Finished driving text-generator for {time()-t}s{self.ERDEFS.END}")

   # Filtering the matrices, in order to produce
   # the necessary files, where gene isoform values
   # have been added.
   # The number of rows should equal to the number of genes expressed.
   
   # to_run is a the queue we wish to run, so as to 
   # parse it into the processing loop. 
    to_run = ((tmp,"tmp"), (count, "count"))
    out_q = Queue()
    procs = list()

    for i in to_run:
     # Note that we are here using the to run
     # tuple that we are iterating through.
     # we are starting processes for each of these. The tricky part is to
     # solve the args issue. 
      p = Process(target=self.filtered_matrix_producer, kwargs={"dataframe": i[0], "out_q":out_q, "name": i[1]})
      p.start()
      print("(FILT) in the process loop")
    print("(FILT)") 
    p.join() 
    output = dict()
    for i in range(len(to_run)):
      data = out_q.get()
      output[data[0]] = data[-1]  
   
    output["tmp"].to_csv(f"{self.output_folder}/{tmp_name_unfilt}_filtered_{self.file_name}.csv")
    output["count"].to_csv(f"{self.output_folder}/{count_name_unfilt}_filtered_{self.file_name}.csv")

  def filtered_matrix_producer(self, **kwargs):
    dataframe = kwargs["dataframe"]
    out_q     = kwargs["out_q"]
    name      = kwargs["name"]
    ID        = id(dataframe) 
    # The dataframe of interest. 
    # It will contain all the compounded values for each gene. 
    print(f"(FILT) Producing filtered dataframes\nID: {ID}")

    # This is the dataframe file that will containt all of the
    # genes and the relevant counts. It is important to note that
    # we are creating this file, and do not append anythinh to it
    # so as to gain in performance and speed. 
    filtered_matrix = DataFrame(data    = [[0.]*len(dataframe.columns)]*len(self.gene_names),
                                index   = self.gene_names, 
                                columns = list(dataframe.columns))

    t = time()

    for i in range(len(dataframe.index.values)):
    # Create the necessary gene-names. 
      if "ERCC" in dataframe.index[i]:
        continue
      
    # Extract the gene-name from the dataframe by index. 
    # Note that we are iterating thriugh the range of the 
    # number of indices, or rows rather, implying we can access them by
    # index number.  
      gene_name = dataframe.index[i].split("|")[-1]
      filtered_matrix.loc[gene_name] = filtered_matrix.loc[gene_name].add(dataframe.iloc[i])

    output = [name, filtered_matrix]
    out_q.put(output)
    print(f"(FILT) Finished producing filtered matrix: {ID}\n")
    #return filtered_matrix



  # This method is used to read the matrices we have produced. 
  # Currently adapted to read from csv's and not tsv's. 
  def read_matrix(self, i):
    df = DataFrame.from_csv(i)
    return df

  # This generator is central to all of the producers.
  # It is used to process the salmon quant.sf files, 
  # performs dictionary lookups and appends the gene names to the ENST's.  
  def text_processing(self, i, parent=None):
    tmp_creator   = lambda data, sample_tmp: DataFrame(data, columns=[sample_tmp], index=data["Name"])
    count_creator = lambda data, sample_count: DataFrame(data, columns=[sample_count], index=data["Name"])
    
    for pos,n in enumerate(i):
      print(f"({parent}) {pos}")
      if pos == 2:
        return StopIteration

      with open(n, "r") as file:

      # This readline method call is problematic, for
      # input files with differing headers. 
        file.readlines(1)

      # Here we would like to save the name of the sample.
      # As the name of the file is always quant.sf, we would 
      # like to get a name that is unique for each sample
      # which can be found at index -2 when we split by "/"
        sample       = n.split("/")[-2]
        sample_tmp   = f"{sample}_tmp"
        sample_count = f"{sample}_count"

        gene_match    = {"Name": [], sample_tmp:[], sample_count:[]}
        no_gene_match = {"Name": [], sample_tmp:[], sample_count:[]}

      # We use k to save all of the rows in the file we are reading. 
      # o is used to measure the time taken for the run. 
        k = file.readlines()

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
        
      # Sorting them, such that we may get them in order, by virtue of their alphabetical 
      # order. 
        tmp     = tmp.sort_index()
        count   = count.sort_index()
        pstmp   = pstmp.sort_index()
        pscount = pscount.sort_index()

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

#####################
### Global Methods###
#####################

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
  process.run_prf()
