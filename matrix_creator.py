# This will produce matrices from specified salmon quant outputs.
# (cc) 2017 Ali Rassolie
# Karolinska Institutet
 
__version__ = "0.2.5 GIMLET"
__doc__ =""" 
This script has been adapted to be used with the following versions:
 snakemake: 3.13.3
 cutadapt: 1.14
 multiqc: 1.0
 fastqc: 0.11.5
 salmon: 0.8.2
 conda: 4.3.21
"""

# We will need to import several methods from fp, 
# among others the opener. 
import paired_end.jeff2.scripts.file_name_producer as fp
import paired_end.jeff2.scripts.ERROR_DEFS as ERROR_DEFS
import re, os
import matplotlib.pyplot as plt
from pandas import DataFrame
from pandas import concat
from collections import Counter 
from collections import OrderedDict
from collections import deque
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

    if directory:
      # The output folder name
      self.output_folder = output_folder

      # Importing error defs
      self.ERDEFS = ERROR_DEFS.Definitions()
   
      # Finding out what the name of the sample is, by splitting the directory
      self.file_name = directory.split(delimiter)[-1] 

      # Using the opener to get config paths.
      _, SYMLINK_DIR, self.ANNOTATION_PATH = self.opener(config, term="ANNOTATION_PATH")
      files        = expander(directory)
      #self.quants  = [ i for i in files if self.isquant(i) ]
      self.quants = files
    else:
      pass
 
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
    run_filt.join()
    run_unfilt.join()

  def run_prf(self):
    self.annotation()
    print(f"{self.ERDEFS.OKGREEN}")
    run = Process(target=self.filtered_matrix)
    run.start()
    run.join()

  def run_pruf(self):
    self.annotation()
    print(f"{self.ERDEFS.BLUE}")
    run = Process(target=self.unfiltered_matrix)
    run.start()
    run.join()

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

   # There seems to be an issue with regard to the necessity of 
   # flushing the queue before moving on, due to its size. 
   # Docs suggests various solutions. 
    out_q  = Queue() # In order to save the returned dataframes from indep. processes.
    procs  = list()  # Necessary to join the processes later
    output = dict()  # Will be used to store the data put in Queue 

    for i in to_run:
     # Note that we are here using the to run
     # tuple that we are iterating through.
     # we are starting processes for each of these. The tricky part is to
     # solve the args issue. 
      p = Process(target=self.filtered_matrix_producer, kwargs={"dataframe": i[0], "out_q":out_q, "name": i[1]})
      procs.append(p)
      p.start()
      print(f"{self.ERDEFS.OKGREEN}(FILT) In the process loop{self.ERDEFS.END}")

    print(f"{self.ERDEFS.OKGREEN}(FILT) Outside of the process loop{self.ERDEFS.END}")

   # There was an issue here, where the output dict had been referenced before assignment. 
   # However, no error was raised in the process, which was the cause of the headache.  
    for i in range(len(procs)):
      print(f"{self.ERDEFS.OKGREEN}(FILT) Retrieving from Queue{self.ERDEFS.END}")
      data = out_q.get() # Blocking
      output[data[0]] = data[-1]
    print(f"{self.ERDEFS.OKGREEN}(FILT) Finished retrieving{self.ERDEFS.END}")
    print(output)

    for proc in procs:
      proc.join() # Blocking
      print("sdasd")
    print(f"{self.ERDEFS.OKGREEN}(FILT) Joined the processes{self.ERDEFS.END}" )
    print(f"{self.ERDEFS.OKGREEN}(FILT) Incorporating the results{self.ERDEFS.END}")
    print(f"{self.ERDEFS.OKGREEN}(FILT) Saving to csv...{self.ERDEFS.END}")

    output["tmp"].to_csv(f"{self.output_folder}/{tmp_name_unfilt}_filtered_{self.file_name}.csv")
    output["count"].to_csv(f"{self.output_folder}/{count_name_unfilt}_filtered_{self.file_name}.csv")

  def filtered_matrix_producer(self, **kwargs):
    dataframe = kwargs["dataframe"]
    out_q     = kwargs["out_q"]
    name      = kwargs["name"]
    ID        = id(dataframe)
    # The dataframe of interest. 
    # It will contain all the compounded values for each gene. 
    print(f"{self.ERDEFS.OKGREEN}(FILT) Producing filtered dataframes; id: {ID}{self.ERDEFS.END}")

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
    print(f"{self.ERDEFS.OKGREEN}(FILT) Putting shit in the output dict{self.ERDEFS.END}")
    out_q.put(output)
    print(f"{self.ERDEFS.OKGREEN}(FILT) finished putting the shit{self.ERDEFS.END}")


  def read_matrix(self, i):
    # This method is used to read the matrices we have produced.
    # Currently adapted to read from csv's and not tsv's.
    df = DataFrame.from_csv(i)
    return df

  def csv_to_tsv(self, folder="temp_data"):
    # Converts csv's to tsv's due to some 
    # down-stream coding issue. 
    print("[GLOBAL] Running the csv-to-tsv converter")
    fileq = deque(os.listdir(folder))
    
    # This will be checking whether the q is empty or not
    # and will continue if it ain't empty. 
    while fileq:
      
      check_file = fileq.pop()
      path = f"{folder}/{check_file}"
      print(f"(GLOBAL) Converting {path}")
      with open(path, "r") as file:
        data = file.read()

      data = data.replace(",", "	")
      with open(path.replace("csv","tsv"), "w") as file:
        file.write(data)
        
              
  def unifier(self, folder, delimiter="/"):
   # Will remove all of the names from gene-names
   # that were not paired with the 
 
    print("[GLOBAL] Loading queue...")
    fileq = deque([ i for i in os.listdir(folder) if "_filtered_" in i and ".csv" in i ])
    

    while fileq:
        
      filtered   = fileq.pop()
      unfiltered = filtered.replace("filtered", "unfiltered")


      filtered_path   = f"{folder}{delimiter}{filtered}"
      unfiltered_path = f"{folder}{delimiter}{unfiltered}"

      print(f"[GLOBAL] Opening {filtered_path} and {unfiltered_path}")

      filtered_data   = self.arr(filtered_path)
      unfiltered_data = self.arr(unfiltered_path) 
      unfiltered_data = list(set([ i[0].split("|")[-1] for i in unfiltered_data[1:]] ))

      print("[GLOBAL] Created data lists")

      names = filtered_data[0]
      print("[GLOBAL] Loaded gene-names")

      filtered_data = filtered_data[1:]
      unified_data  = []
      print("[GLOBAL] Identifying outliers")
      for pos,i in enumerate(filtered_data):

        if i[0] in unfiltered_data:
          unified_data.append(i)

      new_file = f"{ filtered_path.replace('.csv','') }_unified.tsv"

      print(f"[Globa] writing to {new_file}_unified.tsv")
      with open(f"{new_file}", "w") as file:
        
        names = "\t".join(names)
        names = f"{names}\n"
        file.write(names)

        for i in unified_data:
          string_to_write = "\t".join(i)
          string_to_write = f"{string_to_write}"
          file.write(string_to_write) 

  def arr(self, i, delimiter=",", func=None):
    # Reads input and spits out a 2D array.

    with open(i, "r") as file:
      data = file.readlines()
    data = [i.split(delimiter) for i in data]

    return data
    
    

 
  def text_processing(self, i, parent=None):
    # This generator is central to all of the producers.
    # It is used to process the salmon quant.sf files,
    # performs dictionary lookups and appends the gene names to the ENST's.
    tmp_creator   = lambda data, sample_tmp: DataFrame(data, columns=[sample_tmp], index=data["Name"])
    count_creator = lambda data, sample_count: DataFrame(data, columns=[sample_count], index=data["Name"])
    
    for pos,n in enumerate(i):
      print(f"({parent}) {pos}")

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

#      if pos == 0:
#        return StopIteration

def expander(top, bol_=lambda i: "quant.sf" in i, delimiter="/"):
  walk_gen = os.walk(top)
  print(top)

  paths = [ i for i in walk_gen for n in i if bol_(n) ]
  #print(paths)

  abspaths = [ f"{i[0]}{delimiter}quant.sf" for i in paths]
  print(f"[GLOBAL] Going to process {len(abspaths)} quant files")
  return abspaths

if __name__ == "__main__":
  process = MatrixCreator()
  process.unifier("temp_data")







 
