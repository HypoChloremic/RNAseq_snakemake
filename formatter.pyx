from time import time, sleep

def addcolumns(i="data.csv", delimiter=","):
  print(f"[GLOBAL] Opening {i}")

  with open(i) as file:
    data = file.read()

  print(f"[GLOBAL] Splitting data by newline")
  data = data.split("\n")

  # We need to do this to prevent any errors that may rise
  # when the first row of data is parsed. 
  transcript_files = data[0]
  data = data[1:]
  print(data[:10])


  print("[GLOBAL] Sorting data")
  data.sort(key=lambda i: i.split(",")[0].split("|")[-1])
  print(data[:10]) 

  split_data = [i.split(delimiter) for i in data]

  print("[GLOBAL] Getting gene names")
  gene_names = [ i.split(delimiter)[0].split("|")[-1] for i in data ]


  print("[GLOBAL] separating names from values")
  data_float = []
  for i in split_data:
    float_it = []
    for k in i[1:]:
      float_it.append(float(k))
    data_float.append(float_it) 



  prev_name   = None
  prev_values = None
  final_value = []
  final_name  = [] 
  for name,values in zip(gene_names, data_float):
    if name is "":
      continue
    if prev_name == None:
      prev_name   = name
      prev_values = values

    elif prev_name != None:
      if prev_name == name:
        temp = []
        for z,k in zip(prev_values, values):
          temp.append(z+k)
     
        prev_values = list(temp)

      elif prev_name != name:
        final_name.append(prev_name)
        final_value.append(prev_values)
        prev_name = name
        prev_values = values

  return final_value, final_name








      
