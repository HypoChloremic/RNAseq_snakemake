from os import walk

class Expander:
  def __init__(self):
    pass

  # def expander(self, top, bol_=lambda i: ["quant.sf"] in i, delimiter="\\"):
  #   walk_gen = walk(top)
  #   paths = [ i for i in walk_gen if bol_(i) ]
  #   abspaths = [ f"{i[0]}{delimiter}quant.sf" for i in paths]
  #   print(abspaths)
  #   print()

  def expander(self, top, bol_=lambda i: "quant.sf" in i, delimiter="\\"):
    walk_gen = walk(top)
    paths = [ i for i in walk_gen for n in i if bol_(n) ]
    abspaths = [ f"{i[0]}{delimiter}quant.sf" for i in paths]
    print(abspaths)
    print()


if __name__ == '__main__':
  parent_dir = (r"C:\Users\Ali Rassolie\Documents\GitHub\RNAseq_snakemake\paired_end\output", 
                r"C:\Users\Ali Rassolie\Documents\GitHub\RNAseq_snakemake\single_end\output")
  run = Expander()
  for i in parent_dir:
    run.expander(top=i)