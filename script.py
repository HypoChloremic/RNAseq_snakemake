def isgz(i):
  bool = re.search(r"\w\.gz", i)

  try:
    assert bool
    return True

  except AssertionError:
    return False

def file_expander(n=3):
  

