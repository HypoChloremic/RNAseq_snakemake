import re
def isgz(i):
  param = re.search(r"\w\.gz$", i)
  try:
    assert param
    return True
  except AssertionError:
    return False

print(isgz("1231.gz."))