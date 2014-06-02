import sys
from skimage import io

def getImage(sliceID):
   # Use synthetic data
   if sys.argv[1]:
      return False
   else:
      if sliceID < 1 or sliceID > len(sys.argv) - 2:
         print "IO Error: slice ID out of bounds"
         return False
      return io.imread(sys.argv[sliceID + 1])

def hasImage(sliceID):
   if sys.argv[1]:
      if sliceID < 1 or sliceID > sys.argv[2]:
         return False
      else:
         return True
   else:
      if sliceID < 1 or sliceID > len(sys.argv) - 2:
         return False
      else:
         return True
