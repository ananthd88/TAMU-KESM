import sys
import Synthesis
from skimage import io


iseries = None

synthetic = False
imageFiles = False

def setSynthetic(val):
   global synthetic
   synthetic = val

def setImageFiles(val):
   global imageFiles
   imageFiles = val
   
def getImage(sliceID):
   # Use synthetic data
   global iseries
   global synthetic
   global imageFiles
   
   if synthetic:
      if not iseries:
         iseries = Synthesis.generate(synthetic)
      
      if sliceID < 1 or sliceID > iseries.getNumImages():
         return False
      return iseries.images[sliceID - 1]
      return False
   else:
      #if sliceID < 1 or sliceID > len(sys.argv) - 2:
      #   print "IO Error: slice ID out of bounds"
      #   return False
      #return io.imread(sys.argv[sliceID + 1])
      if sliceID < 1 or sliceID > len(imageFiles):
         print "IO Error: slice ID out of bounds"
         return False
      return io.imread(imageFiles[sliceID - 1])


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

def showAll():
   iseries.showAll()
