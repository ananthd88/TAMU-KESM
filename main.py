import sys
import math
import Timer
import Slice
import IO
import argparse

from collections import deque
import numpy as np
import matplotlib.pyplot as plt

from skimage import data, filter, io
from skimage import img_as_ubyte
from scipy import ndimage

from skimage.morphology import watershed
from skimage.feature import peak_local_max

from mpl_toolkits.mplot3d import Axes3D


# Case 1, 03941_0_74_62.jpg.gif top vessel 
# Explored correctly in full
#sliceID = 1
#seed = (22, 45) 
# testforspecialMerge(): distiI, distjJ <= 10.0 or 20.0
# testforspecialMerge(): distij, distIJ <= 20.0 or 35.0
# searchScope = 35.0
# traceBoundaryVessels = False --> 1 - 118
# traceBoundaryVessels = True  --> 1 - 122


# Case 2, 03941_0_74_62.jpg.gif central vessel
# Explored correctly in full
#sliceID = 1
#seed = (81, 25) 
# testforspecialMerge(): distiI, distjJ <= 20.0
# testforspecialMerge(): distij, distIJ <= 35.0
# searchScope = 35.0
# traceBoundaryVessels = False --> 1 - 84
# traceBoundaryVessels = True  --> 1 - 156 --------> Demo moving from 120 to 119 the edge vessel merge not being detected

# Case 3, 04007_0_74_62.jpg.gif, to test if network can be explored from any 
# seed point.
# Explored correctly in full
#sliceID = 63
#seed = (74, 53) 
# testforspecialMerge(): distiI, distjJ <= 10.0 or 20.0
# testforspecialMerge(): distij, distIJ <= 20.0 or 35.0
# searchScope = 35.0
# traceBoundaryVessels = False --> 1 - 118
# traceBoundaryVessels = True  --> 1 - 122

# Case 4, 04156_0_74_62.jpg.gif, another full trace
# Explored correctly.
#sliceID = 210
#seed = (28, 68) 
# searchScope = 35.0
# testforspecialMerge(): distiI, distjJ <= 10.0
# testforspecialMerge(): distij, distIJ <= 20.0
   # traceBoundaryVessels = False --> 158 - 224
   # traceBoundaryVessels = True  --> Error (Duplicates)
# testforspecialMerge(): distiI, distjJ <= 20.0
# testforspecialMerge(): distij, distIJ <= 35.0
   # traceBoundaryVessels = False --> 158 - 224
   # traceBoundaryVessels = True  --> Error (Duplicates)

# Synthetic case 1,
#sliceID = 1
#seed = (10, 10)
# searchScope = 35.0
# testforspecialMerge(): distiI, distjJ <= 10.0
# testforspecialMerge(): distij, distIJ <= 20.0
   # traceBoundaryVessels = False --> 158 - 224
   # traceBoundaryVessels = True  --> Error (Duplicates)
# testforspecialMerge(): distiI, distjJ <= 20.0
# testforspecialMerge(): distij, distIJ <= 35.0
   # Incomplete exploration

# Case 5, 03997_0_74_62.jpg.gif, Another point.
# Not explored correctly!!!
#sliceID = 53
#seed = (10, 80) 

# Case , Slice 43 to test backtracing
#sliceID = 1
#seed = (93, 49)

# Case , Slice 44 to debug a specialmerge that shouldn't be
#sliceID = 1
#seed = (83.8, 76.9)


#oldCentroid = False # 03961_0_74_62.jpg.gif
vessel = Slice.Vessel()
isInteractive = True

def BFSSetVisitedAs(queue, visited):
   while queue:
      vessel = queue.popleft()
      vessel.setVisitedAs(visited)
      for next in vessel.getNextVessels():
         if next.visited != visited:
            queue.append(next)
      for prev in vessel.getPreviousVessels():
         if prev.visited != visited:
            queue.append(prev)
      
def writeToVTKFile(vtkFilename, queue, sliceID, seed):
   coordinates = []
   cells = []
   
   while queue:
      vessel = queue.popleft()
      vessel.setVisitedAs(True)
      for next in vessel.getNextVessels():
         cells.append((vessel.coordinates[-1][3], next.coordinates[0][3]))
         if not next.visited:
            queue.append(next)
      for previous in vessel.getPreviousVessels():
         cells.append((previous.coordinates[-1][3], vessel.coordinates[0][3]))
         if not previous.visited:
            queue.append(previous)
      for i in range(len(vessel.coordinates)):
         coordinates.append([vessel.coordinates[i][3],    \
                             vessel.coordinates[i][1][0], \
                             vessel.coordinates[i][1][1], \
                             vessel.coordinates[i][0]])   
         if i != 0:
            cells.append((vessel.coordinates[i-1][3], vessel.coordinates[i][3]))
   
   
   vtkFile = open(vtkFilename, 'w')
   vtkFile.write("# vtk DataFile Version 1.0\n" + \
                 "Test vessel\n" + \
                 "ASCII\n\n" + \
                 "DATASET UNSTRUCTURED_GRID\n" + \
                 "POINTS " + str(len(coordinates)) + " float\n")
   coordinates.sort(key=lambda k: k[0])

   for coordinate in coordinates:
      # The X ([1]) & Y ([2]) coordinates are inverted since that is the format
      # convention for the vtk file.
      vtkFile.write(str(int(coordinate[2])) + " " + \
                    str(int(coordinate[1])) + " " + \
                    str(int(coordinate[3])) + "\n")
   vtkFile.write("\nCELLS " + str(len(cells)) + " " + str(len(cells) * 3) + \
                 "\n")
   
   for cell in cells:
      vtkFile.write("2 " + str(cell[0]) + " " + str(cell[1]) + "\n")
   
   vtkFile.write("\nCELL_TYPES " + str(len(cells)) + "\n")
   
   for cell in cells:
      vtkFile.write("4\n")
   
   vtkFile.close()
   

def BFS(queue, sliceID, seed):
   global isInteractive
   
   if sliceID == 1:
      first = True
   else:
      first = False
   colors = ['r', 'g', 'b', 'k', 'm', 'y']
   lencolors = len(colors)
   fig = plt.figure()
   ax = fig.gca(projection='3d')
   minz = 270
   maxz = 0
   
   i = 0
   while queue:
      vessel = queue.popleft()
      if vessel.isExplored():
         continue
      if first:
         vessel.beginVessel('out-of-scope')
         first = False
      if not vessel.endExplored:
         if sliceID and seed:
            vessel.exploreToEnd(sliceID, seed)
            sliceID = False
         else:
            prevCoordinate = vessel.getPreviousCoordinate()
            vessel.exploreToEnd(prevCoordinate[0] + 1, None)         
      if not vessel.beginningExplored:
         if sliceID and seed:
            vessel.exploreToBeginning(sliceID, seed)
            sliceID = False
         else:
            vessel.setExploreInReverse(True)
            prevCoordinate = vessel.getPreviousCoordinate()
            vessel.exploreToBeginning(prevCoordinate[0] - 1, None)
      exits = []
      for next in vessel.getNextVessels():
         if not next.isExplored():
            if len(next.coordinates) == 0:
               print "One of the next vessels, do not have any coordinates."
            exits.append(next)
            #queue.append(next)
      for previous in vessel.getPreviousVessels():
         if not previous.isExplored():
            if len(previous.coordinates) == 0:
               print "One of the previous vessels, do not have any coordinates."
            exits.append(previous)
            #queue.append(previous)
      #print "Vessel Color: " + colors[vessel.vesselID % lencolors]
      vessel.printVessel()
      
      if isInteractive:
         j = 1
         print "Vessels leading along +ve z-axis:"
         for next in vessel.getNextVessels():
            if not next.isExplored():
               coordinate = next.getPreviousCoordinate()
               print "\t" + str(j) + ". " + "Slice - " + str(coordinate[0]) + \
                     " - (" + str(coordinate[1][0]) + ", " + \
                     str(coordinate[1][1]) + ")"
               j += 1
         print "Vessels leading along -ve z-axis:"
         for prev in vessel.getPreviousVessels():
            if not prev.isExplored():
               coordinate = prev.getPreviousCoordinate()
               print "\t" + str(j) + ". " + "Slice - " + str(coordinate[0]) + \
                     " - (" + str(coordinate[1][0]) + ", " + \
                     str(coordinate[1][1]) + ")"
               j += 1
         if j != 1:
            choice = None
            while not choice:
               choice = int(raw_input("Enter the index of the vessel youd like " +\
                                      "to explore -- "))
               if choice < 1 or choice >= j:
                  choice = None
                  continue
            queue.append(exits[choice - 1])
               
               # TODO:::
               
      else:
         for exit in exits:
            queue.append(exit)
      j = 0
      for coordinate in vessel.coordinates:
         if maxz < coordinate[0]:
            maxz = coordinate[0]
         if minz > coordinate[0]:
            minz = coordinate[0]
         ax.scatter(coordinate[1][0], coordinate[1][1], coordinate[0], \
                    zdir='z', c=colors[vessel.vesselID % lencolors],)
         vessel.coordinates[j].append(i)
         j += 1
         i += 1
   print "Range of network: " + str(minz) + " - " + str(maxz)
   ax.legend()
   ax.set_xlim3d(0, 120)
   ax.set_ylim3d(0, 100)
   ax.set_zlim3d(minz, maxz)
   ax.set_xlabel('X Axis')
   ax.set_ylabel('Y Axis')
   ax.set_zlabel('Z Axis')
   plt.show()
   
   return vessel

   
def exploreVessel(vessel):
   value = vessel.exploreToEnd(1, seed) # Case 1, 2
   #value = vessel.exploreToBeginning(43, seed) # Case 3
   if value:
      print "Vessel exploration was successful"
   else:
      print "Vessel exploration was unsuccessful"
   vessel.printVessel()

def main():
   parser = argparse.ArgumentParser(prog = "Tracer", description='Tracer Program')
   parser.add_argument('-e', action='store_true', help = 'Enable tracing of vessels that touch the image boundary')
   parser.add_argument('-i', action='store_true', help = 'Enable interactive mode')
   parser.add_argument('--version', action='version', version='The Tracer 1.0')
   parser.add_argument('--images', metavar = '<image files>', type = str, nargs = '+', required = False,
                      help = 'Image files')                   
   parser.add_argument('-s', metavar='<Synthetic data no.>', type = int, nargs = 1, choices = [1], required = False,
                      help = 'Use one of the synthetic data')
   parser.add_argument('-c', metavar='<coordinate>', type = float, nargs = 3, required = True,#choices = [1, 2],
                      help = 'Coordinates (Slice No, x, y) | -c <sliceID> <x> <y>')
   parser.add_argument('-r', metavar='<Search Radius>', type = float, nargs = 1, required = False,
                      help = 'Search radius')
   parser.add_argument('-v', metavar='<Threshold for validation>', type = float, nargs = 1, required = False,
                      help = 'Threshold for validation')
   parser.add_argument('-w', action='store_true', help = 'Enable warnings')
   
   args = parser.parse_args(sys.argv[1:])
   global isInteractive
   
   if args.r:
      Slice.setSearchScope(args.r[0])
   else:
      Slice.setSearchScope(35.0)
      
   if args.v:
      Slice.setThresholdValidation(args.v[0])
   else:
      Slice.setThresholdValidation(0.36)
      
   if args.e:
      Slice.setTraceBoundaryVessels(True)
   else:
      Slice.setTraceBoundaryVessels(False)
   
   if args.w:
      Slice.setEnableWarnings(True)
   else:
      Slice.setEnableWarnings(False)
   
   if args.i:
      isInteractive = args.i
   else:
      isInteractive = False
   
   if args.c:
      sliceID = int(args.c[0])
      seed = (args.c[1], args.c[2])
   
   if args.s:
      IO.setSynthetic(args.s[0])
      vtkFilename = "synth_" + str(sliceID) + "_" + str(seed[0]) + "_" + \
                    str(seed[1]) + ".vtk"
   else:
      if not args.images:
         print "Please provide list of images for non-synthetic cases."
         parser.print_usage()
         return
      else:
         IO.setSynthetic(False)
         IO.setImageFiles(args.images)
         vtkFilename = str(sliceID) + "_" + str(seed[0]) + "_" + \
                    str(seed[1]) + ".vtk"
      
   
   global vessel
   
   queue = deque()
   queue.append(vessel)
   vessel = BFS(queue, sliceID, seed)
   
   queue = deque()
   queue.append(vessel)
   BFSSetVisitedAs(queue, False)
   
   queue = deque()
   queue.append(vessel)
   writeToVTKFile(vtkFilename, queue, sliceID, seed)
   
main()
