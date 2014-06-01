import sys
import math
import Timer
import Slice

from collections import deque
import numpy as np
import matplotlib.pyplot as plt

from skimage import data, filter, io
from skimage import img_as_ubyte
from scipy import ndimage

from skimage.morphology import watershed
from skimage.feature import peak_local_max

from mpl_toolkits.mplot3d import Axes3D


# Load picture and detect edges
# 03972_0_74_62.jpg.gif - Image with a streak
# 03941_0_74_62.jpg.gif - Standard image with no abnormalities except a partial 
#  circle
# 03980_0_74_62.jpg.gif - A very abnormal streak
# 03988_0_74_62.jpg.gif - An ellipse splitting
# 03989_0_74_62.jpg.gif - An ellipse split
# 03960_0_74_62.jpg.gif - A cirple just before splitting
# 03961_0_74_62.jpg.gif - A cirple just after splitting

# Case 1, 03941_0_74_62.jpg.gif top vessel 
# Explored correctly in full
#sliceID = 1
#seed = (22, 45) 

# Case 2, 03941_0_74_62.jpg.gif central vessel
# Explored correctly in full
#sliceID = 1
#seed = (81, 25) 

# Case 3, 04007_0_74_62.jpg.gif, to test if network can be explored from any 
# seed point.
# Explored correctly in full
#sliceID = 63
#seed = (74, 53) 

# Case 4, 04156_0_74_62.jpg.gif, another full trace
# Explored correctly.
sliceID = 210
seed = (28, 68) 

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



def BFS(queue, sliceID, seed):
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
   while queue:
      vessel = queue.popleft()
      #print "Exploring vessel: " + str(vessel.vesselID)
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
            prevCoordinate = vessel.getPreviousCoordinate()
            vessel.exploreToBeginning(prevCoordinate[0] - 1, None)
      for next in vessel.getNextVessels():
         if not next.isExplored():
            if len(next.coordinates) == 0:
               print "One of the next vessels, do not have any coordinates."
            queue.append(next)
      for previous in vessel.getPreviousVessels():
         if not previous.isExplored():
            if len(previous.coordinates) == 0:
               print "One of the previous vessels, do not have any coordinates."
            queue.append(previous)
      vessel.printVessel()
      for coordinate in vessel.coordinates:
         if maxz < coordinate[0]:
            maxz = coordinate[0]
         if minz > coordinate[0]:
            minz = coordinate[0]
         ax.scatter(coordinate[1][0], coordinate[1][1], coordinate[0], \
                    zdir='z', c=colors[vessel.vesselID % lencolors])
   ax.legend()
   ax.set_xlim3d(0, 120)
   ax.set_ylim3d(0, 100)
   ax.set_zlim3d(minz, maxz)
   ax.set_xlabel('X Axis')
   ax.set_ylabel('Y Axis')
   ax.set_zlabel('Z Axis')
   plt.show()
   
def exploreVessel(vessel):
   value = vessel.exploreToEnd(1, seed) # Case 1, 2
   #value = vessel.exploreToBeginning(43, seed) # Case 3
   if value:
      print "Vessel exploration was successful"
   else:
      print "Vessel exploration was unsuccessful"
   vessel.printVessel()

   #for next in vessel.getNextVessels():
   #   print
   #   print
   #   prevCoordinate = next.getPreviousCoordinate()
   #   value = next.exploreToEnd(prevCoordinate[0] + 1, None)
   #   if value:
   #      print "Vessel exploration was successful"
   #   else:
   #      print "Vessel exploration was unsuccessful"
   #   next.printVessel()
   
queue = deque()
queue.append(vessel)
BFS(queue, sliceID, seed)


#exploreVessel(vessel)
