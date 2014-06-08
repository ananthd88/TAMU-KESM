from skimage import filter, io
from scipy import ndimage
from collections import deque
from skimage.segmentation import clear_border
from skimage.morphology import label
from skimage.measure import regionprops
from skimage.morphology import watershed
from skimage.feature import peak_local_max
import matplotlib.pyplot as plt
import numpy
import math
import sys
import IO

vesselID = 1
searchScope = 35.0
thresholdDetection = 0.0
thresholdValidation = 0.36
enableWarnings = False
traceBoundaryVessels = False
coordinatesSet = set()

def setTraceBoundaryVessels(val):
   global traceBoundaryVessels
   traceBoundaryVessels = val
   
def setEnableWarnings(val):
   global enableWarnings
   enableWarnings = val
   

def setSearchScope(val):
   global searchScope
   searchScope = val   

def setThresholdValidation(val):
   global thresholdValidation
   thresholdValidation = val   

class Slice:
   def __init__(self, sliceID):
      self.image = IO.getImage(sliceID)
      self.sliceID = sliceID
      self.edges = None
      self.filtered = None
      self.regions = None
      self.maxI = self.image.shape[0]
      self.maxJ = self.image.shape[1]
   def getImage(self):
      return self.image
   def getSliceID(self):
      return self.sliceID
   def getFiltered(self):
      return self.filtered
   def filterOut(self, threshold):
      self.filtered = numpy.zeros((self.maxI, self.maxJ))
      for i in range(self.maxI):
         for j in range(self.maxJ):
            if self.image[i][j] < threshold:
               self.filtered[i][j] = False
            else:
               self.filtered[i][j] = True
      return self.filtered
   def findEdges(self, sigmaValue):
      # sigma=2 slightly better than sigma=3 (2 gives smoother edges)
      # sigma=1 helpful in case of splitting circles or ellipses.
      # high_threshold=500 completely eliminates edge detection
      self.edges = filter.canny(self.image, sigma=sigmaValue, low_threshold=10, high_threshold=50)
      return self.edges
         
   # Should be called only after findEdges()
   # Add edges along the boundary for regions that touch the boundary
   def correctEdges(self):
      flag = False
      for i in range(self.maxI):
         if self.image[i][0] > 128:
            flag = True
         else:
            flag = False
         self.edges[i][0] = flag
      for i in range(self.maxI):
         if self.image[i][self.maxJ - 1] > 128:
            flag = True
         else:
            flag = False
         self.edges[i][self.maxJ - 1] = flag
      for j in range(self.maxJ):
         if self.image[0][j] > 128:
            flag = True
         else:
            flag = False
         self.edges[0][j] = flag
      for i in range(self.maxI):
         if self.image[self.maxI - 1][j] > 128:
            flag = True
         else:
            flag = False
         self.edges[self.maxI - 1][j] = flag
      return self.edges
   def isClosedEdge(self, edge):
      blah = 0
      for point in edge:
         if self.isOnBoundingEdge(point):
            return False
         i = point[0]
         j = point[1]
         neighbors = [(i-1, j-1), \
                      (i-1, j)  , \
                      (i-1, j+1), \
                      (i  , j-1), \
                      (i  , j+1), \
                      (i+1, j-1), \
                      (i+1, j)  , \
                      (i+1, j+1)]
         first = neighbors[0]
         flag = self.edges[first[0]][first[1]]
         count = 0
         for i in range(1, 8):
            (j, k) = neighbors[i]
            if (self.edges[j][k] and not flag) or \
               (not self.edges[j][k] and flag):
               flag = not flag
               count += 1
               if count == 3:
                  break
         if count < 3:
            print " found a point with count < 3 at point no. " + str(blah) + " (" + str(point[0]) + ", " + str(point[1]) + ")",
            return False
         blah += 1
      return True
   def getEdges(self):
      return self.edges
   
   def findAllRegions(self):
      cleared = self.filtered.copy()
      #clear_border(cleared)
      label_image = label(cleared)
      #borders = numpy.logical_xor(self.filtered, cleared)
      #label_image[borders] = 0
      return regionprops(label_image, ['Area', 'BoundingBox', 'Centroid'])
   def findEnclosingEdges(self, point):
      topEdge = set()
      bottomEdge = set()
      leftEdge = set()
      rightEdge = set()
      print "Point given to findEnclosingEdges = " + str(point[0]) + ", " + str(point[1])
      i = point[0]
      while i >= 0:
         if self.edges[i][point[1]]:
            self.edgeCollector(topEdge, i, point[1])
            break
         i -= 1
      
      i = point[0] + 1
      while i < self.edges.shape[0]:
         if self.edges[i][point[1]]:
            self.edgeCollector(bottomEdge, i, point[1])
            break
         i += 1
      
      j = point[1]
      while j >= 0:
         if self.edges[point[0]][j]:
            self.edgeCollector(leftEdge, point[0], j)
            break
         j -= 1
      
      j = point[1] + 1
      while j < self.edges.shape[1]:
         if self.edges[point[0]][j]:
            self.edgeCollector(rightEdge, point[0], j)
            break
         j += 1
      
      flag = 0
      if topEdge == rightEdge:
         flag |= 1
      if rightEdge == bottomEdge:
         flag |= 2
      if bottomEdge == leftEdge:
         flag |= 4
      if leftEdge == topEdge:
         flag |= 8
      
      return (topEdge, rightEdge, bottomEdge, leftEdge, flag)
      
   def findBoundingRectangle(self, edgePoints):      
      mini = self.image.shape[0] 
      minj = self.image.shape[1]
      maxi = 0
      maxj = 0
      for point in edgePoints:
         i, j = point
         if i < mini:
            mini = i
         if i > maxi:
            maxi = i
         if j < minj:
            minj = j
         if j > maxj:
            maxj = j
      return ( (mini, minj), (maxi, maxj) )
      
   def edgeCollector(self, points, seedi, seedj):
      #points = set()
      
      candidates = deque()
      candidates.append((seedi, seedj))
      returnValue = False
      while True:
         try:
            point = candidates.popleft()
         except IndexError:
            break
         else:
            if point not in points and self.edges[point[0]][point[1]]:
               points.add(point)
               returnValue = True
               if point[0] != 0:
                  candidates.append((point[0] - 1, point[1]))
                  if point[1] != 0:
                     candidates.append((point[0] - 1, point[1] - 1))
                  if point[1] != self.maxJ:
                     candidates.append((point[0] - 1, point[1] + 1))
               if point[0] != self.maxI:
                  candidates.append((point[0] + 1, point[1]))
                  if point[1] != 0:
                     candidates.append((point[0] + 1, point[1] - 1))
                  if point[1] != self.maxJ:
                     candidates.append((point[0] + 1, point[1] + 1))
               if point[1] != 0:
                  candidates.append((point[0], point[1] - 1))
               if point[1] != self.maxJ:
                  candidates.append((point[0], point[1] + 1))
      return returnValue
   def findCentroid(self, points):
      oldCentroid = (0.0, 0.0)
      i = 0
      if len(points) == 0:
         print "No edge points"
      
      for point in points:
         newCentroid = (oldCentroid[0] + (point[0] - oldCentroid[0])/(i+1), oldCentroid[1] + (point[1] - oldCentroid[1])/(i+1))
         i += 1
         oldCentroid = newCentroid
      return newCentroid
   def findCentroidOfRegion(self, startI, startJ):
      edgePoints = self.edgeCollector(startI, startJ, self.edges)
      centroid = self.findCentroid(edgePoints)
      return centroid
   def distance(self, point1, point2):
      a = point1[0]
      b = point1[1]
      c = point2[0]
      d = point2[1]
      return math.sqrt((a-c)*(a-c) + (b-d)*(b-d))
   def isOnBoundingEdge(self, point):
      #return point[0] == 0 or point[0] == self.maxI - 1 or point[1] == 0 or point[1] == self.maxJ - 1
      return point[0] == 0 or point[0] == 1 or \
             point[0] == self.maxI or point[0] == self.maxI - 1 or \
             point[1] == 0 or point[1] == 1 or \
             point[1] == self.maxJ or point[1] == self.maxJ - 1
      
class vesselCrossSection:
   def __init__(self, periphery, centroid, sliceID):
      self.periphery = periphery
      self.centroid = centroid
      self.sliceID = sliceID
class regionOfInterest(Slice):
   def __init__(self, slice1, startI, startJ, lengthI, lengthJ):
      self.startI = startI
      self.startJ = startJ
      self.lengthI = lengthI
      self.lengthJ = lengthJ
      self.filtered = slice1.filtered
      self.edges = slice1.getEdges()
      self.maxI = slice1.maxI
      self.maxJ = slice1.maxJ
      self.filtered = self.filtered[startI:startI+lengthI, startJ:startJ+lengthJ]
      self.edges = self.edges[startI:startI+lengthI, startJ:startJ+lengthJ]
      self.sliceID = slice1.sliceID
      self.image = None
      self.filled = None
      
      # TODO: Complete the function.
   def getROI(self):
      return self.image
   def fillRegion(self):
      self.filled = ndimage.binary_fill_holes(self.edges)
      return self.filled
   
   
class Vessel:
   def __init__(self):
      global vesselID
      self.vesselID = vesselID
      vesselID += 1
      self.coordinates = deque()
      self.next = set()
      self.previous = set()
      self.beginningExplored = False
      self.endExplored = False
      self.exploreInReverse = False # If the vessel has to be explored in reverse
      self.endReason = None
      self.beginningReason = None
      self.visited = False
      
      # Auxiliary variables, to be cleared after the vessel has been fully 
      # explored.
      self.onBoundary = False
      self.slice = None
      self.prevRegions = None
      self.prevRegionIndex = -1
      self.regionIndex = -1
      self.matches = None
      self.matchedPrevRegions = None
      self.unmatchedPrevRegions = None
      self.matchedRegions = None
      self.unmatchedRegions = None
      self.closestRegions = None
      self.closestPrevRegions = None
      self.currentSliceID = None
      
      # TODO: Maybe delay these two computations till the time the vessel is 
      # being processed
      self.regions = None
      self.distances = None
      
   def setVisitedAs(self, visited):
      self.visited = visited
   def setExploreInReverse(self, val):
      self.exploreInReverse = val
   def isExplored(self):
      return self.beginningExplored and self.endExplored
   def addCoordinate(self, sliceID, region):
      global coordinatesSet
      if self.isExplored():
         print "Vessel already fully explored. Cannot add coordinates."
         return
      if self.beginningExplored and self.exploreInReverse:
         print "Vessel beginning already explored. Cannot add coordinates."
         return
      if self.endExplored and not self.exploreInReverse:
         print "Vessel end already explored. Cannot add coordinates."
         return
      coordinate = (sliceID, int(region['centroid'][0]), \
                    int(region['centroid'][1]))
      if coordinate in coordinatesSet:
         print "Duplicate coordinate:"
         print "Vessel ID: " + str(self.vesselID)
         print "Slice ID: " + str(self.currentSliceID)
         print "Coordinate: " + str(region['centroid'][0]) + ", " + str(region['centroid'][1])
         print "Error. Exiting."
         # TODO: Duplicate coordinates are not ignored, just a warning
         coordinatesSet.add(coordinate)
         #sys.exit(0)
      else:
         coordinatesSet.add(coordinate)
      if self.exploreInReverse:
         self.coordinates.appendleft([sliceID, region['centroid'], \
                                      region['area']])
      else:
         self.coordinates.append([sliceID, region['centroid'], \
                                  region['area']])      
   def addNext(self, other):
      self.next.add(other)
   def addPrevious(self, other):
      self.previous.add(other)
   def getCoordinates(self):
      return self.coordinates
   def getNextVessels(self):
      return self.next
   def getPreviousVessels(self):
      return self.previous
   def getPreviousCoordinate(self):
      if self.exploreInReverse:
         return self.coordinates[0]
      else:
         return self.coordinates[-1]
   def endVessel(self, reason):
      self.endExplored = True
      self.endReason = reason
      
   def beginVessel(self, reason):
      self.beginningExplored = True
      self.beginningReason = reason
      
   def clearVariables(self):
      self.slice = None
      self.prevRegions = None
      self.regions = None
      self.distances = None
      self.matches = None
      self.prevRegionIndex = -1
      self.regionIndex = -1
      self.matchedPrevRegions = None
      self.unmatchedPrevRegions = None
      self.matchedRegions = None
      self.unmatchedRegions = None
      self.closestRegions = None
      self.closestPrevRegions = None
      self.onBoundary = False
      self.currentSliceID = None

   def splitVessel(self, sliceID, indexes):
      next1 = Vessel()
      next1.slice = self.slice
      next1.addCoordinate(sliceID, self.regions[indexes["next"][0]])
      next1.prevRegions = self.regions
      next1.prevRegionIndex = indexes["next"][0]
      
      next2 = Vessel()
      next2.slice = self.slice
      next2.addCoordinate(sliceID, self.regions[indexes["next"][1]])
      next2.prevRegions = self.regions
      next2.prevRegionIndex = indexes["next"][1]
      
      if self.exploreInReverse:
         next1.endVessel('merge')
         next2.endVessel('merge')
         next1.addNext(self)
         next2.addNext(self)
         self.beginVessel('merge')
         self.addPrevious(next1)
         self.addPrevious(next2)
      else:
         next1.beginVessel('split')
         next2.beginVessel('split')
         next1.addPrevious(self)
         next2.addPrevious(self)
         self.endVessel('split')
         self.addNext(next1)
         self.addNext(next2)
      next1.exploreInReverse = self.exploreInReverse
      next2.exploreInReverse = self.exploreInReverse
      self.clearVariables()
   
   def specialMergeVessels(self, sliceID, indexes):
      next1 = Vessel()
      next2 = Vessel()
      other = Vessel()
      if self.exploreInReverse:
         self.addPrevious(next1)
         other.addPrevious(next1)
         other.addPrevious(next2)
         next1.addNext(self)
         next1.addNext(other)
         next2.addNext(other)
         self.beginVessel('specialmerge')
         other.beginVessel('specialmerge')
         next1.endVessel('specialmerge')
         next2.endVessel('specialmerge')

      else:      
         self.addNext(next1)         
         other.addNext(next1)
         other.addNext(next2)         
         next1.addPrevious(self)
         next1.addPrevious(other)         
         next2.addPrevious(other)
         self.endVessel('specialmerge')
         other.endVessel('specialmerge')
         next1.beginVessel('specialmerge')
         next2.beginVessel('specialmerge')

      other.exploreInReverse = not self.exploreInReverse
      next1.exploreInReverse = self.exploreInReverse
      next2.exploreInReverse = self.exploreInReverse
      
      next1.slice = self.slice
      next1.addCoordinate(sliceID, self.regions[indexes["next"][0]])
      next1.prevRegions = self.regions
      next1.prevRegionIndex = indexes["next"][0]
      
      next2.slice = self.slice
      next2.addCoordinate(sliceID, self.regions[indexes["next"][1]])
      next2.prevRegions = self.regions
      next2.prevRegionIndex = indexes["next"][1]
      
      other.addSlice(sliceID - 1)
      other.addCoordinate(sliceID, self.prevRegions[indexes["previous"][1]])
      other.prevRegions = self.prevRegions
      other.prevRegionIndex = indexes["previous"][1]
      
      self.clearVariables()
      
   def mergeVessels(self, sliceID, indexes):
      next = Vessel()
      other = Vessel()      
      other.exploreInReverse = not self.exploreInReverse

      next.exploreInReverse = self.exploreInReverse
      next.slice = self.slice
      next.addCoordinate(sliceID, self.regions[indexes["next"]])
      next.prevRegions = self.regions
      next.prevRegionIndex = indexes["next"]

      
      other.prevRegions = self.prevRegions
      other.prevRegionIndex = indexes["previous"][1]

      if self.exploreInReverse:
         next.addNext(self)
         next.addNext(other)
         self.addPrevious(next)
         other.addPrevious(next)
         self.beginVessel('split')         
         next.endVessel('split')
         other.addSlice(sliceID + 1)
         other.addCoordinate(sliceID + 1, \
                             self.prevRegions[indexes["previous"][1]])
         other.beginVessel('split')
      else:
         next.addPrevious(self)
         next.addPrevious(other)
         self.addNext(next)
         other.addNext(next)
         self.endVessel('merge')         
         next.beginVessel('merge')
         other.addSlice(sliceID - 1)
         other.addCoordinate(sliceID - 1, \
                             self.prevRegions[indexes["previous"][1]])
         other.endVessel('merge')
      
      self.clearVariables()
   
   def addSlice(self, sliceID):
      self.slice = Slice(sliceID)
      self.slice.filterOut(192) # 192
      self.slice.findEdges(2)
      self.slice.correctEdges()
      self.regions = self.slice.findAllRegions()
      
   def distance(self, point1, point2):
      x1 = point1[0]
      y1 = point1[1]
      x2 = point2[0]
      y2 = point2[1]
      return math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))
      
   def areSamePoints(self, point1, point2):
      return self.distance(point1, point2) <= 2.0
      
   def findAllDistances(self):
      self.distances = numpy.zeros((len(self.prevRegions), len(self.regions)))
      for i in range(len(self.prevRegions)):
         for j in range(len(self.regions)):
            self.distances[i][j] = self.distance(self.prevRegions[i]['centroid'], self.regions[j]['centroid'])         
   def findClosestRegion(self, prevRegionIndex):
      mindist = 10000.0
      minj = -1
      for j in range(self.distances.shape[1]):
         if mindist > self.distances[prevRegionIndex][j]:
            mindist = self.distances[prevRegionIndex][j]
            minj = j # mini corresponds to the closest region in the previous 
                     # slice --> Let this be the candidate
      mindist = 10000.0
      mini = -1
      for i in range(self.distances.shape[0]):
         if mindist > self.distances[i][minj]:
            mindist = self.distances[i][minj]
            mini = i # minj corresponds to the region closest to the 'candidate'
                     # in the current slice
      if mini == prevRegionIndex:
         return minj
      return -1
   
   # Should be called only if self.prevRegions exists
   def findMatches(self):
      self.findAllDistances()
      self.matches = []
      for i in range(len(self.prevRegions)):
         self.matches.append(self.findClosestRegion(i))
      #print str(len(self.matches)) + " regions in previous slice."
         
      self.matchedPrevRegions = set()
      self.unmatchedPrevRegions = set()
      self.matchedRegions = set()
      self.unmatchedRegions = set() # TODO: Easy way to do this, set operations
      
      i = 0      
      for prevRegion in self.prevRegions:
         #print str(i) + ": " + str(self.prevRegions[i]['Centroid'][0]) + ", " + str(self.prevRegions[i]['Centroid'][1]),
         j = self.matches[i]
         if j >= 0:
            self.matchedRegions.add(j)
            self.matchedPrevRegions.add(i)
            #print " : " + str(self.regions[j]['Centroid'][0]) + ", " + str(self.regions[j]['Centroid'][1])
         else:
            self.unmatchedPrevRegions.add(i)            
         i += 1
      self.regionIndex = self.matches[self.prevRegionIndex]
      #if self.regionIndex >= 0:
      #   print "Previous cross section of vessel was matched to a " +\
      #        "cross section in this slice."      

      for i in range(len(self.regions)):
         if i not in self.matchedRegions:
            self.unmatchedRegions.add(i)
      
   def orderRegions(self, point):
      distances = {}
      i = 0
      for region in self.regions:
         distances[i] = self.distance(region['Centroid'], point)
         i += 1
      self.closestRegions = [k for (k,v) in sorted(distances.items(), key=lambda (k, v): v)]
   def orderPrevRegions(self, point):
      distances = {}
      i = 0
      for region in self.prevRegions:
         distances[i] = self.distance(region['Centroid'], point)
         i += 1
      self.closestPrevRegions = [k for (k,v) in sorted(distances.items(), key=lambda (k, v): v)]
   
   def warnCloseLimit(self, val, lessthan, limit, tolerance, variableName, string):
      global enableWarnings
      if not enableWarnings:
         return
      if lessthan:
         # if val < limit, it is ignored
         # limit - tolerance < val < limit
         if limit - tolerance <= val and val < limit: 
            print "Warning: " + string + variableName + " = " + str(val) + \
                  ": close to limit of " + str(limit) + ": SliceID = " + \
                  str(self.currentSliceID) + ": Vessel ID = " + \
                  str(self.vesselID)
      else:
         # if val > limit, it is ignored
         # limit < val < limit + tolerance
         if limit < val and val <= limit + tolerance:
            print "Warning: " + string + variableName + " = " + str(val) + \
                  ": close to limit of " + str(limit) + "\nSliceID = " + \
                  str(self.currentSliceID) + "\nVessel ID = " + \
                  str(self.vesselID)
   def testForMerge1(self):
      #  i      j
      #   .    |
      #    .  |
      #     .|
      #     J
      global searchScope
      global thresholdValidation
      global coordinatesSet
      watchSlice = False
      if watchSlice and self.currentSliceID == watchSlice:
         print "Inside testForMerge1(): vessel: " + str(self.vesselID) + \
               " slice: " + str(self.currentSliceID) + " regionIndex = " + \
               str(self.regionIndex)
      i = self.prevRegionIndex
      if self.regionIndex != -1:
         if watchSlice and self.currentSliceID == watchSlice:
            print "Inside testForMerge1(): Point1"      
         return False
      # No need to check for initial area change since we could not project
      # i to I
      for k in range(1, len(self.closestPrevRegions)):
         j = self.closestPrevRegions[k]
         J = self.matches[j]
         if J == -1 or (self.currentSliceID, \
                        int(self.regions[J]['centroid'][0]), \
                        int(self.regions[J]['centroid'][1])) in coordinatesSet:
            continue
         distij = self.distance(self.prevRegions[i]['centroid'], \
                                self.prevRegions[j]['centroid'])
         if watchSlice and self.currentSliceID == watchSlice:
            print "Inside testForMerge1(): Point2: distij = " + str(distij)         
         self.warnCloseLimit(distij, False, searchScope, 5.00, "distij", "Inside testForMerge1(): ")
         if distij >= searchScope:
            return False
         
         oldArea1 = self.prevRegions[i]['area']
         oldArea2 = self.prevRegions[j]['area']
         newArea  = self.regions[J]['area']
         areaChange = math.fabs(float(newArea - oldArea1 - oldArea2)) / \
                      float(oldArea1 + oldArea2)
         if watchSlice and self.currentSliceID == watchSlice:
            print "Inside testForMerge1(): Point3: areaChange = " + \
                  str(areaChange)
            self.printRegion(i, True)
            self.printRegion(j, True)
            self.printRegion(J, False)
         # Case 2: I: Slice 117 (60.7708333333, 28.75) 144.0
         #         i: Slice 118 (55.6470588235, 36.0784313725) 51.0
         #         j: Slice 118 (63.4, 20.9272727273) - 55.0
         #         areaChange = 0.358490566038
         self.warnCloseLimit(areaChange, False, thresholdValidation, 0.05, "areaChange", "Inside testForMerge1(): ")
         if areaChange <= thresholdValidation:
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForMerge1(): Point4: Success"
         
            results = { \
                       "previous": (i, j), \
                       "next"    : (J) \
                      }
            return results
      if watchSlice and self.currentSliceID == watchSlice:
         print "Inside testForMerge1(): Point5 : No merge1"
      return False
      
   def testForMerge2(self):
      #  i      j
      #   \    .
      #    \  .
      #     \.
      #     I
      global searchScope
      global thresholdValidation
      global coordinatesSet
      watchSlice = False
      i = self.prevRegionIndex
      I = self.regionIndex
      if watchSlice and self.currentSliceID == watchSlice:
         print "Inside testForMerge2(): vessel: " + str(self.vesselID) + \
               " slice: " + str(self.currentSliceID)
      
      oldArea1 = self.prevRegions[i]['area']
      newArea = self.regions[I]['area']
      areaChange1 = float(newArea - oldArea1)/float(oldArea1)
      if watchSlice and self.currentSliceID == watchSlice:
         print "Inside testForMerge2(): areaChange1 = " + str(areaChange1)
      
      self.warnCloseLimit(areaChange1, True, 0.30, 0.05, "areaChange1", "Inside testForMerge2(): ")
      if areaChange1 >= 0.30:
         for k in range(1, len(self.closestPrevRegions)):
            j = self.closestPrevRegions[k]
            if self.exploreInReverse:
               addendum = 1
            else:
               addendum = -1
            if (self.currentSliceID + addendum, \
                self.prevRegions[j]['centroid'][0], \
                self.prevRegions[j]['centroid'][1]) in coordinatesSet:
               continue
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForMerge2(): Point1"
               self.printRegion(i, True)
               self.printRegion(j, True)
               self.printRegion(I, False)
         
            distij = self.distance(self.prevRegions[i]["centroid"], \
                                   self.prevRegions[j]["centroid"])
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForMerge2(): Point2: distij = " + str(distij)
            
            self.warnCloseLimit(distij, False, searchScope, 5.00, "distij", "Inside testForMerge2(): ")
            if distij >= searchScope:
               if watchSlice and self.currentSliceID == watchSlice:
                  print "Inside testForMerge2(): i & j are too far: " + \
                        str(distij)            
               return False         
            if self.matches[j] < 0:
               oldArea2 = self.prevRegions[j]['area']
               areaChange2 = math.fabs(float(newArea - oldArea1 - oldArea2)) / \
                             float(oldArea1 + oldArea2)
               # TODO: Is this limit too broad?
               if watchSlice and self.currentSliceID == watchSlice:
                  print "Inside testForMerge2(): Point4: areaChange2 = " + str(areaChange2)
               
               self.warnCloseLimit(areaChange2, False, thresholdValidation, 0.05, "areaChange2", "Inside testForMerge2(): ")
               if areaChange2 <= thresholdValidation:
                  if watchSlice and self.currentSliceID == watchSlice:
                     print "Inside testForMerge2(): Point5: Found Merge2"
            
                  result = {"previous": (i, j), \
                            "next"    : (I)  \
                           }
                  return result
      if watchSlice and self.currentSliceID == watchSlice:
         print "Inside testForMerge2(): Point6: No merge2"
         
      return False

   # Special merge: 
   # Previous slice - Two adjacent regions a & b, b about to split
   # Current slice - Two adjacent regions A & B, part of b split off and merged 
   # with A
   def testForSpecialMerge(self):
      watchSlice = False
      oldArea = self.getPreviousCoordinate()[2]
      newArea = self.regions[self.regionIndex]['area']
      areaChange = float(newArea - oldArea)/float(oldArea)
      if watchSlice and self.currentSliceID == watchSlice:
         print "Inside testForSpecialMerge(): Vessel = " + str(self.vesselID) +\
               " sliceID = " + str(self.currentSliceID) + " areaChange = " + \
               str(areaChange)
      self.warnCloseLimit(areaChange, True, 0.30, 0.05, "areaChange", "Inside testForSpecialMerge(): ")
      if areaChange >= 0.30:
         #print "Inside testforspecialmerge(): areaChnage > 0.30"
         if len(self.closestPrevRegions) < 2 or len(self.closestRegions) < 2:
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForSpecialMerge(): len(closestRegions) or " +\
                     "len(closestPrevregions) == 0"      
            return False
         i = self.closestPrevRegions[0]
         j = self.closestPrevRegions[1]
         I = self.closestRegions[0]
         J = self.closestRegions[1]
         distiI = self.distance(self.prevRegions[i]["centroid"], \
                                self.regions[I]["centroid"])
         self.warnCloseLimit(distiI, False, 20.0, 5.0, "distiI", "Inside testForSpecialMerge(): ")
         if distiI >= 20.0: # 10.0 (case 1, 3, 4), 20.0 (case 1, 2, 3, 4, syth01 (once random))
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForSpecialMerge(): i & I are too far: " + \
                     str(distiI)
            return False
         
         distjJ = self.distance(self.prevRegions[j]["centroid"], \
                                self.regions[J]["centroid"])
                                
         self.warnCloseLimit(distjJ, False, 20.0, 5.0, "distjJ", "Inside testForSpecialMerge(): ")
         if distjJ >= 20.0: # 10.0, 20.0, greater the limit, less constrained condition
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForSpecialMerge(): j & J are too far: " + \
                     str(distjJ)
            return False
         distij = self.distance(self.prevRegions[i]["centroid"], \
                                self.regions[j]["centroid"])
         self.warnCloseLimit(distiI, False, 35.0, 15.0, "distij", "Inside testForSpecialMerge(): ")
         if distij >= 35.0: # 20.0, 35.0, greater the limit, less constrained condition
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForSpecialMerge(): i & j are too far: " + \
                     str(distij)
            return False
         distIJ = self.distance(self.regions[I]["centroid"], \
                                self.regions[j]["centroid"])
         self.warnCloseLimit(distIJ, False, 35.0, 15.0, "distIJ", "Inside testForSpecialMerge(): ")
         if distIJ >= 35.0: # 20.0, 35.0 greater the limit, less constrained condition
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForSpecialMerge(): I & J are too far: " + \
                     str(distIJ)
            return False
         diffijIJ = math.fabs(self.distance(self.prevRegions[i]["centroid"], \
                                            self.prevRegions[j]["centroid"]) - \
                              self.distance(self.regions[I]["centroid"], \
                                            self.regions[j]["centroid"]))
         self.warnCloseLimit(diffijIJ, False, 10.0, 10.0, "diffijIJ", "Inside testForSpecialMerge(): ")
         if diffijIJ >= 10.0: #  greater the limit, less constrained condition
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForSpecialMerge(): (i & j) & (I & J) are too different: " + \
                     str(diffijIJ)
            return False
         
         areaiI = self.prevRegions[i]['area'] - self.regions[I]['area']
         areajJ = self.prevRegions[j]['area'] - self.regions[J]['area']
         if (areaiI > 0 and areajJ > 0) or (areaiI < 0 and areajJ < 0):
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForSpecialMerge(): Change in area is not a " +\
                     "zero sum: areaiI = " + str(areaiI) + ": areajJ = " + \
                     str(areajJ)
            return False
         self.warnCloseLimit(areaiI + areajJ, False, 20.0, 10.0, "areaiI + areajJ", "Inside testForSpecialMerge(): ")
         if areaiI + areajJ >= 20.0:
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForSpecialMerge(): Sum of gain and loss in " +\
                     "area is greater than threshold of () :" + \
                     str(areaiI + areajJ)
            return False
            
         areaij = self.prevRegions[i]['area'] + self.prevRegions[j]['area']
         areaIJ = self.regions[I]['area'] + self.regions[J]['area']
         diffTotalArea = math.fabs(float(areaIJ - areaij))/float(areaij)
         self.warnCloseLimit(diffTotalArea, False, 0.15, 0.05, "diffTotalArea", "Inside testForSpecialMerge(): ")
         #if diffTotalArea >= 15.0:
         if diffTotalArea >= 0.15:
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForSpecialMerge(): Percentage change in " +\
                     "total area is greater than threshold of 15%: " + \
                     str(diffTotalArea)
            return False
         if watchSlice and self.currentSliceID == watchSlice:
            print "Inside testForSpecialMerge(): Successfully detected " +\
                  "specialmerge."
         result = {"previous": (i, j), \
                   "next"    : (I, J)  \
                  }
         return result
      if watchSlice and self.currentSliceID == watchSlice:
         print "Inside testForSpecialMerge(): specialmerge not detected "
   
      return False

   def testForSplit(self):
      #     i
      #    | .
      #   |   .
      #  |     .
      # I      J
      
      global searchScope
      global thresholdValidation
      global coordinatesSet
      watchSlice = False
      i = self.prevRegionIndex
      I = self.regionIndex
      
      
      if watchSlice and self.currentSliceID == watchSlice:
         print "Inside testForSplit(): vessel: " + str(self.vesselID) + \
               " slice: " + str(self.currentSliceID)
         
      # Extract the area from the last coordinate
      oldArea  = self.prevRegions[i]['area']
      newArea1 = self.regions[I]['area']
      areaChange1 = float(newArea1 - oldArea)/float(oldArea)
      if watchSlice and self.currentSliceID == watchSlice:
         print "Inside testForSplit(): Point1 : areaChange1 = " + \
               str(areaChange1)
      # If the area decreased by more than 30% --> Potential split
      # Check for J iff this test is passed
      # Case 2: i: (32.4, 44.5916666667) - 240.0
      #         I: (41.6648648649, 44.3081081081) - 185.0
      #         J: (21.8484848485, 48.0202020202) - 99.0
      #         areaChange = -0.229166666667 ==> <= 0.22
      
      self.warnCloseLimit(areaChange1, False, -0.22, 0.05, "areaChange1", "Inside testForSplit(): ")
      if areaChange1 <= -0.22:
         for k in range(1, len(self.closestRegions)):
            J = self.closestRegions[k]
            if (self.currentSliceID, int(self.regions[J]['centroid'][0]), \
                int(self.regions[J]['centroid'][1])) in coordinatesSet:
               continue
               
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForMerge2(): Point1.5"
               self.printRegion(i, True)
               self.printRegion(I, False)
               self.printRegion(J, False)
         
            distIJ = self.distance(self.regions[I]["centroid"], \
                                   self.regions[J]["centroid"])
            if watchSlice and self.currentSliceID == watchSlice:
               print "Inside testForSplit(): Point2 : distIJ = " + str(distIJ)
            self.warnCloseLimit(distIJ, False, searchScope, 5.0, "distIJ", "Inside testForSplit(): ")
            if distIJ >= searchScope: # 20.0, 35.0
               if watchSlice and self.currentSliceID == watchSlice:
                  print "Inside testForSplit(): Point3 : distIJ >= searchScope "
               return False
            if J in self.unmatchedRegions:
               newArea2 = self.regions[J]['area']
               areaChange2 = math.fabs(float(newArea1 + newArea2 - oldArea))/ \
                             float(oldArea)
               if watchSlice and self.currentSliceID == watchSlice:
                  print "Inside testForSplit(): Point4 : areaChange2 = " + \
                        str(areaChange2)
               self.warnCloseLimit(areaChange2, False, thresholdValidation, 0.05, "areaChange2", "Inside testForSplit(): ")
               if areaChange2 <= thresholdValidation: # TODO: Is this limit too broad?
                  results = { \
                             "previous": (i), \
                             "next"    : (I, J)\
                            }
                  if watchSlice and self.currentSliceID == watchSlice:
                     print "Inside testForSplit(): Point5 : Success"
                  return results
      if watchSlice and self.currentSliceID == watchSlice:
         print "Inside testForSplit(): Point6 : No split."               
      return False
      
   def printRegion(self, index, isPrev):
      if isPrev:
         print "Prev Region: " + str(index) + " - (" + str(self.prevRegions[index]['Centroid'][0]) + ", " + str(self.prevRegions[index]['Centroid'][1]) + ") - " + str(self.prevRegions[index]['Area'])
      else:
         print "Region: " + str(index) + " - (" + str(self.regions[index]['Centroid'][0]) + ", " + str(self.regions[index]['Centroid'][1]) + ") - " + str(self.regions[index]['Area'])

   def printTrace(self):
      if self.prevRegions:
         for i in range(len(self.prevRegions)):
            print str(i) + " - (" + str(self.prevRegions[i]['Centroid'][0]) + ", " + str(self.prevRegions[i]['Centroid'][1]) + ") - " + str(self.prevRegions[i]['Area']) + " --> ",
            j = self.matches[i]
            if j != -1:
               print str(j) + " - (" + str(self.regions[j]['Centroid'][0]) + ", " + str(self.regions[j]['Centroid'][1]) + ") - " + str(self.regions[j]['Area'])
            else:
               print
      
   def actualRegionProximity(self, region1, region2):
      (c1x, c1y) = region1['Centroid']
      (c2x, c2y) = region2['Centroid']
      
      if c1x == c2x:
         if c1y == c2y:
            return 0.0
         else:
            if c1y < c2y:
               less = c1y
               more = c2y
            else:
               less = c2y
               more = c1y
            for j in range(less, more + 1):
               if self.filtered[c1x][j] == False:
                  (x1, y1) = (c1x, j)
                  break
            for j in range(more, less - 1, -1):
               if self.filtered[c1x][j] == False:
                  (x2, y2) = (c1x, j)
                  break
      elif c1y == c2y:
         if c1x < c2x:
            less = c1x
            more = c2x
         else:
            less = c2x
            more = c1x
         for i in range(less, more + 1):
            if self.filtered[i][c1y] == False:
               (x1, y1) = (i, c1y)
               break
         for i in range(more, less - 1, -1):
            if self.filtered[i][c1y] == False:
               (x2, y2) = (i, c1y)
               break
      else:
         if c1x < c2x:
            less = c1x
            more = c2x
         else:
            less = c2x
            more = c1x
         for i in range(less, more + 1):
            j = float(c1y - c2y)*float(i - c1x)/float(c2x - c1x) + c1y
            if self.filtered[i][j] == False:
               (x1, y1) = (i, j)
               break
         for i in range(more, less - 1, -1):
            j = float(c1y - c2y)*float(i - c1x)/float(c2x - c1x) + c1y
            if self.filtered[i][j] == False:
               (x2, y2) = (i, j)
               break
      return self.distance((x1, y1), (x2, y2))
      
   # Pre-condition: sliceID & odlCentroid should be valid and appropriate.
   def extendVessel(self, sliceID, oldCentroid):
      global traceBoundaryVessels
      self.currentSliceID = sliceID
      if not IO.hasImage(sliceID):
         if self.exploreInReverse:
            self.beginVessel('out-of-scope')
         else:
            self.endVessel('out-of-scope')
         self.clearVariables()
         return False

      self.addSlice(sliceID)
      image1 = self.slice.getImage()
      image2 = self.slice.getFiltered()
      image3 = self.slice.getEdges()
      
      candidatePoint = oldCentroid
      
      if len(self.regions) == 0:
         # TODO: if there are no regions in the image
         print "Warning: No regions could be identified in the slice.\n" +\
               "Skipping slice no." + str(sliceID)
         return False
      
      self.orderRegions(candidatePoint)
      
      if self.prevRegions:
         # Match the regions in the previous slice to regions to regions in this
         # slice.
         self.findMatches()
         self.orderPrevRegions(candidatePoint)
         # If the cross sectional region of the vessel in the previous slice
         # could not be matched to a cross sectional region in this slice,
         # NOTE: This particulr case of merging has to be tested here itself.
         testResult = self.testForMerge1()
         if testResult:
            self.mergeVessels(sliceID, testResult)
            return False
         elif self.regionIndex < 0:
            if self.exploreInReverse:
               self.beginVessel('vanished')
            else:
               self.endVessel('vanished')
            self.clearVariables()
            return False
      else:
         self.regionIndex = self.closestRegions[0]
      
      newCentroid = self.regions[self.regionIndex]['centroid']
      newArea = self.regions[self.regionIndex]['area']
      (mini, minj, maxi, maxj) = self.regions[self.regionIndex]['bbox']
      
      # If the vessel moves beyond the edges of the slice, end it.
      if not traceBoundaryVessels and \
         (self.slice.isOnBoundingEdge((mini, minj)) or \
          self.slice.isOnBoundingEdge((maxi, maxj))):
         # TODO: Should this coordinate be added to the queue?
         if self.exploreInReverse:
            self.beginVessel('out-of-scope')
         else:
            self.endVessel('out-of-scope')
         self.clearVariables()
         return False


      returnValue = True
      childConfirmed = False
      childCentroid = None
      
      # Area analysis to determine if vessel split or merged.
      #if len(self.coordinates) > 0:
      if self.prevRegions:
         testResult = self.testForSplit()
         if testResult:
            childCentroid = self.regions[testResult["next"][1]]['Centroid']
            self.splitVessel(sliceID, testResult)
            childConfirmed = True
            return False # Return false if the vessel ended.               
      
         testResult = self.testForMerge2()
         if testResult:
            self.mergeVessels(sliceID, testResult)
            return False
         testResult = self.testForSpecialMerge()
         if testResult:
            self.specialMergeVessels(sliceID, testResult)
            return False
         
         ## Another case here
         #elif self.matches[indexesPrev[2]] < 0:
         #   #print "The 2nd region closest to this region in the previous slice" \
         #   #      " could not be matched to a region in this slice. Merge " \
         #   #      "detected."
         #   returnValue = False
         
         
      if returnValue:
         self.addCoordinate(sliceID, self.regions[self.regionIndex])
         self.prevRegionIndex = self.regionIndex
         self.prevRegions = self.regions
      
      if sliceID < 0:
         plt.figure(figsize=(12, 4))
         
         image4 = image2[mini - 5:maxi + 5, minj - 5:maxj + 5]
         
         distance = ndimage.distance_transform_edt(image4)
         local_maxi = peak_local_max(distance, indices=False, footprint=numpy.ones((3, 3)),
                                     labels=image4)
         markers = ndimage.label(local_maxi)[0]
         labels = watershed(-distance, markers, mask=image4)
         
         
         image3[newCentroid[0]][newCentroid[1]] = True
         if childConfirmed:
            image3[childCentroid[0]][childCentroid[1]] = True
         plt.subplot(141)
         plt.imshow(image1, cmap=plt.cm.gray)
         plt.axis('off')
         plt.title("Slice " + str(sliceID), fontsize=20)

         plt.subplot(142)
         plt.imshow(image2, cmap=plt.cm.gray)
         plt.axis('off')
         plt.title('Filtered', fontsize=20)

         plt.subplot(143)
         plt.imshow(image3, cmap=plt.cm.gray)
         plt.axis('off')
         plt.title('Edges', fontsize=20)

         plt.subplot(144)
         plt.imshow(labels, cmap=plt.cm.gray)
         plt.axis('off')
         plt.title('Redundant', fontsize=20)
         
         plt.show()
      
      return True
   # Pre-condition:  sliceID should always be valid and appropriate.
   #                 seed should be valid if this is a first exploration.
   #           
   def exploreToEnd(self, sliceID, seed):
      if seed == None:
         if len(self.coordinates) == 0:
            print "A seed point required to initiate the vessel."
            return False
         else:
            seed = self.getPreviousCoordinate()[1]
      if not IO.hasImage(sliceID):
         print "slice ID out of bounds."
         return False
      flag = True
      while(flag):
         flag = self.extendVessel(sliceID, seed)
         # TODO: take care of case where no coordinates?
         seed = self.getPreviousCoordinate()[1]
         sliceID += 1
         if not IO.hasImage(sliceID):
            self.endVessel("out-of-scope")
            self.clearVariables()
            return True
      return True
   def exploreToBeginning(self, sliceID, seed):
      #if seed:
      #   print "exploreToBeginning: " + str(sliceID) + ": (" + str(seed[0]) + ", " + str(seed[1]) + ")"
      #else:
      #   print "exploreToBeginning: " + str(sliceID)
      if seed == None:
         if len(self.coordinates) == 0:
            print "A seed point required to initiate the vessel."
            return False
         else:
            seed = self.getPreviousCoordinate()[1]
      if not IO.hasImage(sliceID):
         print "slice ID out of bounds."
         return False
      flag = True
      self.setExploreInReverse(True)
      while(flag):
         flag = self.extendVessel(sliceID, seed)
         # TODO: take care of case where no coordinates?
         seed = self.getPreviousCoordinate()[1]
         sliceID -= 1
         if not IO.hasImage(sliceID):
            self.beginVessel("out-of-scope")
            self.clearVariables()
            return True
      return True
   def printVessel(self):
      print "Vessel ID: " + str(self.vesselID)
      print "Beginning explored: " + str(self.beginningExplored)
      print "Beginning reason: " + str(self.beginningReason)
      print "End explored: " + str(self.endExplored)
      print "End reason: " + str(self.endReason)
      print "Length of vessel: " + str(len(self.coordinates)) + " slices"
      print "Range of vessel: " + str(self.coordinates[0][0]) + " - " + \
            str(self.coordinates[-1][0])
      print "Slice ID\tCentroid\t\tCross section area"
      for coordinate in self.coordinates:
         print str(coordinate[0]).rjust(3) + "\t\t(" + str("%.2f" % coordinate[1][0]).rjust(6) + ", " +\
               str("%.2f" % coordinate[1][1]).rjust(6) + ")\t" + str(coordinate[2]).rjust(3)
      print "Next vessels:"
      for next in self.next:
         print "\tVessel " + str(next.vesselID),
         if next.isExplored():
            print " - Explored"
         else:
            print " - To be explored"
      if not self.next:
         print "\tNone"
      print "Previous vessels:"
      for prev in self.previous:
         print "\tVessel " + str(prev.vesselID),
         if prev.isExplored():
            print " - Explored"
         else:
            print " - To be explored"
      if not self.previous:
         print "\tNone"
      print
