from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random
import math
import numpy
import scipy

class Segment:
   def __init__(self, origin, end):
      self.origin = origin
      self.end = end
      self.length = int(math.fabs(end[3] - origin[3]))
      self.coordinates = []
      c1x = origin[0]
      c1y = origin[1]
      c1r = origin[2]
      c1z = origin[3]
      c2x = end[0]      
      c2y = end[1]
      c2r = end[2]
      c2z = end[3]
      
      xincrement = float(c2x - c1x)/float(self.length)
      yincrement = float(c2y - c1y)/float(self.length)
      rincrement = float(c2r - c1r)/float(self.length - 1)
      zincrement = float(c2z - c1z)/float(self.length)
      
      for i in range(self.length):
         x = float(c1x) + i*xincrement
         y = float(c1y) + i*yincrement
         r = float(c1r) + i*rincrement
         z = origin[3] + i*zincrement
         if i != 0 and i != self.length - 1:
            #temp = 0
            x = random.uniform(x - 1, x + 1)
            y = random.uniform(y - 1, y + 1)
            #r = random.uniform(r - 0.5, r + 0.5)
         self.coordinates.append((x, y, r, z))

   
   def getZLimits(self):
      # TODO: Right now, origin[2] should be less than end[2]
      return (self.origin[2], self.end[2])
      
class ImageSeries:
   def __init__(self, xLimit, yLimit, zLimit):
      self.xLimit = xLimit
      self.yLimit = yLimit
      self.zLimit = zLimit
      self.images = [0]*zLimit
      for i in range(zLimit):
         self.images[i] = numpy.zeros((xLimit, yLimit))
      self.segments = []
         
   def getNumImages(self):
      return len(self.images)
      
   def addSegment(self, segment):
      self.segments.append(segment)
      (z1, z2) = segment.getZLimits()
      if z2 < z1:
         temp = z1
         z1 = z2
         z2 = temp
      #print "Segment:"
      for (x, y, r, z) in segment.coordinates:
         x = int(x)
         y = int(y)
         z = int(z)
         #print "(" + str(x) + ", " + str(y) + ", " + str(r) + ", " + str(z) + ")"
         for i in range(x - int(r) - 1, x + int(r) + 1):
            for j in range(y - int(r) - 1, y + int(r) + 1):
               if i >= 0 and i < self.xLimit and j >= 0 and j < self.yLimit and\
                  (i - x)*(i - x) + (j - y)*(j - y) <= r*r:
                  self.images[z][i][j] = 255
                  
         
   def printAllSegments(self):
      i = 0
      for segment in self.segments:
         print "Segment: " + str(i)
         for coordinate in segment.coordinates:
            print "(" + str(coordinate[0]) + ", " + str(coordinate[1]) + ", " + str(coordinate[3]) + ") - " + str(coordinate[2])
         print
         i += 1
         
   def showAll(self):
      plt.figure(figsize=(4, 4))
      i = 0
      for image in self.images:
         plt.subplot(111)
         plt.imshow(image, cmap=plt.cm.gray)
         plt.axis('off')
         plt.title("Image: " + str(i), fontsize=20)   
         plt.show()
         i += 1
   def plotAll(self):
      fig = plt.figure()
      ax = fig.gca(projection='3d')
      colors = ['r', 'g', 'b', 'k', 'm', 'y']
      
      i = 0
      for segment in self.segments:
         for coordinate in segment.coordinates:
            ax.scatter(coordinate[0], coordinate[1], coordinate[3], zdir='z', \
                       c=colors[i % 6])
         i += 1
                       
      ax.legend()
      ax.set_xlim3d(0, self.xLimit)
      ax.set_ylim3d(0, self.yLimit)
      ax.set_zlim3d(0, self.zLimit)
      ax.set_xlabel('X Axis')
      ax.set_ylabel('Y Axis')
      ax.set_zlabel('Z Axis')
      plt.show()
   
   def saveAll(self):
      i = 1
      for image in self.images:
         scipy.misc.imsave(str(i) +'.jpg', image)
         i += 1
      
         
def generate(synthetic):
   points = []
   if synthetic == 1:
      iseries = ImageSeries(120, 100, 66)
      # Convoluted network
      # Vessel 1
      points.append((10, 10, 4, 0))    # Red
      points.append((30, 20, 5, 8))
      
      # Vessel 2
      points.append((30, 20, 3, 8))    # Green
      points.append((45, 20, 5, 21))
      points.append((45, 20, 5, 21))   # Blue
      points.append((50, 20, 4, 34))
      
      # Vessel 3
      points.append((30, 20, 4, 8))    # Black
      points.append((40, 50, 5, 15))
      
      # Vessel 4
      points.append((40, 50, 3, 15))   # Magenta
      points.append((30, 60, 6, 29))
      
      # Vessel 5
      points.append((40, 50, 4, 15))   # Yellow
      points.append((70, 70, 3, 30))
      
      # Vessel 6
      points.append((50, 80, 5, 8))
      points.append((65, 60, 5, 15))
      points.append((65, 60, 5, 15))
      points.append((65, 40, 5, 25))
      points.append((65, 40, 5, 25))   # Red
      points.append((50, 20, 3, 34))
      
      # Vessel 7
      points.append((50, 20, 5, 34))   # Green
      points.append((70, 10, 5, 42))
      points.append((70, 10, 5, 42))   # Blue
      points.append((90, 20, 5, 50))
      points.append((90, 20, 5, 50))   # Black
      points.append((100, 30, 5, 58))
      points.append((100, 30, 5, 58))  # Magenta
      points.append((100, 50, 5, 66))
      
      # Vessel 8
      #points.append((75, 40, 3, 4))    # Yellow
      #points.append((100, 50, 5, 8))
      
      # Vessel 9
      #points.append((75, 40, 4, 4))    # 
      #points.append((70, 60, 5, 0))
      
      # Vessel 10
      points.append((70, 70, 5, 30))   # Red
      points.append((75, 85, 4, 38))
      points.append((75, 85, 4, 38))   # Green
      points.append((73, 87, 5, 50))
      
      # Vessel 11
      points.append((110, 90, 3, 20))  # Blue
      points.append((70, 70, 4, 30))
      
      
      # Vessel 12
      #points.append((73, 87, 4, 50))   # Red
      #points.append((40, 85, 5, 64))
      
      # Vessel 13
      #points.append((73, 87, 3, 50))   # Green
      #points.append((75, 50, 5, 62))
      
      
      
   
   # Simple Split
   # Origin and end point for Vessel 1 (x, y, radius, z); z is the slice number
   #points.append((20, 90, 6, 0))
   #points.append((70, 55, 8, 16))
   
   # Origin and end point for Vessel 2
   #points.append((60, 54, 4.5, 16))
   #points.append((80, 30, 4.5, 24))
   
   # Origin and end point for Vessel 3
   #points.append((75, 59, 5, 16))
   #points.append((100, 95, 4, 24))
   
   
   # Simple Merge
   # Vessel 1
   #points.append((100, 95, 4, 0))   # Green
   #points.append((72, 53, 5, 8))
   # Vessel 2
   #points.append((90, 30, 4.5, 0))
   #points.append((70, 54, 5, 8))
   # Vessel 3
   #points.append((70, 55, 8, 8))
   #points.append((20, 90, 6, 24))
   
   
   for i in range(len(points)/2):
      iseries.addSegment(Segment(points[2*i], points[2*i + 1]))
   
   #iseries.showAll()
   print "Series generated"
   #iseries.showAll()
   #iseries.printAllSegments()
   iseries.saveAll()
   return iseries
   
#generate()
