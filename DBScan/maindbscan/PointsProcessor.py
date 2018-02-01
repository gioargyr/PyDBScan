'''
Created on Dec 19, 2017

@author: indiana
'''

import numpy
from osgeo import gdal
from shapely.geometry import LineString
from shapely.geometry import MultiPoint
from shapely.geometry import Point
from shapely.geometry import Polygon
import sys, os
import time


class PointsProcessor:
    
    def __init__(self, imageFilePath):
        
        # DB Scan constants/thresholds 
        self.pixelValueThreshold    = 1.0     #2.5
        self.eps                    = 0.00035
        self.minPts                 = 4
        
        self.allCorePoints      = []
        self.coreToReachables   = {}
        self.checkedCorePoints  = []
        self.corePointsInClusters   = []
        self.clusters           = []
        
        start = time.time()
        self.driver(imageFilePath)
        print("\nRuning time:\t" + str(time.time() - start))

    
    """
        The driver method runs the current class.
    """    
    def driver(self, inpFile = None):
        
        pixelPointsAboveThres = self.readImage(inpFile)
        print("pixelPointsAboveThres =\t" + str(len(pixelPointsAboveThres)))
        
        self.findCoreAndReachables(pixelPointsAboveThres)
        print("allCorePoints =\t" + str(len(self.allCorePoints)))
        print("coreToReachables =\t" + str(len(self.coreToReachables)))
        
        self.clusteringCorePoints()
        print("corePointsInClusters =\t" + str(len(self.corePointsInClusters)))
        
        self.addingReachablesToClusters()
        print("clusters =\t" + str(len(self.clusters)))
        
        self.printingClusters(os.path.dirname(inpFile))
        
        ## DEBUGGING EVERYTHING:
        i = 0
        n = 0
        onlyReachables = []
        clusterToReacheables = {}
        allPointsOfInterest = []
        for cl in self.clusters:
            i = i + 1
            n = n + len(cl)
            k = 0
            for point in cl:
                allPointsOfInterest.append(point)
                if point not in self.allCorePoints:
                    onlyReachables.append(point)
                    k = k + 1
            clusterToReacheables[str(i)] = onlyReachables
            onlyReachables = []
            
            print("CL" + str(i) + " has " + str(len(cl)) + " points.\t" + str(k) + " of them are plain reachables.")
                    
        
        print("all points in clusters =\t" + str(n))
        
        print("allPointsOfInterest are = " + str(len(allPointsOfInterest)))
#         multiAllPointsOfInterest = MultiPoint(allPointsOfInterest)
#         print("\nall points of interest:\n" + str(multiAllPointsOfInterest) + "\n")
        
#############################  END OF __init__  #################################

    """
        Find Core Points that are close to each Core Point (same eps as in Reachables)
        Use corePointToCloseCorePoints dictionary and create legitimate clusters of Core Points.
    """
    def printingClusters(self, outDir, printCHPolygons = True, printMPoints = False):
        
        clustersMultipoints = []
        clustersCHPolygons  = []
        
        for cluster in self.clusters:
            clustersMultipoints.append(MultiPoint(cluster))
        
        for mps in clustersMultipoints:
            clustersCHPolygons.append(mps.convex_hull)
              
        outFile = os.path.join(outDir, "DBScanClusters_new.txt")
        with open(outFile, "w") as out:
            
            if printCHPolygons:
                out.write("Clusters as Polygons(convex_hull):\n")
                for polygon in clustersCHPolygons:
                    out.write(str(polygon) + "\n")
            
            if printMPoints:
                out.write("Clusters as MultiPoints:\n")
                for mp in clustersMultipoints:
                    out.write(str(mp) + "\n")
            
            if not printCHPolygons and not printMPoints:
                print("No printing to output file defined!")
            
    
    """
        Find Core Points that are close to each Core Point (same eps as in Reachables)
        Use corePointToCloseCorePoints dictionary and create legitimate clusters of Core Points.
    """
    def addingReachablesToClusters(self):
        
        reachables = []
        completedCluster = []
        
        for cluster in self.corePointsInClusters:
            for corePoint in cluster:
                reachables = self.coreToReachables.get(str(corePoint))
                
                for point in reachables:
                    if point not in completedCluster:
                        completedCluster.append(point)
            
            self.clusters.append(completedCluster)
            completedCluster = []
                    

    """
        Find Core Points that are close to each Core Point (same eps as in Reachables)
        Use corePointToCloseCorePoints dictionary and create legitimate clusters of Core Points.
    """
    def clusteringCorePoints(self, kati = None):
        
        # At first find Core Points that are close to each Core Point.
        closeCorePoints = []
        corePointToCloseCorePoints = {}
        
        for corePoint in self.allCorePoints:
            for point in self.allCorePoints:
                if corePoint.distance(point) <= self.eps:
                    closeCorePoints.append(point)
            
            corePointToCloseCorePoints[str(corePoint)] = closeCorePoints
            closeCorePoints = []
            
        ## Uncomment if you want to view how many ReachablePoints and ReachableCorePoints each CorePoint has.
        #for p in listOfCorePoints:
        #    print(str(p) + "reachable points and reachable core-points:")
        #    print(str(len(coreToReachable[str(p)])))
        #    print(str(len(coreToCoreReachable[str(p)])))
        
        # Use the corePointToCloseCorePoints dictionary to create legit clusters of Core Points.
        legitClusterOfCorePoints = []
        i = 1
        for corePoint in self.allCorePoints:
            if corePoint not in self.checkedCorePoints:
                cluster = []
                closeCorePoints = []
                closeCorePoints = corePointToCloseCorePoints.get(str(corePoint))
                
                for closePoint in closeCorePoints:
                    if closePoint not in cluster:
                        cluster.append(closePoint)
                
                self.checkedCorePoints.append(corePoint)
                legitClusterOfCorePoints = self._checkCluster(cluster, corePointToCloseCorePoints)
                self.corePointsInClusters.append(legitClusterOfCorePoints)
                
                print(str(i) + "\t" + str(len(legitClusterOfCorePoints)) + "\t" + str(len(self.checkedCorePoints)))
                i = i + 1
    
    
    """
        checkCluster is checking the validation of a cluster of Core points
        Used recursively from clusteringCorePoints only.
    """
    def _checkCluster(self, inputCluster, corePointToCloseCorePoints):
        
        outputCluster   = []
        isContained     = True
        
        for point in inputCluster:
            if point not in self.checkedCorePoints:
                isContained = False
                break
        
        if isContained:
            return inputCluster
        
        else:
            outputCluster = list(inputCluster)  # practically, it copies one to another. Faster than copy!
            for clusteredPoint in inputCluster:
                if clusteredPoint not in self.checkedCorePoints:
                    closeCorePoints = []
                    closeCorePoints = corePointToCloseCorePoints.get(str(clusteredPoint))
                    for closePoint in closeCorePoints:
                        if closePoint not in outputCluster:
                            outputCluster.append(closePoint)
                    self.checkedCorePoints.append(clusteredPoint)
                    
            return self._checkCluster(outputCluster, corePointToCloseCorePoints)
            

    """
        Read image, filter pixels with values above threshold.
        (Return filtered pixels as list of points with geocoordinates.)
    """
    def readImage(self, imageFilePath):
        
        imageAsData     = gdal.Open(imageFilePath)
        pixelsAsNDArray = imageAsData.ReadAsArray()
        [cols, rows]    = pixelsAsNDArray.shape
        
        xoff, a, b, yoff, d, e = imageAsData.GetGeoTransform()
        
        pixelPointsAboveThres = []
        
        for i in range(rows):
            for j in range(cols):
                if pixelsAsNDArray[j][i] < -self.pixelValueThreshold or pixelsAsNDArray[j][i] > self.pixelValueThreshold:
                    x = a * i + b * j + xoff
                    y = d * i + e * j + yoff
                    pointInMap = Point(x, y)
                    pixelPointsAboveThres.append(pointInMap)
                    
        ## If you want to view all points in map:
#         multiPointsAboveThres = MultiPoint(pixelPointsAboveThres)
#         print("\nall (" + str(len(pixelPointsAboveThres)) + ") points above threshold:\n" + str(multiPointsAboveThres) + "\n")
                    
        return pixelPointsAboveThres


    """
        Check all points of interest(input list) and create a dictionary where:
        Key => Core Point, Value => Reachable Points
        (Return the dictionary)
    """
    def findCoreAndReachables(self, pointsOfInterest):
        
        reachablePoints = []
        
        for possibleCorePoint in pointsOfInterest:
            for point in pointsOfInterest:
                if possibleCorePoint.distance(point) <= self.eps:
                    reachablePoints.append(point)
            
            if len(reachablePoints) > self.minPts:
                self.allCorePoints.append(possibleCorePoint)
                self.coreToReachables[str(possibleCorePoint)] = reachablePoints
            
            reachablePoints = []
        
        ## Uncomment if you want to view how many ReachablePoints and ReachableCorePoints each CorePoint has.
        #for p in listOfCorePoints:
        #    print(str(p) + "reachable points and reachable core-points:")
        #    print(str(len(coreToReachable[str(p)])))
        #    print(str(len(coreToCoreReachable[str(p)])))
        


if __name__ == "__main__":
    if len(sys.argv) == 2:
        print("Going to read input from file: " + sys.argv[1])
        pointsProcessor = PointsProcessor(sys.argv[1]);
    else:
        print("No input file defined.")
        sys.exit(2)