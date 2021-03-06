'''
Created on Dec 19, 2017

@author: indiana
'''

from osgeo import gdal
#from shapely.geometry import LineString
from shapely.geometry import MultiPoint
from shapely.geometry import Point
from shapely.geometry import Polygon
import sys, os
import time


class PointsProcessor:
    
    def __init__(self, imageFilePath, outFileName, pixelValueThreshold = 3.1, eps = 0.000269, minPts = 5):
        
        # DB Scan constants/thresholds 
        self.pixelValueThreshold    = float(pixelValueThreshold)
        self.eps                    = float(eps)
        self.minPts                 = int(minPts)
        
        
        self.allCorePoints      = []
        self.coreToReachables   = {}
        self.checkedCorePoints  = []
        self.corePointsInClusters   = []
        self.clusters           = []
        
        self.loglike = os.path.join(os.path.dirname(imageFilePath), "log-like.txt")
        
        start = time.time()
        self.driver(outFileName, imageFilePath)
        #print("\nRunning time:\t" + str(time.time() - start))
        with open(self.loglike, "a") as out:
            out.write("\nRunning time:\t" + str(time.time() - start))

    
    """
        The driver method runs the current class.
    """    
    def driver(self, outFileName, inpFile = None):
        
        pixelPointsAboveThres = self.readImage(inpFile)
        with open(self.loglike, "a") as out:
            out.write("\npixelPointsAboveThres =\t" + str(len(pixelPointsAboveThres)))
        #print("\npixelPointsAboveThres =\t" + str(len(pixelPointsAboveThres)))
        
        self.findCoreAndReachables(pixelPointsAboveThres)
        with open(self.loglike, "a") as out:
            out.write("\nallCorePoints =\t" + str(len(self.allCorePoints)))
            out.write("\ncoreToReachables =\t" + str(len(self.coreToReachables)))
        #print("allCorePoints =\t" + str(len(self.allCorePoints)))
        #print("coreToReachables =\t" + str(len(self.coreToReachables)))
        
        self.clusteringCorePoints()
        with open(self.loglike, "a") as out:
            out.write("\ncorePointsInClusters =\t" + str(len(self.corePointsInClusters)))
        #print("corePointsInClusters =\t" + str(len(self.corePointsInClusters)))
        
        self.addingReachablesToClusters()
        with open(self.loglike, "a") as out:
            out.write("\nclusters =\t" + str(len(self.clusters)))
        #print("clusters =\t" + str(len(self.clusters)))
        
        self.printingClusters(os.path.dirname(inpFile), outFileName)
        
        ## DEBUGGING EVERYTHING:
        deb_time = time.time()
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
            
            with open(self.loglike, "a") as out:
                out.write("\nCL" + str(i) + " has " + str(len(cl)) + " points.\t" + str(k) + " of them are plain reachables.")
            #print("CL" + str(i) + " has " + str(len(cl)) + " points.\t" + str(k) + " of them are plain reachables.")
                    
        with open(self.loglike, "a") as out:
            out.write("\n\nall points in clusters =\t" + str(n))
            out.write("\nallPointsOfInterest are = " + str(len(allPointsOfInterest)))
            out.write("\n\ndeb-time = " + str(time.time() - deb_time))
        #print("all points in clusters =\t" + str(n))
        
        #print("allPointsOfInterest are = " + str(len(allPointsOfInterest)))
#         multiAllPointsOfInterest = MultiPoint(allPointsOfInterest)
#         print("\nall points of interest:\n" + str(multiAllPointsOfInterest) + "\n")
        
#############################  END OF __init__  #################################

    """
        Find Core Points that are close to each Core Point (same eps as in Reachables)
        Use corePointToCloseCorePoints dictionary and create legitimate clusters of Core Points.
    """
    def printingClusters(self, outDir, outFileName, printCHPolygons = True, printMPoints = False):
        
        clustersMultipoints = []
        clustersCHPolygons  = []
        
        for cluster in self.clusters:
            clustersMultipoints.append(MultiPoint(cluster))
        
        for mps in clustersMultipoints:
            clustersCHPolygons.append(mps.convex_hull)
              
        outFile = os.path.join(outDir, outFileName)
        with open(outFile, "w") as out:
            
            if printCHPolygons:
#                 out.write("Clusters as Polygons(convex_hull):\n")
                for polygon in clustersCHPolygons:
                    if isinstance(polygon, Polygon):
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
                
                with open(self.loglike, "a") as out:
                    out.write("\n" + str(i) + "\t" + str(len(legitClusterOfCorePoints)) + "\t" + str(len(self.checkedCorePoints)))
                #print(str(i) + "\t" + str(len(legitClusterOfCorePoints)) + "\t" + str(len(self.checkedCorePoints)))
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
           
        while len(pixelPointsAboveThres) < 400 and self.pixelValueThreshold > 1:
            pixelPointsAboveThres = []
            self.pixelValueThreshold = self.pixelValueThreshold - 0.1
            for i in range(rows):
                for j in range(cols):
                    if pixelsAsNDArray[j][i] < -self.pixelValueThreshold or pixelsAsNDArray[j][i] > self.pixelValueThreshold:
                        x = a * i + b * j + xoff
                        y = d * i + e * j + yoff
                        pointInMap = Point(x, y)
                        pixelPointsAboveThres.append(pointInMap)
            
            print(str(len(pixelPointsAboveThres)))
            print(str(self.pixelValueThreshold))
                    
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
                #print(str(len(self.allCorePoints)))
                self.coreToReachables[str(possibleCorePoint)] = reachablePoints
                #print(str(len(self.coreToReachables)))
            
            reachablePoints = []
        
        ## Uncomment if you want to view how many ReachablePoints and ReachableCorePoints each CorePoint has.
        #for p in listOfCorePoints:
        #    print(str(p) + "reachable points and reachable core-points:")
        #    print(str(len(coreToReachable[str(p)])))
        #    print(str(len(coreToCoreReachable[str(p)])))
        


if __name__ == "__main__":
    if len(sys.argv) == 3:
        print("Going to read input from file: " + sys.argv[1])
        print("Going to write output to file: " + sys.argv[2])
        print("No parameters for DBScan defined. Going to use defaults.")
        pointsProcessor = PointsProcessor(sys.argv[1], sys.argv[2]);
    elif len(sys.argv) == 6:
        print("Going to read input from file: " + sys.argv[1])
        print("Going to write output to file: " + sys.argv[2])
        print("User defined pixelValueThreshold =\t" + sys.argv[3])
        print("User defined eps =\t\t\t" + sys.argv[4])
        print("User defined minPts =\t\t\t" + sys.argv[5])
        pointsProcessor = PointsProcessor(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]);     
    else:
        print("No valid arguments.")
        sys.exit(2)