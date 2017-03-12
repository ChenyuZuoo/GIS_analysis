# -*- coding: utf-8 -*- 
from __future__ import division
import arcpy
import numpy as np
import math


# Read mFaceSet data into coordList 
coordList = []

for i in range(0, len(mFaceSet)):
  coordList[i][0][0] = Face[i].p1.x
  coordList[i][0][1] = Face[i].p1.y
  coordList[i][0][2] = Face[i].p1.z
  coordList[i][1][0] = Face[i].p2.x
  coordList[i][1][1] = Face[i].p2.y
  coordList[i][1][2] = Face[i].p2.z
  coordList[i][2][0] = Face[i].p3.x
  coordList[i][2][1] = Face[i].p3.y
  coordList[i][2][2] = Face[i].p3.z

# Create empty Point and Array objects
point = arcpy.Point()
array = arcpy.Array()

triaList = []

for tria in coordList:
    # For each coordinate pair, set the x,y,z properties and add to the Array object
    for coordPair in tria:
        point.X = coordPair[0]
        point.Y = coordPair[1]
        point.Z = coordPair[2]
        array.add(point)

    # Add the first point of the array in to close off the polygon
    array.add(array.getObject(0))

    # Create a Polygon object based on the array of points
    polygon = arcpy.Polygon(array,None,True,True)

    # Clear the array for future use
    array.removeAll()

    # Append to the list of Polygon objects
    triaList.append(polygon)

arcpy.CopyFeatures_management(triaList, OutFC)