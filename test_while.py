# -*- coding: utf-8 -*- 
from __future__ import division
# import arcpy
import numpy as np
import math

EPS = 1e-8

# Points = [[0,0,0],[3,1,0],[0,3,0],[2,0,0],[0,4,6],[-2,0,7],[0,6,9],[1,-5,3]]

#############################################################
#################### Define Class ###########################
#############################################################
class Point(object):
  def __init__(self, x, y, z):
    self.x = x
    self.y = y
    self.z = z
    self.flag = False

class Edge(object):
  def __init__(self,p1,p2):
    self.p1 = p1
    self.p2 = p2
    self.flag = False

class Face(object):
  def __init__(self, p1, p2, p3):
    self.p1 = p1
    self.p2 = p2
    self.p3 = p3
    ### surface varibles ax+by+cz+d=0
    self.a = ( (p2.y-p1.y)*(p3.z-p1.z)-(p2.z-p1.z)*(p3.y-p1.y) )
    self.b = ( (p2.z-p1.z)*(p3.x-p1.x)-(p2.x-p1.x)*(p3.z-p1.z) )
    self.c = ( (p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x) )
    self.d = ( 0-(self.a*p1.x+self.b*p1.y+self.c*p1.z))
    ### is inside hull
    if self.a < 0:
      self.a = -self.a
      self.b = -self.b
      self.c = -self.c
      self.d = -self.d
    # edges in this face
    self.e1 = Edge(p1, p2)
    self.e2 = Edge(p2, p3)
    self.e3 = Edge(p3, p1)

    self.flag = False
  

p0 = Point(0, 0, 0)
p1 = Point(1, 0, 0)
p2 = Point(1, 1, 0)
p3 = Point(0, 1, 0)
p4 = Point(0, 0, 1)
p5 = Point(1, 0, 1)
p6 = Point(1, 1, 1)
p7 = Point(0, 1, 1)
p8 = Point(0.5, 0.5, 0.5)
p9 = Point(2, 0.9, 0.5)
p10 = Point(2, 2, 0.5)
p11 = Point(0.5, 0.5, 1)
p12 = Point(0.5, 0.5, -1)
p13 = Point(2, 0.8, 0.4)
p14 = Point(3, 0.8, 0.4)
p15 = Point(3,3,0)

L1 = Edge(p1, p2)
L2 = Edge(p2, p3)
# L1_rev = Edge(p2, p1)
# L1_rep = Edge(p1, p2)
L3 = Edge(p3, p4)

f1 = Face(p1, p2, p3)
f2 = Face(p1, p2, p6)
f3 = Face(p2, p3, p6)
f4 = Face(p1, p3, p6)

# read all the points info into memory
def readPoints(inFC):
    mPointsSet = []
    fieldList = ["SHAPE@XYZ"]
    with arcpy.da.SearchCursor(inFC, fieldList) as cur:
        for row in cur:
          point = Point(row[0][0],row[0][1],row[0][2])
          mPointsSet.append(point)
    del cur
    return mPointsSet


#####################################################################
############### basic scitific computation functions ################
#####################################################################

### Distance from point to point 
def distPnt2Pnt(p1,p2):
  return math.sqrt((p1.x-p2.x)**2 + (p1.y-p2.y)**2 + (p1.z - p2.z)**2)


### Distance from point to line
### import point, start and end point of line
def distPnt2Line(p,l1,l2):
  lineVec = l1.x-l2.x, l1.y-l2.y,l1.z-l2.z
  pntVec = l1.x-p.x, l1.y-p.y,l1.z-p.z
  lineLen = distPnt2Pnt(l1,l2)
  # check if the distance is the perpendicular length from the point to the line 
  # or is the distance from the point to the start or end of the line 
  t = (l1.x-l2.x)/lineLen*(l1.x-p.x)/lineLen + (l1.y-l2.y)/lineLen*(l1.y-p.y)/lineLen + (l1.z-l2.z)/lineLen*(l1.z-p.z)/lineLen
  if t < 0.0:
    t = 0.0
  elif t > 1.0:
    t= 1.0
  nearest = (l1.x-l2.x)*t, (l1.y-l2.y)*t, (l1.z-l2.z)*t 
  dist = math.sqrt((nearest[0]-pntVec[0])**2 + (nearest[1]-pntVec[1])**2 + (nearest[2] - pntVec[2])**2)
  return dist


### Distance from point to face
def distPnt2Face(p, face):
  return math.fabs(face.a*p.x + face.b*p.y + face.c*p.z + face.d)/math.sqrt(face.a*face.a + face.b*face.b + face.c*face.c)


### Check if points are noncollinear  
def isNonCollinear(p0,p1,p2):
  flag = True
  if ((2*((p0.x-p2.x)*(p1.x-p0.x)+(p0.y-p2.y)*(p1.y-p0.y)+(p0.z-p2.z)*(p1.z-p0.z)))**2 -4*((p0.x-p2.x)**2+(p0.y-p2.y)**2+(p0.z-p2.z)**2)*((p1.x-p0.x)**2+(p1.y-p0.y)**2+(p1.z-p0.z)**2)==0):
    flag = False
  return flag
     
### Check if points are noncoplanar
def isNonCoplaner(p, face):
  flag = True
  if disPoint2Face(p, face) < EPS : 
    flag = False
  return (flag)

### find extreme points in the direction
def extremePoints(pointSet):
    # data = []
    mmPoints = []
    max_x = -1
    max_y = -1
    max_z = -1
    min_x = 999999999
    min_y = 999999999
    min_z = 999999999

    for i in range(0, len(pointSet)):
        # print pointSet[i].x
        if (pointSet[i].x < min_x):
            min_x = pointSet[i].x
            index_xmin = i
        if (pointSet[i].x > max_x):
            max_x = pointSet[i].x
            index_xmax = i
        if (pointSet[i].y < min_y):
            min_y = pointSet[i].y
            index_ymin = i
        if (pointSet[i].y > max_y):
            max_y = pointSet[i].y
            index_ymax = i
        if (pointSet[i].z < min_z):
            min_z = pointSet[i].z
            index_zmin = i
        if (pointSet[i].z > max_z):
            max_z = pointSet[i].z
            index_zmax = i

    mmPoints = [pointSet[index_xmin], pointSet[index_xmax], pointSet[index_ymin], pointSet[index_ymax], pointSet[
        index_zmin], pointSet[index_zmax]]
    return mmPoints
  
def furthestPnt2Pnt(pntSet):
  ### Build up a basic triangle of the tetrahedron
  ### Fisrt, find two most distant points to build up a line of the triangle
  maxP2PDist = -1
  
  ### Check if there are more than two max-min points in the dataset   
  if len(pntSet) > 2 :
    for i in range(0,len(pntSet)):
      for j in range(i+1, len(pntSet)):
        # Calculate distance point to point 
        p2PDist = math.sqrt((pntSet[i].x-pntSet[j].x)**2 + (pntSet[i].y-pntSet[j].y)**2 + (pntSet[i].z - pntSet[j].z)**2)
        if (p2PDist > maxP2PDist):
          maxP2PDist = p2PDist
          index1 = i
          index2 = j
  
  ### in extreme cases, there might be only two extremes points among the data set, with a shape like spindle
  else:
    index1 = 0
    index2 = 1
#   print ("Index1 and index2")
#   print (index1, index2)
  return pntSet[index1], pntSet[index2]

### Find a distant point to the line
# inFC : mmPoints
def furthestPnt2Line (pntSet,p1,p2):
  maxP2LDist = -1
  for k in range(0, len(pntSet)):
    p2LDist = distPnt2Line(pntSet[k],p1,p2)
    if isNonCollinear(pntSet[k],p1,p2):
      if (p2LDist>maxP2LDist):
        maxP2LDist = p2LDist
        index3 = k
#   print ("index3")
#   print (index3)
  return pntSet[index3]

### Find the most distant point from a face
def furthestPnt2Face(pntSet,face):
  dist = 0
  for j in range(0, len(pntSet)):
    dist0 = distPnt2Face(pntSet[j], face)
    if dist0 > dist:
      dist = dist0
      index4 = j
    #if (dist == 0) :
      #print("All the points from the point clouds are Coplaner.")
      #break  
#   print ("index4")
#   print (index4)
  return pntSet[index4]    
      

  
def getMiddlePnt(pnt1, pnt2):
  x = (pnt1.x + pnt2.x) / 2
  y = (pnt1.y + pnt2.y) / 2
  z = (pnt1.z + pnt2.z) / 2
  return Point(x, y, z)

# this func used for find the inner point of the initial tetrahedron
def getInnerPnt(face, pnt):
  p = getMiddlePnt(getMiddlePnt(getMiddlePnt(face.p1, face.p2), face.p3), pnt)
  p.flag = True
  return p  
  
  
#####################################################################
################# space logic computation functions #################
#####################################################################

# construct a face with right-handed rule
# input is 4 points, output is a face build up with first three points
# rPnt is a point inside the init hull, so the volume with any face would be negative
def faceFactory(faceP1, faceP2, faceP3, innerPnt):
  face = Face(faceP1, faceP2, faceP3)
  if not checkVisibility(face, innerPnt):
    return face
  else:
    return Face(faceP1, faceP3, faceP2)

# detecting whether two edges are the same
def isSameEdge(edge1, edge2):
  if (edge1.p1.x == edge2.p1.x and edge1.p1.y == edge2.p1.y and edge1.p1.z == edge2.p1.z \
    and edge1.p2.x == edge2.p2.x and edge1.p2.y == edge2.p2.y and edge1.p2.z == edge2.p2.z) \
  or (edge1.p1.x == edge2.p2.x and edge1.p1.y == edge2.p2.y and edge1.p1.z == edge2.p2.z \
    and edge1.p2.x == edge2.p1.x and edge1.p2.y == edge2.p1.y and edge1.p2.z == edge2.p1.z):
    return True
    
### if there's more than one light face, a outer ring is needed to be find
### then using the points on this ring to constract new faces
def getOutRing(edgeList):
    temp = []
    index = []

    for i in range(len(edgeList)):
      for j in range(i+1, len(edgeList)):
        if isSameEdge(edgeList[i], edgeList[j]):
          index.append(i)
          index.append(j)
    
    # print index
    for i in range(len(edgeList)):
      if i not in index:
        temp.append(edgeList[i])

    return temp

# get all the edges of a face set
def getEdgeSet(faceSet):
  temp = []
  for face in faceSet:
    temp.append(face.e1)
    temp.append(face.e2)
    temp.append(face.e3)
  return temp  
  
##########################################
  
  

### Check if the point can see the faces 
### right-hand rule is crutial here
def checkVisibility(face, p):
  
  ax = face.p1.x - p.x
  ay = face.p1.y - p.y
  az = face.p1.z - p.z
  bx = face.p2.x - p.x
  by = face.p2.y - p.y
  bz = face.p2.z - p.z
  cx = face.p3.x - p.x
  cy = face.p3.y - p.y
  cz = face.p3.z - p.z
  
  vis = ax *(by*cz - bz*cy) + ay * (bz*cx - bx*cz)+ az * (bx*cy - by*cx)
  
  # if (vis > 0.5) :
  #   return 1  ### unvisible(inside face)
  # elif (vis<-0.5) :
  #   return -1  ### visible(outside of face)
  # else:
  #   return 0  ### coplannar with hull
  
  if vis > 0:
    return True
  else:
    return False

### Assign each point to a face
def pointsAssignment(faceSet, pointSet):
    tempPointMap = {}
    for point in pointSet:
        for face in faceSet:
            if checkVisibility(face, point):
                if face in tempPointMap:
                    # tempPointMap[face].append(point.__dict__)
                    tempPointMap[face].append(point)
                else:
                    # tempPointMap[face] = [point.__dict__]
                    tempPointMap[face] = [point]
                break
    return tempPointMap
  
# get light faces of a point from all the active faces
def getLightFaces(pnt, faceSet):
  tempSet = []
  for face in faceSet:
    if checkVisibility(face, pnt):
      face.flag = True
      tempSet.append(face)
  return tempSet

# ignore face without points
# def faceSetUpdate(faceSet, hullSet, faceMap):
#     tempFaceSet = []
#     for face in faceSet:
#       if isFinalFace(face,faceMap):
#         hullSet.append(face)
#       else:
#         tempFaceSet.append(face)
#     return tempFaceSet

def lightFacesRem(faceSet):
    tempFaceSet = []
    for face in faceSet:
      if not face.flag:
        tempFaceSet.append(face)
    return tempFaceSet

# check if a face is final hul
# def isFinalFace(face, mPoint2FaceMap):
#   f = False
#   if face.flag:
#     f = True
#     return f
#   else:
#     if not mPoint2FaceMap.has_key(face):
#       f = True
#       return f
#   return f

# check if a point is inside convex hull
def isInnerPoint(point, mFaceSet):
  f = True
  if point.flag:
    return f
  else:
    for face in mFaceSet:
      if checkVisibility(face, point):
        f = False
        return f
  return f

# update active point list
def pointSetUpdate(pointSet, faceSet): 
  tempSet = []
  for point in pointSet:
    if not isInnerPoint(point, faceSet):
      tempSet.append(point)
  return tempSet

# get all the outside points belong to a face
def getPoints(faceMap, face):
  points = faceMap[face]
  # print faceMap.keys()
  return points

# detect the furest point of evert active faces, and build some new faces
def expandHull(faceSet, faceMap, pointSet, innerPnt):
  for face in faceSet:
    if faceMap.has_key(face):
      # 1. find the furthest point in the set of this face
      pointOfFace = getPoints(faceMap, face)
      furpnt = furthestPnt2Face(pointOfFace, face)

      # 2. get all the (light faces) visible faces for this point
      lightFaces = getLightFaces(furpnt, faceSet)

      # 3. construct new faces with furest point and light faces
      edgeset = getEdgeSet(lightFaces)
      outring = getOutRing(edgeset)
      # newFaces = []
      for edge in outring:
        newFace = faceFactory(edge.p1, edge.p2, furpnt, innerPnt)
        faceSet.append(newFace)

      # 4. deactive light face (delete light faces from faceSet)
      # for lightface in lightFaces:
      #   faceSet.remove(lightface)
      faceSet = lightFacesRem(faceSet)

  return faceSet



### active data set
mFaceSet = [] # temp hull
# mHullTriSet = [] # identified as final hull facet element 
mPointSet = [p0, p1, p2, p3, p4, p5, p6, p7]#, p8, p9, p10, p11, p12, p13, p14, p15] # readPoints(tempFC) # Points outside of temp hull
# dic to restore all the assignment of points to face
mPoint2FaceMap = {} # belonging status of active points and active faces

extPntsSet = extremePoints(mPointSet)
pnt1 = furthestPnt2Pnt(extPntsSet)[0]
pnt2 = furthestPnt2Pnt(extPntsSet)[1]
pnt3 = furthestPnt2Line (extPntsSet,pnt1,pnt2)
face1 = Face(pnt1, pnt2, pnt3) # build up a basic triangle for the initial hull
pnt4 = furthestPnt2Face(extPntsSet,face1)
InnerPnt = getInnerPnt(face1, pnt4)
# re-initialise faces stored in a right  hand way
face1 = faceFactory(pnt1, pnt2, pnt3, InnerPnt)
face2 = faceFactory(pnt1, pnt2, pnt4, InnerPnt)
face3 = faceFactory(pnt1, pnt3, pnt4, InnerPnt)
face4 = faceFactory(pnt2, pnt3, pnt4, InnerPnt)

mFaceSet = [face1, face2, face3, face4]
mPointSet = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15] 
mPoint2FaceMap = pointsAssignment(mFaceSet, mPointSet)

# keys = mPoint2FaceMap.keys()

# print mPoint2FaceMap

###############################################
################### init ######################
###############################################
# mPointSet, mFaceSet, mPoint2FaceMap, innerPnt = init(inFC)

###############################################
##################Iteration####################
###############################################
# main function, using global varibles
while (mPoint2FaceMap):
  expandHull(mFaceSet, mPoint2FaceMap, mPointSet, InnerPnt)
  # point set update 
  mPointSet = pointSetUpdate(mPointSet, mFaceSet)
  # faceMap update 
  mPoint2FaceMap = pointsAssignment(mFaceSet, mPointSet)
  # face set update
  # mFaceSet = faceSetUpdate(mFaceSet, hullSet, mPoint2FaceMap)


