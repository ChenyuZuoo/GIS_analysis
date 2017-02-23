# -*- coding: utf-8 -*- 
################################################
### GIS Application Project : 3D Convex Hull ###
################################################

#https://github.com/diwi/QuickHull-3D/blob/master/src/MAIN_quickhull/DwConvexHull3D.java

# add Edge
# add findEdge(edgeList):
# add getMiddlePnt(pnt1, pnt2):
# add getInnerPnt(face, pnt):


import arcpy
import numpy as np
import math
from __future__ import division

EPS = 1e-8

Points = [[0,0,0],[3,1,0],[0,3,0],[2,0,0],[0,4,6],[-2,0,7],[0,6,9],[1,-5,3]]

#############################################################
###########################Define Class######################
#############################################################
class Point(object):
  def __init__(self, x, y, z):
    self.x = x
    self.y = y
    self.z = z
    self.flag = False
    
class Edge(object):
  def __init__(self, p1, p2):
    self.p1 = p1
    self.p2 = p2
    # self.flag = False

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

### get normal vector 已知三点坐标，求法向量

# def getNormal(p1,p2,p3):
#   a = ( (p2.y-p1.y)*(p3.z-p1.z)-(p2.z-p1.z)*(p3.y-p1.y) )
#   b = ( (p2.z-p1.z)*(p3.x-p1.x)-(p2.x-p1.x)*(p3.z-p1.z) )
#   c = ( (p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x) )
#   return Vec3(a,b,c);





#####################################################################
####
#####################################################################

### Distance from point to point 
def distPnt2Pnt(p1,p2):
  return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2] - p2[2])**2)


### Distance from point to line
### import point, start and end point of line
def distPnt2Line(p,l1,l2):
  lineVec = l1[0]-l2[0], l1[1]-l2[1],l1[2]-l2[2]
  pntVec = l1[0]-p[0], l1[1]-p[1],l1[2]-p[2]
  lineLength = distPnt2Pnt(l1,l2)
  lineVec2 = l1[0]-l2[0]/lineLength,l1[1]-l2[1]/lineLength, l1[2]-l2[2]/lineLength
  pntVec2 = l1[0]-p[0]/lineLength,l1[1]-p[0]/lineLength, l1[2]-p[0]/lineLength
  t = dot(line_unitvec, pnt_vec_scaled)    
  if t < 0.0:
    t = 0.0
  elif t > 1.0:
    t= 1.0
  nearest = scale(line_vec, t)
  dist = distance(nearest, pnt_vec)
  return dist

### Distance from point to face
def distPnt2Face(p, face):
  return f_abs(face.a*p.x + face.b*p.y + face.c*p.z + face.d)/sqrt(face.a*face.a + face.b*face.b + face.c*face.c)


### Check if points are noncollinear  
def isCollinear():
# check = True
# a = x[0][0]-x[1][0],x[0][1]-x[1][1],x[0][2]-x[1][2]
# b = x[2][0]-x[1][0],x[2][1]-x[1][1],x[2][2]-x[1][2]
#  if (a[0]/b[0] == a[1]/b[1]):
#    if (a[2]/b[2] == a[1]/b[1]):
#      check = False
#      break
#  return check
# if ((2*((x[0][0]-x[2][0])*(x[1][0]-x[0][0])+(x[0][1]-x[2][1])*
#         (x[1][1]-x[0][1])+(x[0][2]-x[2][2])*(x[1][2]-x[0][2])))**2 
#     -4*((x[0][0]-x[2][0])**2+(x[0][1]-x[2][1])**2+(x[0][2]-x[2][2])**2)*
#     ((x[1][0]-x[0][0])**2+(x[1][1]-x[0][1])**2+(x[1][2]-x[0][2])**2)==0) 
#   flag = True
  flag = False 
  if distPnt2Line(p,l1,l2) < EPS:
    flag = True  
  return (flag) 
     
### Check if points are noncoplanar
def isCoplaner(p, face):
  flag = False
  if disPoint2Face(p, face) < EPS : 
    flag = True
  return (flag)
  
### Find Extreme Points that have maximum or minimum x, y, z
def extremePoints(pointSet):
  # data = []
  mmPoints = []
  max_x = -1
  max_y = -1
  max_z = -1
  min_x = 999999999
  min_y = 999999999
  min_z = 999999999
  
  for i in range(0,len(pointSet)):
    if (pointSet[i][0] < min_x):
      min_x = pointSet[i][0]
      index_xmin = i
    if (pointSet[i][0] > max_x):
      max_x = pointSet[i][0]
      index_xmax = i
    if (pointSet[i][1] < min_y):
      min_y = pointSet[i][1]
      index_ymin = i
    if (pointSet[i][1] > max_y):
      max_y = pointSet[i][1]
      index_ymax = i
    if (pointSet[i][2] < min_z):
      min_z = pointSet[i][2]
      index_zmin = i
    if (pointSet[i][2] > max_z):
      max_z = pointSet[i][2]
      index_zmax = i
      
  mmPoints = [pointSet[index_xmin], pointSet[index_xmax], pointSet[index_ymin], pointSet[index_ymax], pointSet[index_zmin], pointSet[index_zmax]]
  return mmPoints
  
def furthestPnt2Pnt(mmPoints):
  ### Build up a basic triangle of the tetrahedron
  ### Fisrt, find two most distant points to build up a line of the triangle
  maxP2PDist = -1
  
  ### Check if there are more than two max-min points in the dataset   
  if len(mmPoints) > 2 :
    for i in range(0,len(mmPoints)):
      for j in range(i+1, len(mmPoints)):
        # Calculate distance point to point 
        p2PDist = math.sqrt((mmPoints[i][0]-mmPoints[j][0])**2 + (mmPoints[i][1]-mmPoints[j][1])**2 + (mmPoints[i][2] - mmPoints[j][2])**2)
        if (p2PDist > maxP2PDist):
          maxP2PDist = p2PDist
          index1 = i
          index2 = j
  
  ### in extreme cases, there might be only two extremes points among the data set, with a shape like spindle
  else:
    index1 = 0
    index2 = 1
  return index1, index2  


### Find a distant point to the line
# inFC : mmPoints
def furthestPnt2Line (inFC,index1,index2):
  maxP2LDist = -1
  for i in range(0, len(inFC)):
    p2LDist = distPnt2Line(inFC[i],inFC[index1],inFC[index2])
    if (p2LDist>maxP2LDist):
      maxP2LDist = p2LDist
      index3 = i
  return index3


### Find the most distant point from a face
def furthestPnt2Face(pntSet,face):
  dist = 0
  for point in pntSet:
    dist0 = disPnt2Face(point, face)
    if dist0 > dist:
      dist = dist0
      disPoint = point
  return point


  for j in range(0, len(pntSet)):
    dist0 = disPnt2Face(pntSet[j], face)
    if dist0 > dist:
      dist = dist0
      index_dist = j
    #if (dist == 0) :
      #print("All the points from the point clouds are Coplaner.")
      #break  
  return pntSet[index_dist]
      
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
        if compareEdge(edgeList[i], edgeList[j]):
          index.append(i)
          index.append(j)
    
    print index
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

### Check if the point can see the faces 
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
  
  if vis > 0.5:
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
            tempPointMap[face].append(point)
          else:
            tempPointMap[face] = [point]
          break
    return tempPointMap

# get light faces of a point from all the active faces
def getLightFaces(pnt, faceSet):
  tempSet = []
  for face in faceSet:
    if checkVisibility(face, pnt):
      tempSet.append(face)
  return tempSet

# ignore face without points
def faceSetUpdate(faceSet, hullSet, faceMap):
    tempFaceSet = []
    for face in faceSet:
      if isFinalFace(face,faceMap):
        hullSet.append(face)
      else:
        tempFaceSet.append(face)
    return tempFaceSet

# check if a face is final hull
def isFinalFace(face, mPoint2FaceMap):
  f = False
  if face.flag:
    f = True
    return f
  else:
    if not mPoint2FaceMap.has_key(face):
      f = True
      return f
  return f

# check if a point is inside convex hull
def isInnerPoint(point, mFaceSet):
  f = True
  if point.flag:
    return f
  else:
    for face in mFaceSet:
      if (checkVisibility(face, point) == -1):
        f = False
        return f
  return f

# update active point list
def pointSetUpdate(mPointsSet): 
  tempSet = []
  for point in mPointsSet:
    if not isInnerPoint:
      tempSet.append(point)
  return tempSet

# get all the outside points belong to a face
def getPoints(faceMap, face):
  points = faceMap[face]
  return points

# detect the furest point of evert active faces, and build some new faces
def expandHull(faceSet, mPoint2FaceMap, pointObjSet):
  for face in faceSet:
    # 1. find the furthest point in the set of this face
    # 2. 检测出该点所有可见的面
    # 3. 建立新的face  调用Point()构建新的面
    # 4. deactive light face
    # 5. 更新mFaceSet


###############################################
###############Initialization##################
###############################################
def init(inFC):
  
  # 1) Read points into pntsSet
  pntsSet = readPoints(inFC)
  
  # 2） compute 4 points, for the initial hull 
  extPntsSet = extremePoints(pntsSet) # find extreme points in x, y, z (maximum and minimum)
  i = furthestPnt2Pnt(extPntsSet[0]) # find two most distant points among extreme points to build up a line segment 
  j = furthestPnt2Pnt(extPntsSet[1])
  pn1 = extPntsSet[i]
  pn2 = extPntsSet[j]
  m = furthestPnt2Line (extPntsSet,i,j) # find the most distant point to the line segment
  pn3 = extPntsSet[m]
  face0 = Face(pn1, pn2, pn3) # build up a basic triangle for the initial hull
  pn4 = furthestPnt2Face(extPntsSet,face0) # find the most distant point to the triangle 
  face1 = Face(pn1, pn2, pn4) # build up the other three faces of the initial hull 
  face2 = Face(pn1, pn3, pn4)
  face3 = Face(pn2, pn3, pn4)
  
  # 1. 读数据
  # 2. 建立初始四面体
  # 3. 分配所有点到临近面
  # 4. 把每个面放入mFaceSet作为激活状态的标记

###############################################
##################Iteration####################
###############################################
while (mFaceSet):
  expandHull()


### active data set
mFaceSet = [] # temp hull
mHullTriSet = [] # identified as final hull facet element 
mPointSet = readPoints(tempFC) # Points outside of temp hull
# dic to restore all the assignment of points to face
mPoint2FaceMap = {} # belonging status of active points and active faces

