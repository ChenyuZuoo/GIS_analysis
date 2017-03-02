# -*- coding: utf-8 -*-
from __future__ import division
# import arcpy
import numpy as np
import math


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
        self.a = ((p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y))
        self.b = ((p2.z - p1.z) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.z - p1.z))
        self.c = ((p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x))
        self.d = (0 - (self.a * p1.x + self.b * p1.y + self.c * p1.z))
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

        self.flag = False;


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

L1 = Edge(p1, p2)
L2 = Edge(p2, p3)
# L1_rev = Edge(p2, p1)
# L1_rep = Edge(p1, p2)
L3 = Edge(p3, p4)

f1 = Face(p1, p2, p3)
f2 = Face(p1, p2, p6)
f3 = Face(p2, p3, p6)
f4 = Face(p1, p3, p6)


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

    vis = ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)
    # print vis

    # if vis > 0.5:
    ### 0.5 as threshold is too bigh here, so change it to strict 0
    ### this value depends on the points scale and distance
    if vis > 0:
        return True  # unvisible
    else:
        return False  # visible


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


def faceFactory(faceP1, faceP2, faceP3, innerPnt):
    face = Face(faceP1, faceP2, faceP3)
    if not checkVisibility(face, innerPnt):
        return face
    else:
        return Face(faceP1, faceP3, faceP2)


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

### Distance from point to face
def distPnt2Face(p, face):
  return math.fabs(face.a*p.x + face.b*p.y + face.c*p.z + face.d)/math.sqrt(face.a*face.a + face.b*face.b + face.c*face.c)

# dis = distPnt2Face(p9, f2)
# print dis

### Find the most distant point from a face
def furthestPnt2Face(pntSet,face):
  dist = 0
  # for pnt in pntSet:
  #   dist = distPnt2Face(pnt, face)
  #   if dist > dist0:
  #     dist0 = dist
  #     temp = pnt
  # return temp
  for j in range(0, len(pntSet)):
    dist0 = distPnt2Face(pntSet[j], face)
    if dist0 > dist:
      dist = dist0
      index4 = j
  # print ("index4")
  # print (index4)
  return pntSet[index4]

# testp = [p9, p10,p14]
# print testp
# pnt = furthestPnt2Face(testp, f2)
# print pnt.__dict__
# get all the outside points belong to a face
def getPoints(faceMap, face):
  # temp = []
  points = faceMap[face]
  # for pnt in points:
  #   temp.append(pnt)
  return points

# get light faces of a point from all the active faces
def getLightFaces(pnt, faceSet):
  tempSet = []
  for face in faceSet:
    if checkVisibility(face, pnt):
      tempSet.append(face)
  return tempSet

def getOutRing(edgeList):
    temp = []
    index = []
    # edgeList = []
    # for face in faceSet:
    #   edgeList.append(face.e1)
    #   edgeList.append(face.e2)
    #   edgeList.append(face.e3)

    for i in range(len(edgeList)):
      for j in range(i+1, len(edgeList)):
        if edgeList[i].__dict__ == edgeList[j].__dict__:
          index.append(i)
          index.append(j)
    
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

# def buildFaces


p = getMiddlePnt(p1, p2)
p = getInnerPnt(f1, p6)
# print p.__dict__
f1 = faceFactory(f1.p1, f1.p2, f1.p3, p)
f2 = faceFactory(f2.p1, f2.p2, f2.p3, p)
f3 = faceFactory(f3.p1, f3.p2, f3.p3, p)
f4 = faceFactory(f4.p1, f4.p2, f4.p3, p)
mPointSet = [p9, p10, p11, p12,p, p13]

mFaceSet = [f1, f2, f3, f4]

mMap = pointsAssignment(mFaceSet, mPointSet)

lightface = getLightFaces(p10, mFaceSet)
# print lightface
edgeset = getEdgeSet(lightface)
# print edgeset
outring = getOutRing(edgeset)
# print outring
print 'mFaceSet'
print mFaceSet
print 'lightface'
print lightface
for face in lightface:
  mFaceSet.remove(face)
print '**************'
print mFaceSet



