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
        self.pointliset = [p1, p2]


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

L1 = Edge(p1, p2)
L2 = Edge(p2, p3)
L1_rev = Edge(p2, p1)
L1_rep = Edge(p1, p2)
L3 = Edge(p3, p4)

f1 = Face(p1, p2, p3)
f2 = Face(p1, p2, p6)
f3 = Face(p2, p3, p6)
f4 = Face(p1, p3, p6)

def compareEdge(edge1, edge2):
  if (edge1.p1.x == edge2.p1.x and edge1.p1.y == edge2.p1.y and edge1.p1.z == edge2.p1.z \
    and edge1.p2.x == edge2.p2.x and edge1.p2.y == edge2.p2.y and edge1.p2.z == edge2.p2.z) \
  or (edge1.p1.x == edge2.p2.x and edge1.p1.y == edge2.p2.y and edge1.p1.z == edge2.p2.z \
    and edge1.p2.x == edge2.p1.x and edge1.p2.y == edge2.p1.y and edge1.p2.z == edge2.p1.z):
    return True

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

edgetest = [L1, L2, L1_rep, L1_rev, L3]
# for e in edgetest:
# 	print e.__dict__
newlist = getOutRing(edgetest)
for i in newlist:
	print i.__dict__

def unique(edgeList):
    temp = []
    index = []

    for i in range(len(edgeList)):
      for j in range(i+1, len(edgeList)):
        if compareEdge(edgeList[i], edgeList[j]):
        # edgeList[i] == edgeList[j]:
          index.append(i)
          index.append(j)
    
    for i in range(len(edgeList)):
      # print i
      if i not in index:
        # print 'not in index'
        # print i
        temp.append(edgeList[i])

    # print index
    return temp




# print L1.__dict__ 
# print L1.pointliset == L1_rev.pointliset.reverse()

# print L1.pointliset
# print L1_rev.pointliset.reverse()
# list1 = [1,2,2,3, 4,4, 5, 6,6]
# print unique(list1)
# print compareEdge(L1, L1_rev)