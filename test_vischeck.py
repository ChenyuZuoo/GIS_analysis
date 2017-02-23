# -*- coding: utf-8 -*- 
from __future__ import division

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

L1 = Edge(p1, p2)
L2 = Edge(p2, p3)
# L1_rev = Edge(p2, p1)
# L1_rep = Edge(p1, p2)
L3 = Edge(p3, p4)

f1 = Face(p1, p2, p3)
f2 = Face(p1, p2, p6)
f3 = Face(p2, p3, p6)

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
  print vis
  
  # if vis > 0.5:   
  ### 0.5 as threshold is too bigh here, so change it to strict 0
  ### might be time comsuming
  if vis > 0:
    return True
  else:
    return False

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

p = getMiddlePnt(p1, p2)
p = getInnerPnt(f1, p6)
# print p.__dict__
print checkVisibility(f1, p)
print checkVisibility(f2, p)
print checkVisibility(f3, p)
print checkVisibility(f2, p9)