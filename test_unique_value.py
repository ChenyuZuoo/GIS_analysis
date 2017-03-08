# -*- coding: utf-8 -*- 

a = ['a', 'b', 'c', 'd', 'a', 'a']
b = [x for x in a if a.count(x) == 1]
# print b

# for x in a:
#   print a.count(x) == 1

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

p1 = Point(2,0,0)
p2 = Point(0.5,0.5,0.5)
p3 = Point(-1, 0,0)
p4 = Point(-2,0,0)
p5 = Point(2, 2, 0.5)
p6 = Point(2, 0.5, 0.5)

L1 = Edge(p1, p2)
L2 = Edge(p2, p3)
L1_rev = Edge(p2, p1)
L1_rep = Edge(p1, p2)
L3 = Edge(p3, p4)
# print L1, L2, L1_rev, L1_rep
# print L1.__dict__ == L1_rep.__dict__
edgeSet = [L1, L2, L1_rep, L1_rep, L3]

def filter_by_uniqueness(edgeList):
    temp = []
    index = []
    for i in range(len(edgeList)):
      for j in range(i+1, len(edgeList)):
        # print i
        # print edgeList[i]
        # print j
        # print edgeList[j]
        if edgeList[i].__dict__ == edgeList[j].__dict__:
          index.append(i)
    # print index
    
    for i in range(len(edgeList)):
      if i not in index:
        temp.append(edgeList[i])

    return temp

edgenew =  filter_by_uniqueness(edgeSet)
for i in edgenew:
  print i.p1.__dict__, i.p2.__dict__
