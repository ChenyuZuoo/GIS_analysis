# -*- coding: utf-8 -*- 
def filter_by_uniqueness(l):
    d = {}
    for elem in l:
        if elem in d:
            d[elem] += 1
        else:
            d[elem] = 1
    uniquel = []
    print d.items()
    for k,v in d.items():
        print k, v
        if v == 1:
            uniquel.append(k)
    return uniquel  
l = ['a', 'b','b','b','b','c']

filter_by_uniqueness(l)