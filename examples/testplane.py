#!/usr/bin/python  -u

import numpy as np
from numpy import array
import argparse
import sys

sys.path.append('..')
import cellmech

d0min=0.8       # min distance between cells
d0max=2.        # max distance connected by links
#######################################################


def generatePoint(L):
	X0=(np.random.rand()-.5)*L 
	Y0=(np.random.rand()-.5)*L 
	Z0=0.
	return array([X0,Y0,Z0])

	

#            MAIN 

p = argparse.ArgumentParser(
		description='Simple plane configuration setup cellmech.')
p.add_argument('-o', action='store', dest='pickleFile', metavar='FILE', 
		default='config.pickle', help='output pickle FILE')
p.add_argument('-b', action='store', dest='bend',  default=10.0,
		type=float, help='bend')
p.add_argument('-t', action='store', dest='twist',  default=1.0,
		type=float, help='twist')
p.add_argument('-L', default=10.0, type=float, help='linear system size')

args=p.parse_args()

N=int(args.L*args.L)           #1 cell / unit area
c=cellmech.Configuration(N)

for i in range(N):
	while True:
		R1=generatePoint(args.L)
		OK=True
		for n in c.nodes:
			d=cellmech.norm(n.getR()-R1,"len")
			if d<d0min : 
				OK=False
				break
		if OK : break
	n=cellmech.node(c,R1)
	n.twist=args.twist
	n.bend=args.bend
	print i

for i,j in cellmech.VoronoiNeighbors( [ n.getR() for n in c.nodes ], d0max ):
	c.nodes[i].addLinkTo(c.nodes[j])

print "---"
cellmech.dumpData(c,args.pickleFile)
