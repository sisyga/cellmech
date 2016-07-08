#!/usr/bin/python  -u

PYTHONUNBUFFERED=1

import sys, argparse 
import numpy as np
import cell
import copy
from sets import Set
import datetime

depthLimit = 3
neighbors  = set()
links      = set()
isNewConf = True

#torques, forces acting on the LINK
def updateMech(l,x) :
	l.e,l.d=cell.norm(x[l.n2.r]-x[l.n1.r])			#update actual link direction and length
	l.n=cell.norm(l.norm1+l.norm2)[0]									# Eq. 7
	norm2=np.cross(l.e,l.n)											# Eq. 9
	if args.d2:
		l.M1=l.n1.bend*np.cross(l.e,cell.R(x[l.n1.phi],l.t1))
		l.M2=l.n2.bend*np.cross(-l.e,cell.R(x[l.n2.phi],l.t2))
	else :
		l.M1=l.n1.bend*np.cross(l.e,cell.R(x[l.n1.phi],l.t1)) +  \
			l.n1.twist*np.cross(l.n,cell.R(x[l.n1.phi],l.norm1))	# Eqs 3, 18, 19
		l.M2=l.n2.bend*np.cross(-l.e,cell.R(x[l.n2.phi],l.t2)) +  \
			l.n2.twist*np.cross(l.n,cell.R(x[l.n2.phi],l.norm2))	# Eqs 5, 18, 19	
	M=l.M1+l.M2
	l.F2=np.dot(M,l.n)/l.d*norm2-np.dot(M,norm2)/l.d*l.n+l.k*(l.d-l.d0)*l.e;  #Eqs (11), (12) and (13) for F2=-F1!
  	l.F1=-l.F2;

def F(x) :
	ret=np.zeros(6*cell.N)	#change in state vector, 6 dof for each unit
	for l in links : updateMech(l,x)
	for n in neighbors :
		n.F=-n.kr*(x[n.r]-n.r0) + n.F0
		n.M=-n.kphi*(x[n.phi]-n.phi0)
	for l in links : 
			l.n1.F=l.n1.F-l.F1
			l.n2.F=l.n2.F-l.F2
			l.n1.M=l.n1.M-l.M1
			l.n2.M=l.n2.M-l.M2
	for n in neighbors :
		ret[n.r]=n.F
		ret[n.phi]=n.M
	return ret

def mechEquilibrium() :
	x = copy.deepcopy(cell.x)
	h=args.dt
	print "h=", h
	i=0
	print "mechEquilibrium START"
	while True :
		k1=F(x)
		Q=np.dot(k1,k1)/cell.N 
		print i,Q
		if i%100 == 0 and args.dump is not None : 
			cell.dumpData("%s/%05d.pickle" % (args.dump,i))
		if ( (i>args.nMax) or (Q < args.qMin) ) : break
		k1*=h
		k2=h*F(x+k1/2)
		k3=h*F(x+k2/2)
		k4=h*F(x+k3)
		x+=(k1+2*k2+2*k3+k4)/6.	
		i+=1
	print "mechEquilibrium END"
	cell.x = copy.deepcopy(x)

def getPrevModification():
	with open(args.f) as f:
		content = f.readlines()[0].split("\t")
		id1 = int(content[0])
		id2 = int(content[1])
		return id1, id2

def getLinks():
	n12 = []
	for n in cell.nodes:
		if hasattr(n, 'depth'):
			n12.append(n)
	if len(n12) != 2:
		isNewConf=False
	return n12

def setNeighbors(n0, depth):
	global links
	global neighbors
	if depth > depthLimit: return
	if hasattr(n0, 'depth') and n0.depth < depth: return
	n0.depth = depth
	neighbors.add(n0)
	for (l,d) in n0.links:
		if d == 1:
			setNeighbors(l.n2, depth+1)			
		else: 
			setNeighbors(l.n1, depth+1)
		links.add(l)			
		
#            MAIN 

p = argparse.ArgumentParser(description='Multicellular simulation.')
p.add_argument('-i', action='store', dest='inputFile', metavar='FILE', 
		required=True, type=cell.checkFile, help='initial configuration data in FILE')
p.add_argument('-o', action='store', dest='outputFile', metavar='FILE',
                required=True, help='final configuration data in FILE')
p.add_argument('-dump', action='store', metavar='DIR', default=None,
		type=cell.checkDir, help='intermediate configurations in DIR')
p.add_argument('-dt', action='store', dest='dt',  default=0.01,
		type=float, help='timestep during mech relaxation ')
p.add_argument('-nmax', action='store', dest='nMax',  default=3000,
		type=int,help='max number of timesteps during mech relaxation ')
p.add_argument('-qmin', action='store', dest='qMin',  default=.001,
		type=float, help='||F|| value to stop mech relaxation ')
p.add_argument('-d2', action = 'store_true', help = 'Two dimension configuration.')
p.add_argument('-old', action = 'store_true', help = 'Boolean: Use old configuration.')
args=p.parse_args()


cell.loadData(args.inputFile)
n12 = getLinks()

# check version of pickle file ( made with new modLink - modified depth of changed link`s nodules )
if isNewConf == True and len(n12) == 2 and not args.old:
	print "New modLink configuration."
	setNeighbors(n12[0], 0)
	setNeighbors(n12[1], 0)
	mechEquilibrium()
else:
	print "Old modLink configuration."
	neighbors = cell.nodes
	links = cell.links
	mechEquilibrium()

cell.dumpData(args.outputFile)
