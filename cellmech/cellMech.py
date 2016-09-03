#!/usr/bin/python  -u

import sys, argparse, os 
import numpy as np
import copy
from sets import Set

#
# see http://stackoverflow.com/questions/1046628/importing-python-modules-from-different-working-directory
#
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import cell 


#torques, forces acting on the LINK
#here we change properties of link objects, parts of cellmech.Configuration object c
def updateLinkForces(l,x) :
	l.e,l.d=cell.norm(x[l.n2.r]-x[l.n1.r])		#actual link direction and length
	l.n=cell.norm(l.norm1+l.norm2,"vec")								# Eq. 7
	norm2=np.cross(l.e,l.n)											# Eq. 9
	if args.d2:
		l.M1=l.n1.bend*np.cross( l.e,cell.R(x[l.n1.phi],l.t1))
		l.M2=l.n2.bend*np.cross(-l.e,cell.R(x[l.n2.phi],l.t2))
	else :
		l.M1=l.n1.bend*np.cross( l.e,cell.R(x[l.n1.phi],l.t1)) +  \
			l.n1.twist*np.cross( l.n,cell.R(x[l.n1.phi],l.norm1))	# Eqs 3, 18, 19
		l.M2=l.n2.bend*np.cross(-l.e,cell.R(x[l.n2.phi],l.t2)) +  \
			l.n2.twist*np.cross( l.n,cell.R(x[l.n2.phi],l.norm2))	# Eqs 5, 18, 19	
	M=l.M1+l.M2
	l.F2=np.dot(M,l.n)/l.d*norm2-np.dot(M,norm2)/l.d*l.n+l.k*(l.d-l.d0)*l.e;  #Eqs (11), (12) and (13) for F2=-F1!
  	l.F1=-l.F2;

def F(x) :
	ret=np.zeros(len(x))	#change in state vector, 6 dof for each unit
	for l in links : updateLinkForces(l,x)
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
	global c
	x = copy.deepcopy(c.x)
	h=args.dt
	print "h=", h
	print "mechEquilibrium START"
	for i in range(args.nMax):
		k1=F(x)
		Q=np.dot(k1,k1)/len(neighbors) 
		print i,Q, max([ x[n.r][0] for n in c.nodes ])
		if i%100 == 0 and args.dumpdir is not None : 
			cell.dumpData("%s/%05d.pickle" % (args.dumpdir,i))
		if Q < args.qMin: break
		k1*=h
		k2=h*F(x+k1/2)
		k3=h*F(x+k2/2)
		k4=h*F(x+k3)
		x+=(k1+2*k2+2*k3+k4)/6.	
	print "mechEquilibrium END"
	c.x = x

def findChangedNodes():
	n12 = []
	for n in c.nodes:
		if hasattr(n, 'dist'):
			n12.append(n)
	return n12

def setDist(n0, dist):
	global links
	global neighbors

	if dist > distLimit: return
	if hasattr(n0, 'dist') and n0.dist < dist: return
	n0.dist = dist
	neighbors.add(n0)
	for l in n0.links:
		if l.n1 == n0:
			setDist(l.n2, dist+1)			
		else: 
			setDist(l.n1, dist+1)
		links.add(l)			
		
#            MAIN 

p = argparse.ArgumentParser(
		description='Moves configuration to mechanical equilibrium.')
p.add_argument('-i', action='store', dest='inputFile', metavar='FILE', 
		required=True, help='initial configuration')
p.add_argument('-o', action='store', dest='outputFile', metavar='FILE',
                required=True, help='final, equilibrium configuration ')
p.add_argument('-dumpdir', action='store', metavar='DIR', default=None,
		help='intermediate configurations in DIR')
p.add_argument('-dt', action='store', dest='dt',  default=0.01,
		type=float, help='timestep for relaxation ')
p.add_argument('-nmax', action='store', dest='nMax',  default=3000,
		type=int,help='max number of relaxation steps')
p.add_argument('-qmin', action='store', dest='qMin',  default=.001,
		type=float, help='target ||F|| value to stop relaxation ')
p.add_argument('-d2', action = 'store_true', help = '2D configuration.')
p.add_argument('-unrestrict', action = 'store_true', help = 'Boolean: Do not restrict solver to nodes near changed links.')
args=p.parse_args()


c=cell.loadData(args.inputFile)

# check version of pickle file ( dist attribute is set for changed nodes )
n12 = findChangedNodes()
distLimit = 3
neighbors  = set()
links      = set()
if len(n12) == 2 and not args.unrestrict:
	setDist(n12[0], 0)
	setDist(n12[1], 0)
else:
	neighbors = c.nodes
	links = c.links

mechEquilibrium()

cell.dumpData(c,args.outputFile)
