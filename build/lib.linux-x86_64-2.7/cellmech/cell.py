#!/usr/bin/python  -u

PYTHONUNBUFFERED=1

import os 
import sys 
import numpy as np
import pickle
import gzip
import argparse
from numpy import array


null=array([0.0,0.0,0.0])
ex=array([1.0,0.0,0.0])
ey=array([0.0,1.0,0.0])
ez=array([0.0,0.0,1.0])

def norm(v):
	n=np.sqrt(np.dot(v,v))
	if n>1e-6 : return [v/n,n]
	return [null,0.]

def R(dphi,v):
    axis,theta=norm(dphi)
    if (theta<1e-6): return v
    a = np.cos(theta/2)
    b,c,d = axis*np.sin(theta/2)
    return np.dot( np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                             [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                             [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]]),v)


nodes = []
neighbors = []
x=array([])
N=0

nNodes=0		#bookkeeping for state vector position 

class node(object):
	def __init__(self,r,r0=None,state=-1):
		global x, nNodes
		self.r=[nNodes,nNodes+1,nNodes+2]; nNodes+=3	#index array
		x[self.r]=r
		if r0 is None: r0=r
		self.r0=r0
		self.phi=[nNodes,nNodes+1,nNodes+2]; nNodes+=3	#index array
		x[self.phi]=null
		self.phi0=null
		self.links=[]
		self.F=null
		self.F0=null
		self.M=null
		self.bend=10.
		self.twist=1.
		self.kr=0.
		self.kphi=0.
		self.state=state
		self.normal = None
		nodes.append(self)

	def findLinkTo(s,n):
		for (l,d) in s.links:
			if ((l.n1==s) & (l.n2==n)) | ((l.n1==n) & (l.n2==s)):
				return [True, l,d]
		return [False, None,None]

	def addLinkTo(s,n,t1=None,t2=None,d0=None,k=None):
		found,l,d=s.findLinkTo(n)
		if not found:
			l=link(s,n,t1,t2,d0,k)
			s.links.append((l,1))
			n.links.append((l,-1))
		return l

	def rmLinkTo(s,n):
		found,l,d=s.findLinkTo(n)
		if found:
			s.links.remove((l,d))
			n.links.remove((l,-d))
			links.remove(l)

	def getR(s):
		#print x[s.r]
		return x[s.r]
		
links = [];

class link(object):
	def __init__(self,n1,n2,t1=None,t2=None,d0=None,k=None,n=None,norm1=None, norm2=None):
		self.n1,self.n2=n1,n2					#nodes at endpoints
		self.e,self.d=norm(x[n2.r]-x[n1.r])			#actual direction (e) and length (d)
		if k is None: k=15.			
		self.k=k;						#spring parameter
		if d0 is None: d0=self.d
		self.d0=d0;						#equilibrium distance
		if t1 is None : self.t1=R(-x[self.n1.phi],self.e) 	#preferred directions 
		else : self.t1=norm(t1)[0]
		if t2 is None : self.t2=R(-x[self.n2.phi],-self.e)
		else : self.t2=norm(t2)[0]
		if n is None : 
			n,q=norm(np.cross(self.e, ez)) 			# n is perpendicular to e
									# n is perpendicular to z (is in the x-y plane)
			if (q<1e-5) : n=norm(np.cross(self.e,ex))[0]	# e || ez   =>    n is perpendicular to x
		self.n=n
		if norm1 is None : norm1=n
		self.norm1=norm1
		if norm2 is None : norm2=n
		self.norm2=norm2
		self.M1, self.M2, self.F1, self.F2 = null, null, null, null
		links.append(self)

	def remove(s):
		s.n1.links.remove((s,1))
		s.n2.links.remove((s,-1))
		links.remove(s)

def dumpData(fname):
	sys.setrecursionlimit(50000)
	with gzip.GzipFile(fname,'wb') as out:
        	pickle.dump(N,out)
        	pickle.dump(x,out)
        	pickle.dump((nodes,links),out)

def loadData(fname):
        global N,x,nodes,links
	try:
        	with gzip.GzipFile(fname,'rb') as inf:
        		N=pickle.load(inf)
        		x=pickle.load(inf)
        		nodes,links=pickle.load(inf)
	except:
        	with open(fname,'rb') as inf:
        		N=pickle.load(inf)
        		x=pickle.load(inf)
        		nodes,links=pickle.load(inf)
       	print >> sys.stderr, "N=", N
       	print >> sys.stderr, len(x), "dof."
        print >> sys.stderr, len(nodes), "nodes."
        print >> sys.stderr, len(links), "links."

def checkDir(s):
        if os.path.exists(s):
                if not os.path.isdir(s):
                        raise argparse.ArgumentTypeError("%s is not a dir." % s)
        else: os.makedirs(s)
        return s

def checkFile(s):
        if os.path.exists(s):
		if s.strip() == '/dev/stdin': 
			return s
                if not os.path.isfile(s):
                        raise argparse.ArgumentTypeError("%s is not a file." % s)
                return s
        else: argparse.ArgumentTypeError("%s does not exists." % s)

