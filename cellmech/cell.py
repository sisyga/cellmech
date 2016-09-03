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

def norm(v,what=None):
	n=np.sqrt(np.dot(v,v))
	if n>1e-6 : 
		u=v/n
	else:
		u,n=null,0.
	if what=='vec':
		return u
	if what=='mag' or what=='len':
		return n
	return [u,n]

def R(dphi,v):
	axis,theta=norm(dphi)
	if (theta<1e-6): return v
	a = np.cos(theta/2)
	b,c,d = axis*np.sin(theta/2)
	return np.dot( np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
							 [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
							 [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]]),v)

class Configuration():
	def __init__(self,num=None):
		self.nodes=[]
		self.links=[]
		self.nNodes=0	#bookkeeping for state vector position 
		if num is not None:
			self.N=num
			self.x=np.zeros(6*num)
		else:
			self.N=0
			self.x=None

def dumpData(configuration,fname):
	sys.setrecursionlimit(50000)
	with gzip.GzipFile(fname,'wb') as out:
		pickle.dump(configuration,out)

def loadData(fname):
	with gzip.GzipFile(fname,'rb') as inf:
		config=pickle.load(inf)

	print >> sys.stderr, "N=", config.N
	print >> sys.stderr, len(config.x), "dof."
	print >> sys.stderr, len(config.nodes), "nodes."
	print >> sys.stderr, len(config.links), "links."
	return config
	

class node(object):
	def __init__(self,config,r,r0=None,state=-1):

		self.config=config	#reference to a specific configuration
													#state vector indices
		if config.nNodes == len(config.x):
			config.N+=1
			config.x.append(np.zeros(6))

		self.r=[config.nNodes, config.nNodes+1, config.nNodes+2]; 	
		config.nNodes+=3	
		config.x[self.r]=r
		if r0 is None: r0=r
		self.r0=r0
													#state vector indices
		self.phi=[config.nNodes, config.nNodes+1, config.nNodes+2]; 
		config.nNodes+=3	
		config.x[self.phi]=null
		self.phi0=null

		self.links=[]
		self.F=null			#net force
		self.F0=null		#external force
		self.M=null			#net torque

		self.bend=10.		#bending rigidity
		self.twist=1.		#torsion spring constant
		self.kr=0.			#spring to anchor point
		self.kphi=0.		#torsion spring to anchor point

		self.state=state	#discrete state variable, like 0, 1, 2, ..
		self.normal = None	#normal vector, if surface node
		config.nodes.append(self)

	def findLinkTo(s,n):
		for l in s.links:
			if ((l.n1==s) & (l.n2==n)) | ((l.n1==n) & (l.n2==s)):
				return l
		return None

	def addLinkTo(s,n,t1=None,t2=None,d0=None,k=None):
		if s.findLinkTo(n) == None:
			l=link(s.config,s,n,t1,t2,d0,k)
			s.links.append(l)
			n.links.append(l)
		return l

	def rmLinkTo(s,n):
		l=s.findLinkTo(n)
		if l != None:
			s.links.remove(l)
			n.links.remove(l)
			s.config.links.remove(l)

	def getR(s):
		return s.config.x[s.r]
		
class link(object):
	def __init__(self,config,n1,n2,t1=None,t2=None,d0=None,k=None,n=None,norm1=None, norm2=None):

		self.config=config				#reference to a specific configuration
		self.n1,self.n2=n1,n2			#nodes at endpoints
										#actual direction (e) and length (d)
		self.e, self.d = norm(config.x[n2.r]-config.x[n1.r])	

		if k is None: k=15.			
		self.k=k;						#spring parameter

		if d0 is None: d0=self.d
		self.d0=d0;						#equilibrium distance
										#preferred directions 
		if t1 is None : self.t1=R(-config.x[self.n1.phi],self.e) 
		else : 			self.t1=norm(t1,"vec")
		if t2 is None : self.t2=R(-config.x[self.n2.phi],-self.e)
		else : 			self.t2=norm(t2,"vec")
		if n is None : 
			n,q=norm(np.cross(self.e, ez)) 			# n is perpendicular to e
							# n is perpendicular to z (l is in the x-y plane)
			if (q<1e-5) : n=norm(np.cross(self.e,ex),"vec")	# e || ez   =>	n is perpendicular to x
		self.n=n
		if norm1 is None : norm1=n
		self.norm1=norm1
		if norm2 is None : norm2=n
		self.norm2=norm2
		self.M1, self.M2, self.F1, self.F2 = null, null, null, null
		config.links.append(self)

	def remove(s):
		s.n1.links.remove(s)
		s.n2.links.remove(s)
		s.config.links.remove(s)
	


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
