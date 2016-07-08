#!/usr/bin/python  -u

import sys 
import argparse 

sys.path.append('..')
import cellmech

leftSelected  = [ ]
rightSelected = [ ]	
percent = 0.1	

# selecteding the points to impact with streching force
def stretch():
	S=[n.getR()[0] for n in c.nodes]
	L=max(S)-min(S);
	print "L: ",L
	side = L / 2.		# coordinates bw += side 
	for n in c.nodes:
		if n.getR()[0] < - side + side * percent:
			leftSelected.append(n)
		if n.getR()[0] > side - side * percent:
			rightSelected.append(n)
										#Fleft*len(left) = Fright*len(right)
	Fleft  = args.mltpl * len(rightSelected)/float(len(leftSelected))	
	Fright = args.mltpl 

	for n in leftSelected:
		n.F0=-Fleft*cellmech.ex;

	for n in rightSelected:
		n.F0=Fright*cellmech.ex;

	for n in leftSelected + rightSelected:
		n.kphi=30.0
		n.bend=30.0
		for l in n.links:
			l.limit = 50.0
		
p = argparse.ArgumentParser(
		description='Adding stretching forces to a configuration.')
p.add_argument('-i', action='store', dest='inputFile', metavar='FILE', 
		required=True, help='input pickled configuration data')
p.add_argument('-o', action='store', dest='outputFile', metavar='FILE', 
		required=True, help='output pickled configuration data')
p.add_argument('-m', action='store', dest='mltpl',  default=1.0,
		type=float, help='Stretching force scaling factor')

args=p.parse_args()

c=cellmech.loadData(args.inputFile)
stretch()
cellmech.dumpData(c,args.outputFile)
