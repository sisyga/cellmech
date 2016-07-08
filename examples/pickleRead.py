#!/usr/bin/python  -u


import os, sys, argparse 

sys.path.append('..')  
import cellmech

p = argparse.ArgumentParser(description='Accessing simulation pickle files.')
p.add_argument('-i', dest='inputFile', metavar='FILE', required=True,
		type=cellmech.checkFile, help='pickled configuration data in FILE')
args=p.parse_args()

c=cellmech.loadData(args.inputFile)

for n in c.nodes:
	print n.F0

