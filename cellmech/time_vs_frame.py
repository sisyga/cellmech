#!/usr/bin/python  -u
PYTHONUNBUFFERED=1

import os 
import sys 
import argparse 
	

def searchLog(file):
	time={}
	dt = 0
	frame = 0
	with open(file, "r") as f:
		for line in f:
			w=line.split()
			dt=w[-1]
			q=w[3].split('.pickle')
			number=q[0].split('.')[-1]
			#print number,dt
			time[int(number)] = dt
	return time

	
#           MAIN 

p = argparse.ArgumentParser(description='Change step-number to time.')
p.add_argument('-l', action='store', dest='logFile', metavar='FILE', 
		help='Log FILE')
p.add_argument('-v', action='store_true', help='verbose (debug)')
args=p.parse_args()

time = searchLog(args.logFile)
t=0;
for i in sorted(time.keys()):
	t=t+float(time[i])
	print i, t

