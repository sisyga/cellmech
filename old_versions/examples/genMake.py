#!/usr/bin/python  -u

PYTHONUNBUFFERED=1

import sys, argparse 
import os

''' Generate makefile for modLink / cellMech stochastic simulation steps.'''

pngs = []

def writeMakefile(n, file):
    print "CELLMECH="+args.cellmechDir
    print
    print firstLineALL(n)
    for i in range(1, n+1):
        print modLinkLines(i, file, n)
        print cellMechLines(i, n)

def listPNG(n):
	for i in range(1, n+1):
		pngs.append(" \\" + "\n\t" + " png/" + str(i).zfill(digits) + ".png")

def firstLineALL(n):
    return ''.join(['ALL:'] + pngs + [' pickles/cellMech.', str(n).zfill(digits), '.pickle.gz', '\n', '\n'] + ['include ${CELLMECH}/png.make', '\n'])

def modLinkLines(i, file, n):
    first = ''
    if(i == 1):
        first = ''.join(['pickles/modLink.', str(i).zfill(digits), '.pickle.gz : ', file, os.linesep, '\t', 'mkdir -p pickles'])
    else:
        first = ''.join(['pickles/modLink.', str(i).zfill(digits), '.pickle.gz : ', 'pickles/cellMech.', str(i - 1).zfill(digits), '.pickle.gz'])
    second = ''.join(['\t', '${CELLMECH}/modLink.py -i  $< -o $@ -c conf'])
    return ''.join([first, '\n', second, '\n'])

def cellMechLines(i, n):
	first = ''.join(['pickles/cellMech.', str(i).zfill(digits), '.pickle.gz : pickles/modLink.', str(i).zfill(digits), '.pickle.gz'])
	d2 = ''
	if args.d2 : 
		d2 = " -d2 "
	unr = ''
	if args.unrestrict : 
		unr = " -unrestrict "
	last = ''.join(['\t', '${CELLMECH}/cellMech.py -dt 0.001 -nmax 10000', d2, unr, ' -i $< -o $@  > /dev/null'])
	return ''.join([ first, '\n', last, '\n'])
        

#            MAIN 

p = argparse.ArgumentParser(
		description='Generate makefile for modLink / cellMech stochastic simulation steps.')
p.add_argument('-d', action='store', dest='cellmechDir', metavar='DIR', 
		default='../cellmech', help='folder containing the cellmech files')
p.add_argument('-i', action='store', dest='inputFile', metavar='FILE', 
		required=True, help='input pickle FILE')
p.add_argument('-n', action='store', dest='n',  default=10,
		type=int, help='number of steps')
p.add_argument('-d2', action = 'store_true', help = 'Two dimensional configuration.')
p.add_argument('-unrestrict', action = 'store_true', help = 'always update all links.')

args=p.parse_args()

digits = len(str(args.n))
listPNG(args.n)
writeMakefile(args.n, args.inputFile)

