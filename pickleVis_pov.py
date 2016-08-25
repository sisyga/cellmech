#!/usr/bin/python  -u

PYTHONUNBUFFERED=1

import os, sys, argparse 

sys.path.append('..')  
import cell
from cell import null, ex, ey, ez

def povBase():
	base =  "global_settings { ambient_light rgb <2.000000, 2.000000, 2.000000> }\n\nbackground { color rgb <0.000000, 0.000000, 0.000000> }\n\nlight_source { <218.217890, 436.435780, 872.871561>\n\tcolor rgb <0.800000, 0.800000, 0.800000>\n}\nlight_source { <-872.871561, -218.217890, -436.435780>\n\tcolor rgb <0.300000, 0.300000, 0.300000>\n}\n\ncamera {\n\tright <-1.0, 0, 0>      //visual uses right-handed coord. system"
	if args.view == 'v':
		base = base + "\n\tlocation <0.000000, 0.000000, " + str(args.range* 1.73) +">\n"
		base = base + "\tsky <0.000000, 1.000000, 0.000000>\n\tlook_at <0.000000, 0.000000, 0.000000>\n\tangle 60.000000\n\trotate <0, 0, 0>\n}\n"
	if args.view == 'd':
		base = base + "\n\tlocation <0.000000, 0.000000, " + str(-args.range* 1.73) +">\n"
		base = base + "\tsky <0.000000, 1.000000, 0.000000>\n\tlook_at <0.000000, 0.000000, 0.000000>\n\tangle 60.000000\n\trotate <0, 0, 0>\n}\n"
	if args.view == 'a':
		base = base + "\n\tlocation <0.000000, " + str(args.range* 1.73) +", 0.00000 >\n"
		base = base + "\tsky <0.000000, 0.000000, 1.000000>\n\tlook_at <0.000000, 0.000000, 0.000000>\n\tangle 60.000000\n\trotate <0, 0, 0>\n}\n"
	if args.view == 'p':
		base = base + "\n\tlocation <0.000000, " + str(-args.range* 1.73) +", 0.00000 >\n"
		base = base + "\tsky <0.000000, 0.000000, 1.000000>\n\tlook_at <0.000000, 0.000000, 0.000000>\n\tangle 60.000000\n\trotate <0, 0, 0>\n}\n"
	if args.view == 'l':
		base = base + "\n\tlocation <" + str(args.range* 1.73) +", 0.0000, 0.00000 >\n"
		base = base + "\tsky <0.000000, 1.000000, 0.000000>\n\tlook_at <0.000000, 0.000000, 0.000000>\n\tangle 60.000000\n\trotate <0, 0, 0>\n}\n"
	if args.view == 'r':
		base = base + "\n\tlocation <" + str(-args.range* 1.73) +", 0.00000, 0.0000 >\n"
		base = base + "\tsky <0.000000, 1.000000, 0.000000>\n\tlook_at <0.000000, 0.000000, 0.000000>\n\tangle 60.000000\n\trotate <0, 0, 0>\n}\n"
	return base

def genSphere(n, radius, color):
	sphere = "sphere {"
	R=n.getR()
	location = "\t<" + str(R[0]) + ", " + str(R[1]) + ", " + str(R[2]) +  ">, " + str(radius)
	color = "\t" + "pigment {color rgbt <" + color + ">}"
	return sphere + "\n" + location + "\n" + color + "\n}\n" 

def genLink(l, radius, color) :
	cylinder = "cylinder {"
	R=l.n1.getR()
	n1pos = "\t<" + str(R[0]) + ", " + str(R[1]) + ", " + str(R[2]) + ">,"
	R=l.n2.getR()
	n2pos = "<" + str(R[0]) + ", " + str(R[1]) + ", " + str(R[2]) + ">,"
	radius = str(radius)  + "\n\t"
	color = "pigment { color rgbt <" + color + ", 0.0> }\n}\n"
	return cylinder + "\n" + n1pos + n2pos + radius + color

def povrayGraph():
	red    = "1.000000, 0.000000, 0.000000"
	blue   = "0.000000, 0.000000, 1.000000"
	green  = "0.000000, 0.500000, 0.000000"
	yellow = "1.000000, 1.000000, 0.000000"
	
	base = povBase()
	color = ""
	print base
	for n in c.nodes: 
		if   n.state < 0.5 : color = red
		elif n.state > 1.5 : color = yellow
		#if n.prolif == 0: color=red
		#elif n.prolif > 0.5: color = yellow
		#elif n.prolif < -0.5 : color = green
		genSphere(n, radius = 0.25, color = color)
		print genSphere(n, radius=0.25, color = color)
	for l in c.links: 
		if l.d < l.d0 : 
			if args.link_rule == 1:
				q=(l.d0-l.d)/l.d0*10
			elif args.link_rule == 2:
				q=l.d0-l.d
		else :
			if args.link_rule == 1:
				q=(l.d-l.d0)/l.d0*10
			elif args.link_rule == 2:
				q=l.d-l.d0
			color= "0, " + str(q) + ", " + str(1-q)
		link = genLink(l = l, radius = 0.05, color = color)
		print link

#            MAIN 

p = argparse.ArgumentParser(description='Visualization of multicellular simulation pickle files.')
p.add_argument('-i', dest='inputFile', metavar='FILE', required=True,
		type=cell.checkFile, help='pickled configuration data in FILE')
p.add_argument('-o', dest='outputFile', metavar='FILE',
                help='capture output window into FILE.xwd')
p.add_argument('-r', dest='range', type=float, help='display range')
p.add_argument('-l', dest='link_rule', type=int, default=1,
		 help='link color rule: 1, 2, ... ')
p.add_argument('-w', dest='view', default='v', help='v: ventral, d:dorsal, a:anterior, p:posterior, l:left, r:right')
args=p.parse_args()

c=cell.loadData(args.inputFile)

povrayGraph()


