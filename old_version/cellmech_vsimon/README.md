# cellmech
<pre>
requires pip, pyhull, numpy
	sudo apt-get install python-pip python-dev
	sudo pip install pyhull
	sudo pip install numpy

	a suitable makefile alternates two programs:
		cellmech/cellmech.py:	calculates mechanical equilibrium
		cellmech/modlink.py:    link addition/removal
	
run the example simulation, see examples/readme.txt

extracting info from pickle files:
	examples/pickleRead.py
		the possible attributes of nodes and links are defined in cellmech/cell.py
		some examples:
			getR()  position
			F		net force 
			F0		external force
			M		net torque
			M0		external torque
	
pickle file visualization:
	cellmech/png.make	calling
		pickleVis_pov.sh 	calling
			pickleVis_pov.py
			povray

	the area shown is set by the content of the povray.range file
		(see examples/povray.range)

	the image resolution is hardcoded to 480x480 in pickleVis_pov.sh

modLink.py:
            reads parameters

            cell.loadData()
            delAttrDist()   removes dist attribute from nodes
            checkLinkX()    removes links intersecting with other links
                intersect()
                    ccw()
            delLinkList()   appends to_del list
            addLinkList()   appends to_add list
            pickEvent()     either removes or adds object
            default_update_d0() or other externally defined function
            cell.dumpData()

</pre>
