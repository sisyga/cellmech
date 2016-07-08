png/%.png : pickles/cellMech.%.pickle.gz povray.range
	mkdir -p png
	../cellmech/pickleVis_pov.sh $< povray.range  > $@
	ln -sf $@ 	newest.png

png/start.png : start.pickle.gz povray.range
	mkdir -p png
	../cellmech/pickleVis_pov.sh $< povray.range  > $@
	ln -sf $@ 	newest.png
