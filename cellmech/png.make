png/%.png : pickles/%.pickle.gz povray.range
	${CELLMECH}/pickleVis_pov.sh $< povray.range  > $@

%.png : %.pickle.gz povray.range
	${CELLMECH}/pickleVis_pov.sh $< povray.range  > $@

png/%.png : pickles/cellMech.%.pickle.gz povray.range
	mkdir -p png
	${CELLMECH}/pickleVis_pov.sh $< povray.range  > $@
	ln -sf $@ 	newest.png

png/start.png : start.pickle.gz povray.range
	mkdir -p png
	${CELLMECH}/pickleVis_pov.sh $< povray.range  > $@
	ln -sf $@ 	newest.png
