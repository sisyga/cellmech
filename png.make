png/%.png : pickles/cellMech.%.pickle.gz
	mkdir -p png
	../cellmech/pickleVis_pov.sh $< povray.range v > $@
	ln -sf $@ 	newest.png

png/start.png : pipe.pickle.gz
	mkdir -p png
	../cellmech/pickleVis_pov.sh $< povray.range v > $@
	ln -sf $@ 	newest.png
