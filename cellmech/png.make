png/%.png : pickles/%.pickle.gz povray.range povray.colorscale
	${CELLMECH}/pickleVis_pov.sh $< povray.range povray.colorscale  > $@

%.png : %.pickle.gz povray.range povray.colorscale
	${CELLMECH}/pickleVis_pov.sh $< povray.range povray.colorscale  > $@

png/%.png : pickles/cellMech.%.pickle.gz povray.range povray.colorscale
	mkdir -p png
	${CELLMECH}/pickleVis_pov.sh $< povray.range povray.colorscale  > $@
	ln -sf $@ 	newest.png

png/start.png : start.pickle.gz povray.range povray.colorscale
	mkdir -p png
	${CELLMECH}/pickleVis_pov.sh $< povray.range povray.colorscale  > $@
	ln -sf $@ 	newest.png

povray.colorscale :
	echo 2.0 > $@

avi : png
	avconv -i png/%03d.png -vcodec libx264 $$(basename $$PWD).avi
