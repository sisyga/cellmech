%.avi : $(wildcard png/*.png)
	avconv -y -i png/%04d.png -vcodec libx264 $$(basename $$PWD).avi
