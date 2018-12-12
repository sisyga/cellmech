CELLMECH=../cellmech

ALL: \
	 png/start.png \
	 png/01.png pickles/cellMech.01.pickle.gz

include ${CELLMECH}/png.make

start.pickle.gz : stretched-plane.pickle.gz
	ln -sf $< $@

pickles/cellMech.01.pickle.gz :  stretched-plane.pickle.gz
	${CELLMECH}/cellMech.py -dt 0.001 -d2 -unrestrict -nmax 10000 -i $< -o $@  


