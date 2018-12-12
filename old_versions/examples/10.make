CELLMECH=../cellmech

ALL: \
	 png/01.png \
	 png/02.png \
	 png/03.png \
	 png/04.png \
	 png/05.png \
	 png/06.png \
	 png/07.png \
	 png/08.png \
	 png/09.png \
	 png/10.png pickles/cellMech.10.pickle.gz

include ${CELLMECH}/png.make

pickles/modLink.01.pickle.gz : stretched-plane.pickle.gz
	mkdir -p pickles
	${CELLMECH}/modLink.py -i  $< -o $@ -c conf

pickles/cellMech.01.pickle.gz : pickles/modLink.01.pickle.gz
	${CELLMECH}/cellMech.py -dt 0.001 -nmax 10000 -i $< -o $@  > /dev/null

pickles/modLink.02.pickle.gz : pickles/cellMech.01.pickle.gz
	${CELLMECH}/modLink.py -i  $< -o $@ -c conf

pickles/cellMech.02.pickle.gz : pickles/modLink.02.pickle.gz
	${CELLMECH}/cellMech.py -dt 0.001 -nmax 10000 -i $< -o $@  > /dev/null

pickles/modLink.03.pickle.gz : pickles/cellMech.02.pickle.gz
	${CELLMECH}/modLink.py -i  $< -o $@ -c conf

pickles/cellMech.03.pickle.gz : pickles/modLink.03.pickle.gz
	${CELLMECH}/cellMech.py -dt 0.001 -nmax 10000 -i $< -o $@  > /dev/null

pickles/modLink.04.pickle.gz : pickles/cellMech.03.pickle.gz
	${CELLMECH}/modLink.py -i  $< -o $@ -c conf

pickles/cellMech.04.pickle.gz : pickles/modLink.04.pickle.gz
	${CELLMECH}/cellMech.py -dt 0.001 -nmax 10000 -i $< -o $@  > /dev/null

pickles/modLink.05.pickle.gz : pickles/cellMech.04.pickle.gz
	${CELLMECH}/modLink.py -i  $< -o $@ -c conf

pickles/cellMech.05.pickle.gz : pickles/modLink.05.pickle.gz
	${CELLMECH}/cellMech.py -dt 0.001 -nmax 10000 -i $< -o $@  > /dev/null

pickles/modLink.06.pickle.gz : pickles/cellMech.05.pickle.gz
	${CELLMECH}/modLink.py -i  $< -o $@ -c conf

pickles/cellMech.06.pickle.gz : pickles/modLink.06.pickle.gz
	${CELLMECH}/cellMech.py -dt 0.001 -nmax 10000 -i $< -o $@  > /dev/null

pickles/modLink.07.pickle.gz : pickles/cellMech.06.pickle.gz
	${CELLMECH}/modLink.py -i  $< -o $@ -c conf

pickles/cellMech.07.pickle.gz : pickles/modLink.07.pickle.gz
	${CELLMECH}/cellMech.py -dt 0.001 -nmax 10000 -i $< -o $@  > /dev/null

pickles/modLink.08.pickle.gz : pickles/cellMech.07.pickle.gz
	${CELLMECH}/modLink.py -i  $< -o $@ -c conf

pickles/cellMech.08.pickle.gz : pickles/modLink.08.pickle.gz
	${CELLMECH}/cellMech.py -dt 0.001 -nmax 10000 -i $< -o $@  > /dev/null

pickles/modLink.09.pickle.gz : pickles/cellMech.08.pickle.gz
	${CELLMECH}/modLink.py -i  $< -o $@ -c conf

pickles/cellMech.09.pickle.gz : pickles/modLink.09.pickle.gz
	${CELLMECH}/cellMech.py -dt 0.001 -nmax 10000 -i $< -o $@  > /dev/null

pickles/modLink.10.pickle.gz : pickles/cellMech.09.pickle.gz
	${CELLMECH}/modLink.py -i  $< -o $@ -c conf

pickles/cellMech.10.pickle.gz : pickles/modLink.10.pickle.gz
	${CELLMECH}/cellMech.py -dt 0.001 -nmax 10000 -i $< -o $@  > /dev/null

