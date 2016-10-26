
# generating a monolayer of particles as test configuration: 
python testplane.py -o plane.pickle.gz 

# obtaining stretched configuration:
python stretch.py   -i plane.pickle.gz -o stretched-plane.pickle.gz -m 10.0

# creating makefile to run simulation:
python genMake.py   -i stretched-plane.pickle.gz -n 10 > 10.make

# running simulation:
make -f 10.make
