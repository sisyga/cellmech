
python testplane.py -o plane.pickle.gz 
python stretch.py   -i plane.pickle.gz -o stretched-plane.pickle.gz -m 10.0
python genMake.py   -i stretched-plane.pickle.gz -n 10 > 10.make

