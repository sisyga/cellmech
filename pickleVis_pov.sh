T=$(tempfile); 
R=$(cat $2)
view=${3:-v}	#default is v
$(dirname $0)/pickleVis_pov.py -i $1 -l 2 -r $R -w $view> $T; 
povray -d -GD -GR -GS +I$T +H480 +W480 -O-  2> /dev/null
rm $T
