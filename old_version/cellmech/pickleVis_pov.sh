T=$(tempfile); 
R=$(cat $2)
if [ $# -gt 2 ] ; then
	CS=" -cs $(cat $3) "
else
	CS=""
fi
$(dirname $0)/pickleVis_pov.py -i $1 -l 1 -r $R $CS > $T; 
povray -d -GD -GR -GS +I$T +H480 +W480 -O-  2> /dev/null
rm $T
