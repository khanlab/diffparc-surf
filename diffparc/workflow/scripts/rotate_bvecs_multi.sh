#!/bin/bash

# from Washington-University/Pipelines/global/scripts/Rotate_bvecs.sh

#fix grep argument after pattern
unset POSIXLY_CORRECT

if [ "$#" -lt 4 ] ; then 
 echo "Usage: <original bvecs> <rotated (output) bvecs> <affine matrix index 0> ... <affine matrix index n-1> "
 echo ""
 echo "<affine matrix> is a FLIRT affine"
 echo ""
 exit 1;
fi

input=$1
output=$2
shift 2
matrices=( "$@" )
nmatrices=$#

echo "nmatrices $nmatrices"
if [ ! -e ${input} ] ; then
	echo "Source bvecs ${input} does not exist!"
	exit 1
fi

#get number of indices from bvecs:
numbvecs=`cat ${input} | head -1 | tail -1 | wc -w`
echo numbvecs $numbvecs

if [  "$nmatrices" -eq "$numbvecs" ]
then
    echo "correct number of matrices provided"
else
    echo "incorrect number of matrices provided, exiting"
    exit 1
fi

for i in `seq 0 $((nmatrices-1))`
do

  matrix=${matrices[$i]} 
  ii=$((i+1))



if [ ! -e ${matrix} ]; then
	echo "Matrix file ${matrix} does not exist!"
	exit 1
fi

m11=`avscale ${matrix} | grep Rotation -A 1 | tail -n 1| awk '{print $1}'`
m12=`avscale ${matrix} | grep Rotation -A 1 | tail -n 1| awk '{print $2}'`
m13=`avscale ${matrix} | grep Rotation -A 1 | tail -n 1| awk '{print $3}'`
m21=`avscale ${matrix} | grep Rotation -A 2 | tail -n 1| awk '{print $1}'`
m22=`avscale ${matrix} | grep Rotation -A 2 | tail -n 1| awk '{print $2}'`
m23=`avscale ${matrix} | grep Rotation -A 2 | tail -n 1| awk '{print $3}'`
m31=`avscale ${matrix} | grep Rotation -A 3 | tail -n 1| awk '{print $1}'`
m32=`avscale ${matrix} | grep Rotation -A 3 | tail -n 1| awk '{print $2}'`
m33=`avscale ${matrix} | grep Rotation -A 3 | tail -n 1| awk '{print $3}'`

echo "$m11"
numbvecs=`cat ${input} | head -1 | tail -1 | wc -w`
tmpout=${output}$$

    X=`cat ${input} | awk -v x=${ii} '{print $x}' | head -n 1 | tail -n 1 | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}' `
    Y=`cat ${input} | awk -v x=${ii} '{print $x}' | head -n 2 | tail -n 1 | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}' `
    Z=`cat ${input} | awk -v x=${ii} '{print $x}' | head -n 3 | tail -n 1 | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}' `
    rX=`echo "scale=7;  (${m11} * $X) + (${m12} * $Y) + (${m13} * $Z)" | bc -l`
    rY=`echo "scale=7;  (${m21} * $X) + (${m22} * $Y) + (${m23} * $Z)" | bc -l`
    rZ=`echo "scale=7;  (${m31} * $X) + (${m32} * $Y) + (${m33} * $Z)" | bc -l`

    if [ "${ii}" -eq 1 ];then
	echo $rX > ${output}; 
	echo $rY >> ${output};
	echo $rZ >> ${output}
    else
	cp ${output} ${tmpout}
	(echo $rX;echo $rY;echo $rZ) | paste ${tmpout} - > ${output}
    fi
    
#    let "ii+=1"
done

cat ${output} | awk '{for(i=1;i<=NF;i++)printf("%10.6f ",$i);printf("\n")}' > ${output}_
mv ${output}_ ${output}

rm -f ${tmpout}

