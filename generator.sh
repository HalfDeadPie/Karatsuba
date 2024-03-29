#!/bin/bash

DEGREE_LIMIT=20
COEF_LIMIT=20

if [ -z "$1" ]
then
    degree=$(( ( RANDOM % $DEGREE_LIMIT )  + 1 ))
else
    degree=$1
fi

mkdir data/deg$degree

file1=data/deg$degree/pol_$degree.1.txt
file2=data/deg$degree/pol_$degree.2.txt

touch $file1
touch $file2

echo "$degree" >> $file1
echo "$degree" >> $file2

for i in $( seq 1 $degree ) 
do
	echo -n "$(( ( RANDOM % $COEF_LIMIT )  + 1 )) " >> $file1
	echo -n "$(( ( RANDOM % $COEF_LIMIT )  + 1 )) " >> $file2
done

echo "" >> $file1
echo "" >> $file2

