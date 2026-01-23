#!/bin/bash

# Script for lata comparison used in AlgoQC test cases
# Done because algoqc has non standard post format

lata=$(ls ${1}_lata_*.sauv.lata 2>/dev/null)
# echo $lata
[ "$lata" = "" ] && echo "No lata files generated. Test Failed" && exit 1

# Trouve le nom du cas test:
fullname=$1

name=${fullname#PAR_}

echo "Test case name: ${fullname}"

# untar ref for compare_lata
if [ -f "ref.tar.gz" ] ; then 
    tar xzf "ref.tar.gz"
else
    UPDATE_LATA_REF=1 # if ref does not exist, generate it
fi


# Update ref if asked
if [ "$UPDATE_LATA_REF" != "" ] && [ "$fullname" == "$name" ]
then

	test_path=${TrioCFD_project_directory}/tests/Reference/$(realpath -s --relative-to=$TRUST_TMP/tests $(pwd))

    echo $(pwd)
    echo ${test_path}
    
    echo "Updating ref"
    [ -f ref ] && rm -r ref
    mkdir -p ref
    cp ${name}_lata_*.sauv.lata* ref/
    tar czf ref.tar.gz ref

    # check ref size (max 15 MB, arbitrary)
	max_size=15000000
	ref_size=$(wc -c <"ref.tar.gz")
	if [ $ref_size -ge $max_size ]; then
		echo "Error: reference size too big : $ref_size"
		exit 1
	fi

    cp ref.tar.gz ${test_path}
fi

  [ $? != 0 ] && echo "failed" && exit $?

# Compare to ref
if [ -f ${name}_lata_0.sauv.lata ] ; then 


    for f in $lata 
    do
        echo "Compare $f with reference..."
        compare_lata $f ref/$f
        [ $? != 0 ] && echo "failed" && exit 1
    done    
else
    echo "${name}_lata_0.sauv.lata not found" && exit 1
fi

exit 0