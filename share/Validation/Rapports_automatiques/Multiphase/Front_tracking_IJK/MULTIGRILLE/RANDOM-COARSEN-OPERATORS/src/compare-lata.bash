#!/bin/bash

# Check if compare_lata is available
if ! type compare_lata 1>/dev/null 2>&1; then
    printf "%s\n" "You must source TRUST first!"
    return
fi

ref="1-LEVEL/COARSEN-111/rcu-111.lata"
repocompar1="1-LEVEL/COMPARELATA"
repocompar2="2-LEVEL/COMPARELATA"
repocompar3="3-LEVEL/COMPARELATA"
mkdir -p "${repocompar1}" "${repocompar2}" "${repocompar3}"

function compare-parallel {
    echo "Run 1-LEVEL"
    lesnums=$(find 1-LEVEL/ -mindepth 1 -type d -name "COARSEN*" | sort | sed "s|1-LEVEL/COARSEN-||g")
    parallel "compare_lata \"${ref}\" \"1-LEVEL/COARSEN-{}/rcu-{}.lata\" > \"${repocompar1}/compare-lata-{}.log\"" ::: "${lesnums}"
    (check-diff "${repocompar1}")

    echo "Run 2-LEVEL"
    lesnums=$(find 2-LEVEL/ -mindepth 1 -type d -name "COARSEN*" | sort | sed "s|2-LEVEL/COARSEN-||g")
    parallel "compare_lata \"${ref}\" \"2-LEVEL/COARSEN-{}/rcu-{}.lata\" > \"${repocompar2}/compare-lata-{}.log\"" ::: "${lesnums}"
    (check-diff "${repocompar2}")

    echo "Run 3-LEVEL"
    lesnums=$(find 3-LEVEL/ -mindepth 1 -type d -name "COARSEN*" | sort | sed "s|3-LEVEL/COARSEN-||g")
    parallel "compare_lata \"${ref}\" \"3-LEVEL/COARSEN-{}/rcu-{}.lata\" > \"${repocompar3}/compare-lata-{}.log\"" ::: "${lesnums}"
    (check-diff "${repocompar3}")
}

function compare-loop {
    echo "Run 1-LEVEL"
    lesnums=$(find 1-LEVEL/ -mindepth 1 -type d -name "COARSEN*" | sort | sed "s|1-LEVEL/COARSEN-||g")
    for num in ${lesnums}
    do
        compare_lata "${ref}" "1-LEVEL/COARSEN-${num}/rcu-${num}.lata" > "${repocompar1}/compare-lata-${num}.log"
    done
    (check-diff "${repocompar1}")

    echo "Run 2-LEVEL"
    lesnums=$(find 2-LEVEL/ -mindepth 1 -type d -name "COARSEN*" | sort | sed "s|2-LEVEL/COARSEN-||g")
    for num in ${lesnums}
    do
        compare_lata "${ref}" "2-LEVEL/COARSEN-${num}/rcu-${num}.lata" > "${repocompar2}/compare-lata-${num}.log"
    done
    (check-diff "${repocompar2}")

    echo "Run 3-LEVEL"
    lesnums=$(find 3-LEVEL/ -mindepth 1 -type d -name "COARSEN*" | sort | sed "s|3-LEVEL/COARSEN-||g")
    for num in ${lesnums}
    do
        compare_lata "${ref}" "3-LEVEL/COARSEN-${num}/rcu-${num}.lata" > "${repocompar3}/compare-lata-${num}.log"
    done
    (check-diff "${repocompar3}")
}

function check-diff {
    cd "${1}" || return
    \grep "Number of differences :" ./*.log | \grep -v "0$" > non-zero-diff.log
}

function main {

    # only does the work if this file does not exists
    # allows to skip parsing the lata files during archive mode of the validation report
    # if changing the name, also change it in the notbook, in the cell where this script is called
    FILE_TEST="compare_lata_done"

    if [ ! -f "${FILE_TEST}" ]; then
        if type parallel 1>/dev/null 2>&1; then
            echo "We use the GNU parallel power!"
            compare-parallel
        else
            compare-loop
        fi
    else 
        echo "compare_lata already done, not running again"
    fi

    touch ${FILE_TEST}

}

main
