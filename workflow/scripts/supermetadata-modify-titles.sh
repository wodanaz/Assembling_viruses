#!/usr/bin/env bash
set -x

if [ "$EVMODE" == "campus" ]
then
    # campus mode
    cat results/$PROJECTNAME.csv | sed -r 's/,/\t/g' | sed -r 's/_NC_045512//g' | awk '{print $1 "\t" $2 }' | sed '1d' > variants.tab
fi

if [ "$EVMODE" == "hospital" ]
then
    # hospital mode
    cat results/${PROJECTNAME}_lineages_of_concern.csv | sed -r 's/,/\t/g' | sed -r 's/\|/\t/g' | awk '{print $1 "\t" $3 }' | sed '1d' > variants.tab
fi

# both hospital and campus mode
sed -r 's/.bqsr//g' results/coverage.gatk.tab | sed '1d' | awk '{print $1 "\t" $7 }'  | sort -k1d  > coverage.sort.tab

join -j 1 <(sort coverage.sort.tab) <( sort $DATETAB) | awk '{ print $1 "\t" $2 "\t" $3 }' > metadata1.tab

join -j 1 <(sort variants.tab) <( sort metadata1.tab ) | awk '{ print $1 "$" $2 "$" $3 "$" $4 }' > results/supermetadata.tab

for i in `cat results/supermetadata.tab`; do
    root=`echo ${i} | cut -d'$' -f 1`;
    rootA=`echo ${i} | cut -d'$' -f 2`;
    rootB=`echo ${i} | cut -d'$' -f 3`;
    rootC=`echo ${i} | cut -d'$' -f 4`;
    result=${rootB/.*}
    if [ "$result" -ge 90  ] ; then  # Here I use a filter to only modify fasta files with a coverage higher than 90%
        sed -r "s/_NC_045512/|${rootA}|${rootB}|${rootC}/g" $root.cleaned.fasta > $root.info.fasta ;
    fi ;
done

cat *info.fasta > results/${PROJECTNAME}.fasta
