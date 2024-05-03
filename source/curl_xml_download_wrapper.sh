#!/bin/bash

#aids_list="assay_list.list"
aids_list=$1

while read line; do
    aid=`echo $line | awk '{print $1}'` 
    echo "begin download AID: ${aid}"
    curl https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/${aid}/description/XML > ./assay_descriptions/${aid}.xml
    echo "completed download AID: ${aid}"
    echo ""
    echo ""
done < ${aids_list} > curl_down_aid_desc_xmls.log 2>&1 

