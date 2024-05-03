#!/bin/bash

# usage: ./wrapper.sh   pc_fetch.cgi.xml

# where the 1st argument is the .xml file obtained from
# the https://pubchem.ncbi.nlm.nih.gov/pc_fetch/pc_fetch.cgi
# where I submitted a list of CIDs to put into my query for
# downloading (structures in smiles format)

###########################

fetchfile=$1
#url="https://pubchem.ncbi.nlm.nih.gov/pc_fetch/pc_fetch.cgi"
url="https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"

#xml="pc_fetch_aid686979_0.cgi"
aid=`echo ${fetchfile} | cut -f3 -d"_" | cut -c4-`
aid_chunk=`echo ${fetchfile} | cut -f5 -d"_" | cut -f1 -d"."`

### submit query

# curl -X POST -d @filename http://hostname/resource
curl -X POST -d @${fetchfile} $url > submitquery_${aid}_${aid_chunk}_log.html 2>&1
echo "Query XML: "${fetchfile} >> submitquery_${aid}_${aid_chunk}_log.html


### if result is ready to download, proceed
reqid=`cat submitquery_${aid}_${aid_chunk}_log.html | grep "Waiting_reqid" | cut -d">" -f2 | cut -d"<" -f1`
if [ -z ${reqid} ]; then
    # get fetch path if it is available in submitquery html file and fetch it
    reqid=`cat submitquery_${aid}_${aid_chunk}_log.html | grep fetch/ | cut -f7 -d"/" | cut -f1 -d"."`
    echo "reqid:${reqid} is ready, downloading smiles for pcba-aid${aid}_${aid_chunk}"
    remotefilepath=`grep fetch submitquery_${aid}_${aid_chunk}_log.html | cut -d">" -f2 | cut -d"<" -f1`
    wget $remotefilepath > download_queryresults_${aid}_${aid_chunk}.log 2>&1

### otherwise, poll the query to check progress, wait 40sec, then attempt download
else
    echo "reqid:${reqid} isn't ready, waiting 40 seconds then will query the poll..."
    # note, you will need a template poll.cgi file available for next line
    cat poll.cgi | sed "s/REQID_NUM/${reqid}/g" > poll_${reqid}_${aid}_${aid_chunk}.cgi
    xml=poll_${reqid}_${aid}_${aid_chunk}.cgi
    sleep 40 
    # curl -X POST -d @filename http://hostname/resource
    curl -X POST -d @$xml $url > pollquery_${aid}_${aid_chunk}_log.html 2>&1
    # now clean up after poll request
    #rm $xml
    echo "reqid:${reqid} was polled, results in pollquery_${aid}_${aid_chunk}_log.html"
    remotefilepath=`grep fetch pollquery_${aid}_${aid_chunk}_log.html | cut -d">" -f2 | cut -d"<" -f1`
    wget $remotefilepath > download_queryresults_${aid}_${aid_chunk}.log 2>&1
fi

echo "downloading smiles...."
#sleep 3
#gunzip ${reqid}.txt.gz
#sleep 5
#mv "$reqid".txt pcba-aid${aid}_${aid_chunk}.smi


