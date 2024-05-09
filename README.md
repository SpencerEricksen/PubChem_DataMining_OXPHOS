
# Pipeline for compiling PubChem BioAssay data for compounds from OXPHOS-related assays.

All data and model files are available on our corresponding data repo for this project: [https://zenodo.org/doi/10.5281/zenodo.11003005]:
   - Zenodo version v1 contains data/source for data mining workflow [https://doi.org/10.5281/zenodo.11003006]
   - Zenodo version v2 contains data/source for building machine learning classifiers [https://doi.org/10.5281/zenodo.11153594]

Here, we applied the pipeline to acquire compound testing data from PubChem. Data are pulled from assays likely related to the OXPHOS pathway. This includes a wide variety of assays involving different biochemical targets, whole cells, etc.  

The pipeline begins with a search using key word terms in the 'Assay Description' to identify assays potentially related to OXPHOS pathway. Here is the site for structuring and submitting your query:
https://www.ncbi.nlm.nih.gov/pcassay/advanced

The following query was performed Feb 22, 2022 at 1:57pm

Here is the search:
```
   ( "electron transport chain"[Assay Description] OR
     "mitochondrial complex"[Assay Description] OR
     "mitochondrial respiratory chain"[Assay Description] OR
     "mitochondrial membrane potential"[Assay Description] ) AND
    ( small_molecule[filt] )
```
This returned a list of 8415 AIDs that were downloaded ('pcassay_result.txt'). To download, results:
```
Click "Send to:', 
    Choose Destination: "File", 
    Format:"Summary (text)', 
    Sort by: "Default order"
    Click: "Create File"
```
This search included all chemical screening assays having any of the oxphos-related search terms in their \[Assay Description\]. 


Activate cheminformatics python env with conda package manager.
```
conda activate py38_chem
```

Download pcassay_result.txt as "pcaassay_result_summary_2022-02-22.txt".

Get AIDs from pcassay results--the following produces same output as a direct download of AID list "pcassay_result_UIlist_2022-02-22.txt".
```
grep AID: pcassay_result_summary_2022-02-22.txt | awk '{print $2}' > assay_list.list
```

Download assay data for assays in list.
```
# if you want direct access to scripts in your path
#export PATH="${PWD}/source:"$PATH
mkdir AIDs
nohup python ./source/get_bioassay_data_v1.4.py assay_list.list > get_bioassay_data_v1.4.log 2>&1 &
```
This writes data files for each assay in AID folder. Each assay data file looks like this "./AID/pcba-aidXXXXX.csv".


Get the number of cpds tested in each assay.
```
cat get_bioassay_data_v1.4.log | grep AID | awk '{print $1, $3}' | tr ":" " " | sort -n -k4,4 -r | awk '{print $2","$4}' > mols_tested.csv
```
"mols_tested.csv" looks like this:
```
AID,mol_total
1465,215402
540299,103205
720637,10486
...
```

Merge all assay data, keep relevant data columns, dump as a pickle file. Pickle (.pkl) maintains data types so better than CSV!
```
python ./source/merge_aid_records.py
```
in: "AIDs/pcba-aid*.csv"   
out: "oxphos_merged_aids_records.pkl"


Build xml files to use for bulk query of all unique CIDs in merged assay data set.
```
python ./source/write_cgi_for_cid_fetch_v1.4.py
```
   in: "oxphos_merged_aids_records.pkl"  gets the set of all assayed CIDs   
   out: "pc_fetch_all_0.cgi", "pc_fetch_all_1.cgi" (pc_fetch query files for all CIDs in 2 batches)


Run the SMILES fetches for all unique CIDs. This is done in 2 batches due to a limit on number of CIDs per request.
```
./source/wrapper_fetch_v1.2.sh pc_fetch_all_0.cgi
gunzip 9245063668458039.txt.gz

./source/wrapper_fetch_v1.2.sh pc_fetch_all_1.cgi
gunzip 2267636078691197203.txt.gz
```

Merge all CID smiles into a single file.
```
cat 9245063668458039.txt  2267636078691197203.txt > pcba_oxphos_all_cids.smi
```

Add SMILES to merged AIDs data. This includes original, rdkit-canonicalized, and -desalted SMILES.
```
python ./source/merge_assaydata_smiles.py
```
  in: "oxphos_merged_aids_records.pkl", "fetch_CGIs_bulk/pcba_oxphos_all_cids.smi"   
 out: "oxphos_merged_aids_cids_clnsmi.pkl"


Get assay descriptions and metadata associated with the AIDs (.xml file for each AID).
```
./source/curl_xml_download_wrapper.sh assay_list.list
```
  in: "assay_list.list"    
 out: "./assay_descriptions/XXXXX.xml"    (file for each AID XXXXX)   

Extract relevant fields from the downloaded .xml AID descriptions.
```
for f in `cat assay_list.list`; do python ./source/parse_pcba_AIDs_desc_xml_v1.4.py $f; done > all_assays_desc.csv
```

Clean up the file "all_assays_desc.csv".
```
sed -i '2,${/^AID/d;}' ./assay_descriptions/all_assays_desc.csv
```

Merge assay descriptions into the compound activity dataframe.
```
python ./source/merge_assaydesc_v1.2.py
```
  in: "oxphos_merged_aids_cids_clnsmi.pkl", "./assay_descriptions/all_assays_desc.csv"    
 out: "oxphos_merged_aids_cids_clnsmi_assaydesc.pkl"    


Add flag for compounds from AIDs having ETC-linked terms in assay Names, Titles, or Abstracts.
```
python ./source/add_flag_ETC-linked_assay_desc.py
```
  in: "oxphos_merged_aids_cids_clnsmi_assaydesc.pkl"   
 out: "oxphos_merged_aids_cids_clnsmi_assaydesc_ETCassay.pkl"   


Add flag for compounds from AIDs having an ETC-linked PMID (paper).
```
python ./source/add_flag_ETC-linked_PMID.py
```
  in: "oxphos_merged_aids_cids_clnsmi_assaydesc_ETCassay.pkl"   
 out: "oxphos_merged_aids_cids_clnsmi_assaydesc_ETCassay_pmid.pkl"   


Isolate the molecule set for descriptor generations and pains/tox flags, save dataframe as a pickle (.pkl).
```
python get_unique_pubchem_cids.py
```
  in: "oxphos_merged_aids_cids_clnsmi_assaydesc_ETCassay_pmid.pkl"   
 out: "unique_pubchem_cids.pkl"   


Determine Bemis-Murcko scaffolds and reduced Murcko scaffolds (carbon frames).
```
python ./source/get_scaffolds.py
```
  in: "unique_pubchem_cids.pkl"   
 out: "unique_pubchem_cids_scaffolds.pkl"   


Get rdkit pains flags for compound set. Note: this one takes a while--maybe 20 minutes?
```
nohup python ./source/rdkit_add_pains_flags.py unique_pubchem_cids.pkl unique_pubchem_cids_pains.pkl > rdkit_add_pains_flags.log 2>&1 &
```
  in: "unique_pubchem_cids.pkl"   
 out: "unique_pubchem_cids_pains.pkl"   


Get morgan fingerprint features (ECFP6, length=2048) for compound set.
```
nohup python source/rdkit_get_morgan_fps.py > rdkit_get_morgan_fps.log 2>&1 &
```
  in: "unique_pubchem_cids.pkl"   
 out: "unique_pubchem_cids_fps.pkl"   


Get RDKit molecular descriptors for compound set.
```
nohup python source/rdkit_calc_mol_descriptors.py > rdkit_calc_mol_descriptors.log 2>&1 &
```
  in: "unique_pubchem_cids.pkl"   
 out: "unique_pubchem_cids_desc.pkl"   


Get NPscores (natural product likeness scores) for compound set.
```
python ./source/npscorer_v1.2.py unique_pubchem_cids.pkl unique_pubchem_cids_npscores.pkl
```
  in: "unique_pubchem_cids.pkl"   
 out: "unique_pubchem_cids_npscores.pkl"   


Get "Lipinski" flags for compound set. These PK flags were actually based on our own bespoke parameter conditions.
```
python ./source/rdkit_add_lipinski_flags.py
```
  in: "unique_pubchem_cids.pkl"   
 out: "unique_pubchem_cids_lipinski.pkl"   


Get complete assay results for each molecule in compound set.
```
python ./source/compile_add_bioassay_results.py
```
  in: "oxphos_merged_aids_cids_clnsmi_assaydesc_ETCassay_pmid.pkl", "unique_pubchem_cids.pkl"   
 out: "unique_pubchem_cids_complete_assay_results.pkl"   


Assign activity labels to the compound set.
```
python ./source/assign_cpd_oxphos_activity_labels.py
```
  in: "unique_pubchem_cids_complete_assay_results.pkl"   
 out: "unique_pubchem_cids_complete_assay_results_w_labels.pkl"   


Merge compound data into one dataframe.
```
python ./source/merge_cpd_data.py
```
  in: "unique_pubchem_cids_complete_assay_results_w_labels.pkl", "unique_pubchem_cids_desc.pkl",    
  in: "unique_pubchem_cids_fps.pkl", "unique_pubchem_cids_lipinski.pkl", "unique_pubchem_cids_npscores.pkl",   
  in: "unique_pubchem_cids_pains.pkl", "unique_pubchem_cids_scaffolds.pkl"   
 out: "unique_pubchem_cids_all.pkl"   


Add cluster labels to the active cpds based on HAC using Jaccard distance matrix on ECFP6 fingerprints. Z-cutoff for clusters was 0.750.
```
python ./source/cluster_oxphos_actives_scipy_HAC.py
```
  in: "unique_pubchem_cids_all.pkl"   
 out: "unique_pubchem_cids_all_wclusts.pkl"   


Build training data with labeled actives/inactives.
```
python ./source/build_training_data_sets.py
```

Draw RDKit mol grids for each cluster for inspection of various active chemotypes. Only the active compounds are clustered.
```
python ./source/rdkit_draw_clusters_aln_oxphos_v1.3.py 0.750
```
Draw mol grid for cluster medoids from clusters with 3 or more compounds.
```
python ./source/rdkit_draw_cluster_reps_ge3_ETC_landscape_v1.6.py
python ./source/rdkit_draw_cluster_reps_ge3_ETC_portrait_v1.6.py
```

Draw mol grid for all cluster medoids.
```
python ./source/rdkit_draw_cluster_reps_all_ETC_landscape_v1.6.py
python ./source/rdkit_draw_cluster_reps_all_ETC_portrait_v1.6.py
```
