
import sys
import glob
import xml.etree.ElementTree as ET


aid = sys.argv[1]

def get_descriptions( aid ):
    '''extract the relevant attribute data from the xml record for the assay (AID)'''

    tree = ET.parse( "./assay_descriptions/"+aid+".xml" )
    root = tree.getroot()

    aid_desc_dict = {}

    for desc in root.iter('{http://www.ncbi.nlm.nih.gov}PC-AssayDescription_description_E'):
        try:
            k = desc.text.split(':')[0].strip()
            v = desc.text.split(':')[1].strip()
            aid_desc_dict[k] = v
        except:
            pass

    for desc in root.iter('{http://www.ncbi.nlm.nih.gov}PC-AssayDescription_comment_E'):
        try:
            k = desc.text.split(':')[0].strip()
            v = desc.text.split(':')[1].strip()
            aid_desc_dict[k] = v
        except:
            pass

    for t in [ 'PC-AssayDescription_name', 'PC-ID_id', 'PC-XRefData_pmid', 'PC-XRefData_gene', 'PC-XRefData_protein-gi' ]:
        k = t.split('_')[-1]
        if '{http://www.ncbi.nlm.nih.gov}'+t in [ i.tag for i in root.iter() ]:
            try:
                for desc in root.iter( '{http://www.ncbi.nlm.nih.gov}'+t ):
                    v = desc.text
            except:
                v = None 
                continue
        else:
            v = None

        aid_desc_dict[k] = v

    return aid_desc_dict



d = get_descriptions(aid)

# going with the '|' delimiter to avoid issues with commas in the abstracts
print('AID|pmid|DOI|Year|ChEMBL_Target_Name|ChEMBL_Target_Type|Target_ChEMBL_ID|Confidence|Relationship_Type|Title|name|Abstract')

keys_list = ['id','pmid','DOI','Year','ChEMBL Target Name','ChEMBL Target Type','Target ChEMBL ID','Confidence','Relationship Type','Title','name','Abstract']

string_list = []
for k in keys_list:
    if k in d:
        string_list.append( str(d[k]) )
    else:
        string_list.append("")

string = "|".join( [ '"'+str(i)+'"' for i in string_list] )
print(string)

