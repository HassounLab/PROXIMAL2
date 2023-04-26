# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 12:00:02 2022

@author: fbalzerani
"""

# First we load libraries and modules.
import requests
import pandas as pd
from rdkit import Chem
import sqlite3
import re
from bioservices import KEGG
import pubchempy as pcp
import os

# =============================================================================
# The path should be:
# pwd = './PROXIMAL2/'
# =============================================================================

# Function definition

def ExtractFromPubchem(moleculeName:str) -> str:
    compound = pcp.get_compounds(moleculeName,'name')
    return compound[0].isomeric_smiles
        
def ExtractFromHMDB(HMDB_ID:str, hmdb_server) -> str:
    response = requests.request('GET',hmdb_server+HMDB_ID+'.xml')
    tmp_response = response.text
    pattern = "<smiles>(\S+)</smiles>"
    return re.findall(pattern,tmp_response)[0]
    
    
def ExtractMOLFromMetanetX(molecule_ID:str, metanetX_server) -> (str,str):
    response = requests.request('GET',metanetX_server+"chem_info/"+molecule_ID)
    tmp_response = response.text
    
    pattern_name = '"name": "([\S\s]+)",\s"url'
    pattern_smiles = '"smiles": \["(\S+)"\],'
    name = re.findall(pattern_name,tmp_response)[0]
    smiles = re.findall(pattern_smiles,tmp_response)[0]
    return (name,smiles)

def ExtractRXNFromMetanetX(rxn_ID:str, metanetX_server) -> str:
    response = requests.request('GET',metanetX_server+"equa_info/"+rxn_ID)
    tmp_response = response.text
    # check if we can reach directly the inchi
    pattern = '<tr><td>equation</td><td>([\S\s]+)</td></tr>\n<tr><td>is balanced?'
    tmp = re.findall(pattern,tmp_response)[0]
    pattern_id = "<br>([\S\s]+)"
    tmp_id = re.findall(pattern_id, tmp)[0].split("=")
    ids = [re.findall("MNXM\d+", x) for x in tmp_id]
    names = []
    for each in ids:
        tmp_name = []
        for id in each:
            tmp_name.append(ExtractMOLFromMetanetX(id, metanetX_server)[0])
        names.append(tmp_name)
    pattern_ec = "<tr><td>EC number</td><td>(\S+)</td></tr>\n"
    ec = re.findall(pattern_ec,tmp_response)
    ec = ec[0].split("<br>")
    for i in range(len(ec)):
        while len(ec[i].split(".")) < 4:
            tmp_ec = ec[i].split(".")
            tmp_ec.append("-")
            ec[i] = '.'.join(tmp_ec)
    return (names,ids,list(set(ec)))

def Convert(tup,di):
    di = dict(tup)
    return di

# =============================================================================
# =============================================================================

# Define the input. 
# csv file with the reactions of interest. Tabulator as separator. 
# Must be present "id", "formula", "EC" columns.
reactions = pd.read_csv("./input/XXXX.csv", sep = "\t")

# csv file with the metabolites that can be involved in the reactions. Tabulator as separator.
# The file must have a column per each following information (even if empty):
# name, hmdb, kegg, metanetx
# If the file is empty or a metabolite in the reactions is not present in the file,
# it will be used just the name from the reaction and the eventual information present in RetroRules,
# otherwise it will be excluded

metabolites = pd.read_csv("./input/XXXX.csv", sep = "\t")
# =============================================================================
# =============================================================================

# We specify the root URL.
kegg = KEGG(verbose = False)
hmdb_server = 'http://www.hmdb.ca/metabolites/'
metanetX_server = "https://www.metanetx.org/"

# =============================================================================
# =============================================================================

# Use the RetroRules DB as source of information

# Connession to the DB
conn = sqlite3.connect('./input/retrorules_dump/mvc.db')
c = conn.cursor()

c.execute("SELECT * FROM chemical_species")

# Extract name-ID metabolite relation

metabolites_db = c.fetchall()

# Extract all the information about metabolite in RetroRules

mets_df = pd.DataFrame(index = range(len(metabolites_db)), columns = ['ID', 'Name', 'cpd_ID', 'HMDB_ID', 'KEGG_ID'])
mets_df.loc[:,'ID'] = [x[0] for x in metabolites_db]
mets_df.loc[:,'Name'] = [x[1] for x in metabolites_db]
mets_df.loc[:,'cpd_ID'] = [x[12] for x in metabolites_db]
mets_df.loc[:,'KEGG_ID'] = [x[7] for x in metabolites_db]
mets_df.loc[:,'HMDB_ID'] = [x[6] for x in metabolites_db]

# =============================================================================
# =============================================================================
  
# Extract the cofactors to be removed from the reactions

cofactor = pd.read_csv('./input/cofactors.csv')
cofactor.drop(columns = "inchi", inplace = True)
cofactor.drop_duplicates(inplace = True)

# Add manually some cofactors

to_add = pd.Series(['molecular entity', 'H(+)', 'H2O','chlorite', "O2",
                    'CO2','NADPH','NADP(+)','NADH','NAD(+)','', 'e(-)',
                    'Zn(2+)','phosphate','H2O2','Mg(2+)','Methane','Fe(3+)',
                    'UDP-D-glucose','UDP-L-rhamnose','Co(2+)', 'S',
                    'coenzyme M','CoA','acryloyl-CoA','UDP-alpha-D-galactose',
                    'phosphonate','UDP-alpha-D-xylose','Na(+)','chloride',
                    'Ca(2+)','Chloride','H2','malonyl-CoA','methane',
                    'benzoyl-CoA','formyl-CoA','hydrogen sulfide',
                    'Nicotinamide adenine dinucleotide - reduced',
                    'Nicotinamide adenine dinucleotide','NH4(+)',
                    'NADH-P-OR-NOP','NAD-P-OR-NOP','bromide','diphosphate'])

cofactor = cofactor['name'].append(to_add, ignore_index = True)

# Extract substrates and products of the reactions

reactions['Sub_name'] = [x.split(" -> ")[0].split(" + ") for x in reactions.formula]
reactions['Prod_name'] = [x.split(" -> ")[1].split(" + ") for x in reactions.formula]

# Remove cofactors from substrates and products

to_rem = []
for idx in reactions.index:
    met_del = []
    for x in reactions['Sub_name'][idx]:
        if x in cofactor.values:
            met_del.append(x)
    
    for x in met_del:
        # finalReactions['Sub_ID'][idx].pop(finalReactions['Sub_name'][idx].index(x))
        reactions['Sub_name'][idx].pop(reactions['Sub_name'][idx].index(x))
    if len(reactions['Sub_name'][idx]) == 0:
        to_rem.append(idx)
        continue
    
    met_del = []
    for x in reactions['Prod_name'][idx]:
        if x in cofactor.values:
            met_del.append(x)
    for x in met_del:
        reactions['Prod_name'][idx].pop(reactions['Prod_name'][idx].index(x))
    if len(reactions['Prod_name'][idx]) == 0:
        to_rem.append(idx)
        continue
    
reactions.drop(index = to_rem, inplace = True)

# Extract the list of metabolites involved in the previous reactions

met = set()
for idx  in reactions.index:
    for x in reactions['Sub_name'][idx]:
        met.add(x)
    for x in reactions['Prod_name'][idx]:
        met.add(x)

# Extract the InChI for each molecule

if os.path.isfile("./input/reachableMolecules.csv"):
    os.remove("./input/reachableMolecules.csv")

AllMolecules = pd.DataFrame(columns = ['name','smiles'], index = range(len(met)))
AllMolecules['name'] = list(met)

# We now search for each of the molecules in the list.
not_done = []

for i in AllMolecules.index:
    print("Annotating metabolites... : %d/%d" % (i,len(AllMolecules)-1))
    name_mol = AllMolecules.name[i]
    info_mol = metabolites.loc[metabolites.name.isin([name_mol]),]
    try:
        # Use the name to go through pubchem
        smiles = ExtractFromPubchem(name_mol)
        if smiles == "":
            raise Exception('Empty smiles')
        AllMolecules.loc[i,'smiles'] = smiles
    except:
        try:
            # try the metanetx id
            id = info_mol.metanetx.values[0]
            smiles = ExtractMOLFromMetanetX(id, metanetX_server)[1]
            if smiles == "":
                raise Exception('Empty smiles')
            AllMolecules.loc[i,'smiles'] = smiles
        except:
            try:
                # try the kegg id
                id = info_mol.kegg.values[0]
                text = kegg.get(id,option = 'mol')
                smiles = Chem.MolToSmiles(Chem.MolFromMolBlock(text))
                if smiles == "":
                    raise Exception('Empty smiles')
                AllMolecules.loc[i,'smiles'] = smiles
            except:
                try:
                    # try the hmdb id
                    id = info_mol.hmdb.values[0]
                    smiles = ExtractFromHMDB(id, hmdb_server)
                    if smiles == "":
                        raise Exception('Empty smiles')
                    AllMolecules.loc[i,'smiles'] = smiles
                except:
                    # try to extract info from retrorules
                    fromRetroRules = mets_df[mets_df['Name'] == name_mol]
                    try: 
                        # extract the corresponding name and go through pubchem
                        id = fromRetroRules.Name.values[0]
                        smiles = ExtractFromPubchem(AllMolecules.name[i])
                        if smiles == "":
                            raise Exception('Empty smiles')
                        AllMolecules.loc[i,'smiles'] = smiles
                    except:
                        try:
                            # extract the corresponding HMDB ID and go through HMDB
                            id = fromRetroRules.HMDB_ID.values[0]
                            smiles = ExtractFromHMDB(id, hmdb_server)
                            if smiles == "":
                                raise Exception('Empty smiles')
                            AllMolecules.loc[i,'smiles'] = smiles
                        except:
                            
                            try:
                                # extract the corresponding KEGG ID and go through KEGG
                                id = fromRetroRules.KEGG_ID.values[0]
                                text = kegg.get(id,option = 'mol')
                                smiles = Chem.MolToSmiles(Chem.MolFromMolBlock(text))
                                if smiles == "":
                                    raise Exception('Empty smiles')
                                AllMolecules.loc[i,'smiles'] = smiles
                            except:
                                not_done.append(i)
 
AllMolecules.drop(not_done, inplace = True)

for idx in AllMolecules.index:
    if re.findall("\.", AllMolecules.smiles[idx]):
        AllMolecules.loc[idx,'smiles'] = AllMolecules.smiles[idx].split(".")[0]
        
AllMolecules['inchi'] = [Chem.MolToInchi(Chem.MolFromSmiles(x)) for x in AllMolecules.smiles]
        
AllMolecules.to_csv('./input/reachableMolecules.csv', sep=',', index=False)

to_rem = []
for idx  in reactions.index:
    new = []
    for x in reactions['Sub_name'][idx]:
        try:
            if x in AllMolecules['name'].values:
                new.append(x)
        except KeyError:
            continue
    reactions.loc[idx,'Substrates'] = ' + '.join(new)
    if len(new) == 0:
        to_rem.append(idx)
        continue
    new = []
    for x in reactions['Prod_name'][idx]:
        try:
            if x in AllMolecules['name'].values:
                new.append(x)
        except KeyError:
            continue
    if len(new) == 0:
        to_rem.append(idx)
        continue
    reactions.loc[idx,'Products'] = " + ".join(new)

# Remove reactions that were made of just cofactors (not interesting)

reactions.drop(index = to_rem, inplace = True)

# Generate the entire formula

reactions['namesFormula'] = reactions[['Substrates', 'Products']].agg(' -> '.join, axis=1)

# Save the dataframe
if os.path.isfile("./input/templateReactions.csv"):
    os.remove("./input/templateReactions.csv")
    reactions.to_csv("./input/templateReactions.csv", index = False)
else:
    reactions.to_csv("./input/templateReactions.csv", index = False)

