# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 15:09:19 2022

@author: fbalzerani
"""

import pandas as pd
import re
from bioservices import KEGG
from rdkit import Chem
import os
from rdkit.Chem import AllChem, DataStructs

cofactors = pd.read_csv("../input/cofactors.csv")

m_cof = [Chem.MolFromInchi(x) for x in cofactors.inchi]
f_cof = [AllChem.GetMorganFingerprint(x,2) for x in m_cof]

KEGG_FILENAMES = {
    'reaction': '../data/reaction',
    'rpair': '../data/rpair',
    'enzyme': '../data/enzyme'
}

f = open(KEGG_FILENAMES['reaction'], 'r')
reactions_text = f.read().split("///")

reactions = {}
enzymes = {}
for each in reactions_text:
    if re.findall("RPAIR", each):
        key_rxn = re.findall("R\d+", each)[0]
        tmp_enzyme = list(set(re.findall("\d+\.\d+\.\d+\.\d+", each)))
        if len(tmp_enzyme) == 0:
            continue
        try:
            tmp = re.split("\n",re.findall("RPAIR([\S\s]+)ENZYME", each)[0])
        except IndexError:
            try:
                tmp = re.split("\n",re.findall("RPAIR([\S\s]+)PATHWAY", each)[0])
            except IndexError:
                tmp = re.split("\n",re.findall("RPAIR([\S\s]+)", each)[0])
        values_pairs = [re.findall("C\d+_C\d+", x)[0] for x in tmp if "main" in x]
        if len(values_pairs) == 0:
            continue
        reactions[key_rxn] = list(set(values_pairs))
        enzymes[key_rxn] = tmp_enzyme
            
k = KEGG()

test_reactions = pd.DataFrame(columns = ["id","namesFormula"])
test_metabolites = pd.DataFrame(columns = ['name','id','inchi','smiles'])
pos_rxn = 0
pos_met = 0
for key,value in reactions.items():
    print("Annotating metabolites for each reaction... %d/%d" % (pos_rxn,len(reactions)-1))
    # test_reactions.loc[pos_rxn,"id"] = key
    for each_v in value:
        tmp_formula = each_v.split("_")
        
        flag = True
        # test_reactions.loc[pos_rxn,"namesFormula"] = " -> ".join(tmp_formula)
        for comp in tmp_formula:
            if not comp in test_metabolites['id'].values:
                mol = k.get(comp,'mol')
                test_metabolites.loc[pos_met,'id'] = comp
                test_metabolites.loc[pos_met,'name'] = comp
                try:
                    tmp_inchi = Chem.MolToInchi(Chem.MolFromMolBlock(mol))
                    tmp_smiles = Chem.MolToSmiles(Chem.MolFromMolBlock(mol))
                    if tmp_inchi == "":
                        flag = False
                    else: 
                        f_met = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(tmp_smiles), 2)
                        if [x for x in cofactors.index if DataStructs.DiceSimilarity(f_met,f_cof[x]) > 0.95]:
                            flag = False
                            continue
                        test_metabolites.loc[pos_met,'inchi'] = tmp_inchi
                        test_metabolites.loc[pos_met,'smiles'] = tmp_smiles
                except:
                    flag = False
                pos_met += 1
        if flag:
            test_reactions.loc[pos_rxn,"id"] = key
            test_reactions.loc[pos_rxn,"namesFormula"] = " -> ".join(tmp_formula)
            test_reactions.loc[pos_rxn,'EC'] = enzymes[key]
            pos_rxn += 1
    
test_reactions.dropna(inplace = True)
test_metabolites.dropna(inplace = True)

to_rem = []
for idx in test_reactions.index:
    tmp = test_reactions.loc[idx,'namesFormula'].split(" -> ")
    present = [x in test_metabolites.name.values for x in tmp]
    if not all(present):
        to_rem.append(idx)
    
test_reactions.drop(to_rem, inplace = True)

to_rem = []
for key,_ in enzymes.items():
    if not key in test_reactions.id.values:
        to_rem.append(key)

for rem in to_rem:
    enzymes.pop(rem,None)

test_enzyme = set()
for key,value in enzymes.items():
    for enz in value:
        test_enzyme.add(enz)

if os.path.isfile("./TESTinput/test_metabolites.csv"):
    os.remove("./TESTinput/test_metabolites.csv")
    test_metabolites.to_csv("./TESTinput/test_metabolites.csv", index = False)
else:
    test_metabolites.to_csv("./TESTinput/test_metabolites.csv", index = False)
    
if os.path.isfile("./TESTinput/test_reactions.csv"):
    os.remove("./TESTinput/test_reactions.csv")
    test_reactions.to_csv("./TESTinput/test_reactions.csv", index = False)
else:
    test_reactions.to_csv("./TESTinput/test_reactions.csv", index = False)
    
test_enzymes = pd.DataFrame(data = {"EC": list(test_enzyme)})
    
if os.path.isfile("./TESTinput/test_enzymes.csv"):
    os.remove("./TESTinput/test_enzymes.csv")
    test_enzymes.to_csv("./TESTinput/test_enzymes.csv", index = False)
else:
    test_enzymes.to_csv("./TESTinput/test_enzymes.csv", index = False)

test_enzyme_rxn_relation = pd.DataFrame(data = {'rxnID': [x for x in enzymes.keys()],
                                                'EC' : [' | '.join(x) for x in enzymes.values()]})
if os.path.isfile("./TESTinput/test_enzymes_reactions_relation.csv"):
    os.remove("./TESTinput/test_enzymes_reactions_relation.csv")
    test_enzyme_rxn_relation.to_csv("./TESTinput/test_enzymes_reactions_relation.csv", index = False)
else:
    test_enzyme_rxn_relation.to_csv("./TESTinput/test_enzymes_reactions_relation.csv", index = False)
    
test_mets_interest = test_metabolites.copy()

test_mets_interest.columns = ["name","ID","Inchi","smiles"]
if os.path.isfile("./TESTinput/test_molecules_interest.csv"):
    os.remove("./TESTinput/test_molecules_interest.csv")
    test_mets_interest.to_csv("./TESTinput/test_molecules_interest.csv", index = False)
else:
    test_mets_interest.to_csv("./TESTinput/test_molecules_interest.csv", index = False)

























