# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 13:19:20 2022

@author: fbalzerani
"""

from copy import deepcopy
import pickle

from .Common2 import Mutate, FindNeighbors, modifyAtomList, modifyBondList, ExtractCharge
from rdkit import Chem
from kcfconvoy import converter
import re
import itertools as it
from rdkit.Chem import AllChem, DataStructs

def GenerateProducts(substrateFilename, molecules_of_interest, operators, metabolites):
    # operators: filename, contains operator's information
    # inputList: input atom list
    # products: outputs after applying operators

    # construct a list from the input
    # first column is main atomtype and others are neighbors
    inputList, inputStructure, inputListNumbering, charge, mol_file = \
        ConstructInputListfromKEGG(substrateFilename, molecules_of_interest)
    
    if len(inputList) == 0:
        products = []
    else:
        # for each main and its neighbors, output the possible products
        products = ConstructOutputList(inputList, inputStructure, operators, mol_file, metabolites) 

    return inputList, inputStructure, inputListNumbering, products, charge

def ConstructInputListfromKEGG(substrateFilename, molecules_of_interest):
    inputStructure, charge, mol_file = ReadKCF(substrateFilename, molecules_of_interest)

    Mutate(inputStructure['Atoms'])

    if len(inputStructure['Bonds']) == 0:
        inputList = []
        inputListNumbering = []
    else:
        inputList, inputListNumbering = FindAllRMs(inputStructure)

    return inputList, inputStructure, inputListNumbering, charge, mol_file 

def ConstructOutputList(inputList, inputStructure, operators, query_mol_file, metabolites): 
    if isinstance(operators, str):
        f = open(operators, 'rb')
        data = pickle.load(f)
        operators = data['operators']
        operatorsMainAtom = data['operatorsMainAtom']
        operatorsMainAtomPos = data['operatorsMainAtomPos']
        f.close()
    else:
        # VP: Not tested
        assert(False)
        # this is for testing one operators
        operatorsMainAtom[0] = operators['Reactant']['R']
        operatorsMainAtomPos = 1

    if len(operators) == 0:
        # print('No operators!')
        return []

    # applying operators
    # we can use inputStructure to construct an input list from products
    return ApplyOperators(inputList, operators, operatorsMainAtom, operatorsMainAtomPos, query_mol_file, metabolites)    

def ApplyOperators(inputList, operators, operatorsMainAtom, operatorsMainAtomPos, query_mol_file, metabolites): 

    op = []
    
    #When there are more than one R, all of them should be able to be applied
    #otherwise the reaction is not applicable.
    
    for idx in range(len(operators)):
        f_query = AllChem.GetMorganFingerprint(query_mol_file,2) 
        mol_sub = Chem.MolFromSmiles(metabolites.loc[metabolites.name.isin(\
                                                [operators[idx]['KCF']['compound1']]),'smiles'].values[0])
        f_sub = AllChem.GetMorganFingerprint(mol_sub,2) 
        if not DataStructs.DiceSimilarity(f_query,f_sub) > 0.6: 
            continue 
        tmp_r = list(filter(lambda x: "R" in x, list(operators[idx]['KCF'].keys())))
        covered_Rs = {}
        input_final = {}
        input_BondToAdd = {}
        input_BondToRemove = {}
        input_BondToChange = {}
        for each_r in tmp_r:
            covered_Rs[each_r] = False
            input_final[each_r] = []
            obt_input = []
            for i, inRow in enumerate(inputList):
                
                number = int(re.findall("\d", each_r)[0])
                
                if operators[idx]['Reactant'][each_r] == inRow['R']:
                    tmp_inrow_m = list(filter(lambda x: re.findall("M\d$",x), 
                                              list(inRow.keys())))
                    tmp_op_m = list(filter(lambda x: re.findall("M"+str(number)+"_M\d$",x), 
                                           list(operators[idx]['Reactant'].keys())))
                    check = {}
                    for i_m in range(len(tmp_inrow_m)):
                        check[tmp_inrow_m[i_m]] = False
                        matchedPattern = False
                        if inRow[tmp_inrow_m[i_m]] == operators[idx]['Reactant'][tmp_op_m[i_m]]:
                            tmp_inrow_mn = list(filter(lambda x: re.findall(tmp_inrow_m[i_m]+\
                                                "\w$",x), list(inRow.keys())))
                            tmp_op_mn = list(filter(lambda x: re.findall(tmp_op_m[i_m] + \
                                        "\w$",x),list(operators[idx]['Reactant'].keys())))
                            mn_query = deepcopy(inRow[tmp_inrow_mn[0]])
                            mn_query.sort()
                            mn_op = deepcopy(operators[idx]['Reactant'][tmp_op_mn[0]])
                            mn_op.sort()
                            if mn_query == mn_op:
                                    matchedPattern = True
                        if matchedPattern:
                            check[tmp_inrow_m[i_m]] = True
                else:
                    continue
                
                if all([x == True for x in check.values()]):
                    covered_Rs[each_r] = True
                    obt_input.append(i+1)
            input_final[each_r] = obt_input.copy()
            
        """
        Add the check for eventual bonds added or removed between matched atoms, as well as
        possible bonds that are modified (passing from a double to a single for example: R00129)
        """
        if len(operators[idx]['KCF']['New_arrangement']\
               [operators[idx]['KCF']['compound1']]['ToAdd'])>0:
            count_NB = 0
            all_comb = {}
            for new_bond in operators[idx]['KCF']['New_arrangement']\
                                                    [operators[idx]['KCF']['compound1']]['ToAdd']:
                tmp_r = list(filter(lambda x: "R" in x, list(new_bond.keys())))
                covered_Rs['ToAdd%d' % count_NB] = False # check converage of new bond (NB)
                NB_Rs = {} # check coverage of each atom related to the NB
                input_BondToAdd['ToAdd%d' % count_NB] = {} # to add the input for each r of the NB
                for each_r in tmp_r:
                    NB_Rs[each_r] = False
                    input_BondToAdd['ToAdd%d' % count_NB][each_r] = []
                    
                    number = int(re.findall("\d", each_r)[0])
                    
                    atom_r = [x for x in operators[idx]['KCF']['atom1'] \
                                                                  if x[0] == new_bond[each_r]][0]
                    
                    tmp_op_m = list(filter(lambda x: re.findall("M"+str(number)+"$",x), 
                                            list(new_bond.keys())))
                    atom_m = [x[1] for x in operators[idx]['KCF']['atom1'] \
                                  for y in new_bond[tmp_op_m[0]] if x[0]== y]
                    atom_m.sort()
                    
                    tmp_op_mn = list(filter(lambda x: re.findall(tmp_op_m[0] + \
                                "\w$",x),list(new_bond.keys())))
                    
                    atom_mn = [x[1] for x in operators[idx]['KCF']['atom1'] \
                                      for y in new_bond[tmp_op_mn[0]] \
                                          for z in y if x[0]== z]
                    atom_mn.sort()
                    obt_input = []
                    for i, inRow in enumerate(inputList):
                        
                        if atom_r[1] == inRow['R']:
                            tmp_inrow_m = list(filter(lambda x: re.findall("M\d$",x), 
                                                      list(inRow.keys())))
                            
                            all_inrow_m = [inRow[x] for x in tmp_inrow_m if not inRow[x] == []]
                            all_inrow_m.sort()
                            if atom_m == all_inrow_m:
                                tmp_inrow_mn = list(filter(lambda x: re.findall("M\d\w$",x), list(inRow.keys())))
                                all_inrow_mn = [inRow[x] for x in tmp_inrow_mn if not inRow[x] == []]
                                all_inrow_mn = [y for x in all_inrow_mn for y in x]
                                all_inrow_mn.sort()
                                if atom_mn == all_inrow_mn:            
                                    NB_Rs[each_r] = True
                                    obt_input.append(i+1)                                
                        else:
                            continue
                            
                    input_BondToAdd['ToAdd%d' % count_NB][each_r] = obt_input.copy()
                if all([x == True for x in NB_Rs.values()]):
                    covered_Rs['ToAdd%d' % count_NB] = True
                    keys_toAdd = input_BondToAdd['ToAdd%d' % count_NB].keys()
                    combinations_toAdd = it.product(*(input_BondToAdd['ToAdd%d' % count_NB][Name] for Name in keys_toAdd))
                    all_comb[count_NB] = []
                    for comb in combinations_toAdd:
                        tmp = list(comb)
                        tmp.append(new_bond['order'])
                        all_comb[count_NB].append(tmp)
                count_NB +=1
            input_toAdd = [list(x) for x in it.product(*(all_comb[Name] for Name in all_comb.keys()))]
        else:
            input_toAdd = []
            
        if len(operators[idx]['KCF']['New_arrangement']\
               [operators[idx]['KCF']['compound1']]['ToRemove'])>0:
            count_NB = 0
            all_comb = {}
            for new_bond in operators[idx]['KCF']['New_arrangement']\
                                                    [operators[idx]['KCF']['compound1']]['ToRemove']:
                tmp_r = list(filter(lambda x: "R" in x, list(new_bond.keys())))
                covered_Rs['ToRemove%d' % count_NB] = False # check converage of new bond (NB)
                NB_Rs = {} # check coverage of each atom related to the NB
                input_BondToRemove['ToRemove%d' % count_NB] = {} # to add the input for each r of the NB
                for each_r in tmp_r:
                    NB_Rs[each_r] = False
                    input_BondToRemove['ToRemove%d' % count_NB][each_r] = []
                    
                    number = int(re.findall("\d", each_r)[0])
                    
                    atom_r = [x for x in operators[idx]['KCF']['atom1'] \
                                                                  if x[0] == new_bond[each_r]][0]
                    
                    tmp_op_m = list(filter(lambda x: re.findall("M"+str(number)+"$",x), 
                                            list(new_bond.keys())))
                    atom_m = [x[1] for x in operators[idx]['KCF']['atom1'] \
                                  for y in new_bond[tmp_op_m[0]] if x[0]== y]
                    atom_m.sort()
                    
                    tmp_op_mn = list(filter(lambda x: re.findall(tmp_op_m[0] + \
                                "\w$",x),list(new_bond.keys())))
                    
                    atom_mn = [x[1] for x in operators[idx]['KCF']['atom1'] \
                                      for y in new_bond[tmp_op_mn[0]] \
                                          for z in y if x[0]== z]
                    atom_mn.sort()
                    obt_input = []
                    for i, inRow in enumerate(inputList):
                        
                        if atom_r[1] == inRow['R']:
                            tmp_inrow_m = list(filter(lambda x: re.findall("M\d$",x), 
                                                      list(inRow.keys())))
                            
                            all_inrow_m = [inRow[x] for x in tmp_inrow_m if not inRow[x] == []]
                            all_inrow_m.sort()
                            if atom_m == all_inrow_m:
                                tmp_inrow_mn = list(filter(lambda x: re.findall("M\d\w$",x), list(inRow.keys())))
                                all_inrow_mn = [inRow[x] for x in tmp_inrow_mn if not inRow[x] == []]
                                all_inrow_mn = [y for x in all_inrow_mn for y in x]
                                all_inrow_mn.sort()
                                if atom_mn == all_inrow_mn:            
                                    NB_Rs[each_r] = True
                                    obt_input.append(i+1)                                
                        else:
                            continue
                            
                    input_BondToRemove['ToRemove%d' % count_NB][each_r] = obt_input.copy()
                if all([x == True for x in NB_Rs.values()]):
                    covered_Rs['ToRemove%d' % count_NB] = True
                    keys_toRemove = input_BondToRemove['ToRemove%d' % count_NB].keys()
                    combinations_toRemove = it.product(*(input_BondToRemove['ToRemove%d' % count_NB][Name] for Name in keys_toRemove))
                    all_comb[count_NB] = []
                    for comb in combinations_toRemove:
                        all_comb[count_NB].append(comb)
                count_NB +=1
            input_toRemove = [list(x) for x in it.product(*(all_comb[Name] for Name in all_comb.keys()))]
        else:
            input_toRemove = []
            
        if len(operators[idx]['KCF']['New_arrangement']\
               [operators[idx]['KCF']['compound1']]['ToChange'])>0:
            count_NB = 0
            all_comb = {}
            for new_bond in operators[idx]['KCF']['New_arrangement']\
                                                    [operators[idx]['KCF']['compound1']]['ToChange']:
                tmp_r = list(filter(lambda x: "R" in x, list(new_bond.keys())))
                covered_Rs['ToChange%d' % count_NB] = False # check converage of new bond (NB)
                NB_Rs = {} # check coverage of each atom related to the NB
                input_BondToChange['ToChange%d' % count_NB] = {} # to add the input for each r of the NB
                for each_r in tmp_r:
                    NB_Rs[each_r] = False
                    input_BondToChange['ToChange%d' % count_NB][each_r] = []
                    
                    number = int(re.findall("\d", each_r)[0])
                    
                    atom_r = [x for x in operators[idx]['KCF']['atom1'] \
                                                                  if x[0] == new_bond[each_r]][0]
                    
                    tmp_op_m = list(filter(lambda x: re.findall("M"+str(number)+"$",x), 
                                            list(new_bond.keys())))
                    atom_m = [x[1] for x in operators[idx]['KCF']['atom1'] \
                                  for y in new_bond[tmp_op_m[0]] if x[0]== y]
                    atom_m.sort()
                    
                    tmp_op_mn = list(filter(lambda x: re.findall(tmp_op_m[0] + \
                                "\w$",x),list(new_bond.keys())))
                    
                    atom_mn = [x[1] for x in operators[idx]['KCF']['atom1'] \
                                      for y in new_bond[tmp_op_mn[0]] \
                                          for z in y if x[0]== z]
                    atom_mn.sort()
                    obt_input = []
                    for i, inRow in enumerate(inputList):
                        
                        if atom_r[1] == inRow['R']:
                            tmp_inrow_m = list(filter(lambda x: re.findall("M\d$",x), 
                                                      list(inRow.keys())))
                            
                            all_inrow_m = [inRow[x] for x in tmp_inrow_m if not inRow[x] == []]
                            all_inrow_m.sort()
                            if atom_m == all_inrow_m:
                                tmp_inrow_mn = list(filter(lambda x: re.findall("M\d\w$",x), list(inRow.keys())))
                                all_inrow_mn = [inRow[x] for x in tmp_inrow_mn if not inRow[x] == []]
                                all_inrow_mn = [y for x in all_inrow_mn for y in x]
                                all_inrow_mn.sort()
                                if atom_mn == all_inrow_mn:            
                                    NB_Rs[each_r] = True
                                    obt_input.append(i+1)                                
                        else:
                            continue
                            
                    input_BondToChange['ToChange%d' % count_NB][each_r] = obt_input.copy()
                if all([x == True for x in NB_Rs.values()]):
                    covered_Rs['ToChange%d' % count_NB] = True
                    keys_ToChange = input_BondToChange['ToChange%d' % count_NB].keys()
                    combinations_toAdd = it.product(*(input_BondToChange['ToChange%d' % count_NB][Name] for Name in keys_ToChange))
                    all_comb[count_NB] = []
                    for comb in combinations_toAdd:
                        tmp = list(comb)
                        tmp.append(new_bond['order'])
                        all_comb[count_NB].append(tmp)
                count_NB +=1
            input_toChange = [list(x) for x in it.product(*(all_comb[Name] for Name in all_comb.keys()))]
        else:
            input_toChange = []
        """
        End check for eventual bonds added or removed between matched atoms, as well as
        possible bonds that are modified (passing from a double to a single for example: R00129)
        """
        
        if all([x == True for x in covered_Rs.values()]):
            keys = input_final.keys()
            combinations = it.product(*(input_final[Name] for Name in keys))
            
            if len(input_toAdd) == 0 and len(input_toRemove) == 0 and len(input_toChange) == 0:
                for each_comb in combinations:
                    if len(set(each_comb)) == len(list(each_comb)):
                        entry = operators[idx].copy()
                        entry['ObtainedFromInput'] = list(each_comb)
                        entry['BondToAdd'] = []
                        entry['BondToRemove'] = []
                        entry['BondToChange'] = []
                        op.append(entry)
            elif not len(input_toAdd) == 0 and len(input_toRemove) == 0 and len(input_toChange) == 0:
                for each_comb in combinations:
                    if len(set(each_comb)) == len(list(each_comb)):
                        for each_toAdd in input_toAdd:
                            entry = operators[idx].copy()
                            entry['ObtainedFromInput'] = list(each_comb)
                            entry['BondToAdd'] = each_toAdd
                            entry['BondToRemove'] = []
                            entry['BondToChange'] = []
                            op.append(entry)
            elif not len(input_toAdd) == 0 and len(input_toRemove) == 0 and not len(input_toChange) == 0:
                for each_comb in combinations:
                    if len(set(each_comb)) == len(list(each_comb)):
                        for each_toAdd in input_toAdd:
                            for each_toChange in input_toChange:
                                entry = operators[idx].copy()
                                entry['ObtainedFromInput'] = list(each_comb)
                                entry['BondToAdd'] = each_toAdd
                                entry['BondToRemove'] = []
                                entry['BondToChange'] = each_toChange
                                op.append(entry)
            elif len(input_toAdd) == 0 and not len(input_toRemove) == 0 and len(input_toChange) == 0:
                for each_comb in combinations:
                    if len(set(each_comb)) == len(list(each_comb)):
                        for each_toRemove in input_toRemove:
                            entry = operators[idx].copy()
                            entry['ObtainedFromInput'] = list(each_comb)
                            entry['BondToAdd'] = []
                            entry['BondToRemove'] = each_toRemove
                            entry['BondToChange'] = []
                            op.append(entry)
            elif len(input_toAdd) == 0 and not len(input_toRemove) == 0 and not len(input_toChange) == 0:
                for each_comb in combinations:
                    if len(set(each_comb)) == len(list(each_comb)):
                        for each_toRemove in input_toRemove:
                            for each_toChange in input_toChange:
                                entry = operators[idx].copy()
                                entry['ObtainedFromInput'] = list(each_comb)
                                entry['BondToAdd'] = []
                                entry['BondToRemove'] = each_toRemove
                                entry['BondToChange'] = each_toChange
                                op.append(entry)
            elif not len(input_toAdd) == 0 and not len(input_toRemove) == 0 and len(input_toChange) == 0:
                for each_comb in combinations:
                    if len(set(each_comb)) == len(list(each_comb)):
                        for each_toAdd in input_toAdd:
                            for each_toRemove in input_toRemove:
                                entry = operators[idx].copy()
                                entry['ObtainedFromInput'] = list(each_comb)
                                entry['BondToAdd'] = each_toAdd
                                entry['BondToRemove'] = each_toRemove
                                entry['BondToChange'] = []
                                op.append(entry)
            elif not len(input_toAdd) == 0 and not len(input_toRemove) == 0 and not len(input_toChange) == 0:
                for each_comb in combinations:
                    if len(set(each_comb)) == len(list(each_comb)):
                        for each_toAdd in input_toAdd:
                            for each_toRemove in input_toRemove:
                                for each_toChange in input_toChange:
                                    entry = operators[idx].copy()
                                    entry['ObtainedFromInput'] = list(each_comb)
                                    entry['BondToAdd'] = each_toAdd
                                    entry['BondToRemove'] = each_toRemove
                                    entry['BondToChange'] = each_toChange
                                    op.append(entry)
            elif len(input_toAdd) == 0 and len(input_toRemove) == 0 and not len(input_toChange) == 0:
                for each_comb in combinations:
                    if len(set(each_comb)) == len(list(each_comb)):
                        for each_toChange in input_toChange:
                            entry = operators[idx].copy()
                            entry['ObtainedFromInput'] = list(each_comb)
                            entry['BondToAdd'] = []
                            entry['BondToRemove'] = []
                            entry['BondToChange'] = each_toChange
                            op.append(entry)
    return op
    
def ReadKCF(kcfFilename, molecules_of_interest):
    smiles = molecules_of_interest.loc[molecules_of_interest['name'].isin([kcfFilename]),
                                      'smiles'].values[0]
    mol_file = Chem.MolFromSmiles(smiles)
    
    # Extract charge
    # Extract from the M  CHG field the charges. the 3rd field represents the number of
    # charged atoms. Then there are paired values (1st the position, 2nd the value of the charge)
    # Then extract all the pairs and store in the charge variable
    charge = ExtractCharge(mol_file)
    
    #Define the KCF and extract atoms and bonds
    kcf_file = converter.rdkmol_to_kcf(mol_file)
    
    #Extract the atoms and bonds information of the molecule
    atom_file = re.findall("ATOM\s+\d+\s+([\S\s]+)BOND",kcf_file)[0].strip().split("\n")
    atom_file = [x.strip() for x in atom_file]
    
    a_molblock = Chem.MolToMolBlock(mol_file)
    tmp_a = re.split("\n", a_molblock)
    min_pos = min([x for x in range(len(tmp_a)) if len(re.split("\s+",tmp_a[x].strip())) == 4])
    max_pos = min([x for x in range(len(tmp_a)) if re.findall("M  ",tmp_a[x])])
    tmp_a = tmp_a[min_pos:max_pos]
    
    bond_file = []
    count = 0
    for x in tmp_a:
        count += 1
        tmp = re.split("\s+",x.strip())[:-1]
        if len(tmp[0]) == 4:
            cor1 = tmp[0][:1]
            cor2 = tmp[0][1:]
        elif len(tmp[0]) == 5:
            cor1 = tmp[0][:2]
            cor2 = tmp[0][2:]
        elif len(tmp[0]) == 6:
            cor1 = tmp[0][:3]
            cor2 = tmp[0][3:]
        else:
            cor1 = tmp[0]
            cor2 = tmp[1]
            
        tmp_2 = [str(count), cor1, cor2, tmp[-1]]
        bond_file.append(' '.join(tmp_2))
    
    atoms = modifyAtomList(atom_file)
    bonds = modifyBondList(bond_file)
    molecule = {'Atoms': atoms, 'Bonds': bonds}
    return molecule, charge, mol_file

def FindAllRMs(inputStructure):
    # For each atom, finds the neighbors and put them in a row of inputList

    # Bond: 1st&2nd cols-atom#, 3rd number of bonds
    inputList = []
    for atomIndex, atomType, element, x, y in inputStructure['Atoms']:
        # first col of inputList is main atom
        entry = {'R': atomType}
        try:
            nghrs, tempVar = FindNeighbors(inputStructure['Bonds'], atomIndex)
        except ValueError:
            return [],[]
        
        neighbors = [inputStructure['Atoms'][x - 1][1] for x in nghrs]
        for j in range(len(neighbors)):
            entry['M%d' % (j + 1)] = neighbors[j]
            nghrsN, tempVar = FindNeighbors(inputStructure['Bonds'], nghrs[j])
            neighborsN = [inputStructure['Atoms'][x - 1][1] for x in nghrsN]
            neighborsN.sort()
            entry['M%dN' % (j + 1)] = neighborsN
        for j in range(len(neighbors) + 1, 5):
            entry['M%d' % j] = []
            entry['M%dN' % j] = []
        inputList.append(entry)

    inputList, inputListNumbering = SortFields(inputList, inputStructure)

    if False:
        print('\t'.join(inputList[0].keys()))
        for row in inputList:
            print('\t'.join(str(x) for x in row.values()))
        print()
        for row in inputListNumbering:
            print('\t'.join(str(x) for x in row))
        exit()

    return inputList, inputListNumbering

def SortFields(inputList, inputStructure):
    inputListNumbering = []
    outputList = deepcopy(inputList)
    for i, inRow in enumerate(inputList):
        tempData = []
        for j in range(1, 5):
            if len(inRow['M%d' % j]) == 0 or len(inRow['M%dN' % j]) == 0:
                break
            tempData.append(inRow['M%d' % j] + ''.join(inRow['M%dN' % j]))
        neighboringAtom, neighborsPos = FindNeighbors(inputStructure['Bonds'], i + 1)
        for j, pos in enumerate(neighborsPos):
            tempData[j] += str(inputStructure['Bonds'][pos - 1][3])
        if len(tempData) == 0:
            continue

        sortedMs = list(enumerate(tempData))
        sortedMs.sort(key=lambda x: x[1])
        order = [x[0] for x in sortedMs]

        numRow = [i + 1, 0, 0, 0, 0]
        for j, x in enumerate(order):
            numRow[j + 1] = neighboringAtom[x]
        inputListNumbering.append(numRow)

        for j, k in enumerate(order):
            outputList[i]['M%d' % (j + 1)] = inRow['M%d' % (k + 1)]
            outputList[i]['M%dN' % (j + 1)] = inRow['M%dN' % (k + 1)]
    return outputList, inputListNumbering
