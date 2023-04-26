# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 11:08:18 2022

@author: fbalzerani
"""

from copy import deepcopy
from math import sqrt, cos, sin
from math import pi as m_pi
import re

from .Common2 import FindNeighbors
from rdkit import Chem
Chem.WrapLogs()
import json
import math
# from io import StringIO
# import sys

def GenerateProductStructure(inputStructure, inputListNumbering, product, charge, 
                             charge_Template):
    tmp_m = list(filter(lambda x: "M" in x, list(product['KCF'].keys())))
    if all([len(product['KCF'][each_m]) == 0 for each_m in tmp_m]):
        return None
    
    # adding the added structure to inputStructure
    return BuildNewStructure(product, inputStructure,
                             [inputListNumbering[x - 1] for x in product['ObtainedFromInput']], 
                             charge, charge_Template)

def GenerateMolFiles(targetMoleculeDir, inputStructure, inputListNumbering, products, charge_orig,
                     charge_Template_orig, query_info, makeKcf=False):
    inputListNumbering_orig = deepcopy(inputListNumbering) 
    for i, product in enumerate(products):
        
        charge_Template = deepcopy(charge_Template_orig)
        charge = deepcopy(charge_orig)
        
        inputListNumbering = deepcopy(inputListNumbering_orig) 
        
        modifiedInputStructure, charge, charge_Template = \
            GenerateProductStructure(inputStructure, inputListNumbering,
                                                          product, charge, charge_Template)
        if modifiedInputStructure is None:
            continue

        filename = 'product_%d' % (i + 1)

        if not makeKcf:
            ConstructMolFile(modifiedInputStructure, charge, 
                             charge_Template, targetMoleculeDir + '/' + filename + '.json', product, query_info)
        else: # VP: added KCF export
            ConstructKcfFile(modifiedInputStructure, 
                             targetMoleculeDir + '/' + filename + '.kcf')
            
def BuildNewStructure(pattern, inputStructure, inputNumbering, charge, charge_Template):
    # pattern: pattern applied to R atom
    # inputStructure: atom and bond information for the test molecule
    # inputNumbering: array with first element is R, others are Ms
    err = 1
    # print(pattern['ObtainedFromInput'])
    # print([x[0] for x in inputNumbering])
    modifiedInputStructure = deepcopy(inputStructure)
    charge_Template_orig = deepcopy(charge_Template)
    
    origLength = len(inputStructure['Atoms'])
    # changing the atom type for R
    
    tmp_r = list(filter(lambda x: "R" in x, list(pattern['KCF'].keys())))
    
    # Extract the atom part of the alignment to not add them again
    tmp_atom = set()
    for each_r in tmp_r:
        number = int(re.findall("\d", each_r)[0])
        tmp_al = pattern['KCF']['Align%d' % number]
        for each in tmp_al:
            if not each[0] == -1 and not each[1] == -1:
                tmp_atom.add(each[1])                
    
    already_insert = []
    for each in tmp_atom:
        already_insert.append(pattern['KCF']['atom2'][each -1])
    
    # Define the variable to keep on the atom mapping through different Rs and
    # the variable to take note of the already added bonds to not add them again
    generalMapping = []
    used_bonds = []
    for bond in pattern['KCF']['bond2']:
        if bond[1] in tmp_atom and bond[2] in tmp_atom:
            used_bonds.append(bond)
            
    for idx in range(len(tmp_r)):
        
        charge_Template = deepcopy(charge_Template_orig)
        number = int(re.findall("\d", tmp_r[idx])[0])
        
        ChangeR(modifiedInputStructure, pattern, inputNumbering[idx], number)
        if any(0 in bond for bond in modifiedInputStructure['Bonds']):
            err = 0
    
        # changing Ms
        modifiedInputStructure, charge, inputNumbering[idx], atomTobeRemoved, bond_removed = \
                                                                    ChangeMs(modifiedInputStructure, 
                                                pattern, inputNumbering[idx], charge, number)
        
        if any([b[1] in [x[0] for x in inputNumbering] and \
                b[2] in [x[0] for x in inputNumbering] for b in bond_removed]):
            return None, [], []
        if any([r in atomTobeRemoved for r in [x[0] for x in inputNumbering]]):
            return None, [], []
        # block the process if we removed atoms added along the process
        if any([r > origLength for r in atomTobeRemoved]):
            return None, [], []
        if any(0 in bond for bond in modifiedInputStructure['Bonds']):
            err = 0
                        
        # changing D
        ChangeD(modifiedInputStructure, pattern, inputNumbering[idx], inputStructure, 
                    charge_Template, number, already_insert, generalMapping, used_bonds)
        if any(0 in bond for bond in modifiedInputStructure['Bonds']):
            err = 0
    
        if err == 0 or len(modifiedInputStructure['Bonds']) == 0:
            return modifiedInputStructure, charge, charge_Template # comment C3

        # removing atoms not connected to the structure
        modifiedInputStructure, charge, charge_Template, atomTobeRemoved = RemovePartsOfStructure(modifiedInputStructure, 
                                                inputNumbering[idx][0], charge, charge_Template)
        
        if any([r in atomTobeRemoved for r in [x[0] for x in inputNumbering]]):
            return None, [], []
        # block the process if we removed atoms added along the process
        if any([r > origLength for r in atomTobeRemoved]):
            return None, [], []
    
    if False:
        print('Bonds\n' + '\n'.join(str(x) for x in modifiedInputStructure['Bonds']))
        print('Atoms\n' + '\n'.join(str(x) for x in modifiedInputStructure['Atoms']))
    
    # Add bonds for new arrangement
    if len(pattern["BondToAdd"]) > 0:
        for bond in pattern["BondToAdd"]:
            if bond[0] in [x[0] for x in modifiedInputStructure['Atoms']] and \
                bond[1] in [x[0] for x in modifiedInputStructure['Atoms']]:
                    if not CheckBondPresence(bond, modifiedInputStructure):
                        b0 = max([x[0] for x in modifiedInputStructure['Atoms']]) + 1
                        b1 = bond[0]
                        b2 = bond[1]
                        b3 = bond[2]
                        modifiedInputStructure['Bonds'].append([b0,b1,b2,b3])
            else:
                return None,[],[]
    # Remove bonds for new arrangement
    if len(pattern["BondToRemove"]) > 0:
        to_rem = []
        for bond in pattern["BondToRemove"]:
            if bond[0] in [x[0] for x in modifiedInputStructure['Atoms']] and \
                bond[1] in [x[0] for x in modifiedInputStructure['Atoms']]:
                for bond_present in modifiedInputStructure['Bonds']:
                    if bond[0] == bond_present[1] and bond[1] == bond_present[2]:
                        if not bond_present in to_rem:
                            to_rem.append(bond_present)
                            break
                    if bond[1] == bond_present[1] and bond[0] == bond_present[2]:
                        if not bond_present in to_rem:
                            to_rem.append(bond_present)
                            break
            else:
                return None,[],[]
        for each_b in to_rem:
            modifiedInputStructure['Bonds'].pop(modifiedInputStructure['Bonds'].index(each_b))
            
    # Modify bonds for new arrangement
    if len(pattern["BondToChange"]) > 0:
        for bond in pattern["BondToChange"]:
            if bond[0] in [x[0] for x in modifiedInputStructure['Atoms']] and \
                bond[1] in [x[0] for x in modifiedInputStructure['Atoms']]:
                for bond_present in modifiedInputStructure['Bonds']:
                    if bond[0] == bond_present[1] and bond[1] == bond_present[2]:
                        bond_present[3] = bond[2]
                        break
                    if bond[1] == bond_present[1] and bond[0] == bond_present[2]:
                        bond_present[3] = bond[2]
                        break
            else:
                return None,[],[]
    
    to_rem = []
    for key,value in charge.items():
        n_bonds_input = 0
        n_bonds_modified = 0
        for b in inputStructure['Bonds']:
            if int(key) == b[1] or int(key) == b[2]:
                n_bonds_input += 1
        for b in modifiedInputStructure['Bonds']:
            if int(key) == b[1] or int(key) == b[2]:
                n_bonds_modified += 1
        
        if n_bonds_input == n_bonds_modified:
            continue
        else:
            new_charge = int(value)-(n_bonds_input-n_bonds_modified)
            if new_charge == 0:
                to_rem.append(key)
            else:
                charge[key] = str(new_charge)
    for r in to_rem:
        charge.pop(r)
        
    try:
        modifiedInputStructure, charge, charge_Template = UpdateNumbering(modifiedInputStructure, charge, 
                                                                      charge_Template)
    except KeyError:
        return None,[],[]
    
    return modifiedInputStructure, charge, charge_Template

def ChangeR(oStruct, pattern, RMList, number_R):
    # change the atom type for R
    # the tmp variable is necessary cause of the not renumbering of the RMList
    # nor Atoms in oStruct. Therefore, the RMList[0] do not correspond to the proper
    # position in oStruct['Atoms']
    tmp = [x for x in range(len(oStruct['Atoms'])) if oStruct['Atoms'][x][0] == RMList[0]][0]
    oStruct['Atoms'][tmp][1] = pattern['Product']['R%d' % number_R]

def ChangeMs(oStruct, pattern, RMList, charge, number_R):
    # removing the bond between a removed neighbor and R
    rmListLocal = deepcopy(RMList)
    atomTobeRemoved = []
    bond_removed = []
    for z in range(1, len(rmListLocal)): # start from "1" to skip the R at index "0"
        if rmListLocal[z] != 0:
            if pattern['Product']["M"+str(number_R)+'_M%d' % z] == '*': # note: z doesn't 
                                                                    # need offset of 1 here
                indexRM1 = FindBondBetween2atoms(oStruct['Bonds'], rmListLocal[0], 
                                                 rmListLocal[z])
                if indexRM1 is not None:
                    # the tmp variable is necessary cause of the not renumbering of the bonds
                    # Therefore, the indexRM1 do not correspond to the proper
                    # position in oStruct['Bonds']
                    tmp = [x for x in range(len(oStruct['Bonds'])) if oStruct['Bonds'][x][0] == indexRM1][0]
                    bond_removed.append(oStruct['Bonds'][tmp])
                    del oStruct['Bonds'][tmp]
                    UpdateBondsM(oStruct, len(oStruct['Bonds']), indexRM1)
                Mneighbor1, _ = FindNeighbors(oStruct['Bonds'], rmListLocal[z])
                if len(Mneighbor1) == 0:
                    oStruct, charge, _, tmp_atomTobeRemoved = RemoveAtom(oStruct, [rmListLocal[z]], charge, 
                                                 {})
                    for each_a in tmp_atomTobeRemoved:
                        atomTobeRemoved.append(each_a)
                    
    return oStruct, charge, rmListLocal, atomTobeRemoved, bond_removed

def UpdateBondsM(oStruct, nBonds, indexBond):
    for m in range(indexBond, nBonds + 1):
        oStruct['Bonds'][m - 1][0] = m

def RemoveAtom(oStruct, atomTobeRemoved, charge, charge_Template):
    if len(atomTobeRemoved) == 0:
        return oStruct, charge, charge_Template, atomTobeRemoved
    # remove atoms
    oStruct['Atoms'] = [atom for atom in oStruct['Atoms'] if atom[0] not in atomTobeRemoved]
    # Remove the information about the charge if needed
    to_del = []
    for key,_ in charge.items():
        if int(key) in atomTobeRemoved:
            to_del.append(key)
            
    for value in to_del:
        charge.pop(value)
    # update atom numbers
    for m in atomTobeRemoved:
        indices = FindNeighborsP(oStruct['Bonds'], m)
        oStruct['Bonds'] = [bond for bond in oStruct['Bonds'] if bond[0] not in indices]
    
    return oStruct, charge, charge_Template, atomTobeRemoved

def ChangeD(oStruct, pattern, RMList, iStruct, charge_Template, number_R, already_insert, 
            generalMapping, used_bonds):
    AddNewFunctionalGroups(oStruct, pattern, RMList, iStruct, charge_Template, number_R,
                           already_insert, generalMapping, used_bonds)
    
def RemovePartsOfStructure(oStruct, nodeR, charge, charge_Template):
    visitedAtoms = set()
    atomsToVisit = [nodeR]
    while len(atomsToVisit) > 0:
        atom = atomsToVisit.pop()
        visitedAtoms.add(atom)
        for index, atom1, atom2, order in oStruct['Bonds']:
            if atom1 == atom and atom2 not in visitedAtoms:
                atomsToVisit.append(atom2)
            if atom2 == atom and atom1 not in visitedAtoms:
                atomsToVisit.append(atom1)
    unvisitedAtoms = set([x[0] for x in oStruct['Atoms']]) - visitedAtoms
    oStruct, charge, charge_Template, atomTobeRemoved = RemoveAtom(oStruct, unvisitedAtoms, charge, 
                                                               charge_Template)
    return oStruct, charge, charge_Template, atomTobeRemoved

def AddNewFunctionalGroups(oStruct, pattern, RMList, iStruct, charge_Template, number_R, 
                           already_insert, generalMapping, used_bonds):
    number_Ds = list(filter(lambda x: re.findall("D"+str(number_R)+"_D\d",x), 
                                pattern['Product'].keys()))
    # Don't repeat addition of bonds
    pos_d = []
    for d in range(1,len(number_Ds)+1):
        try:
            pos_d.append(pattern['Product']["D"+str(number_R)+"_D%d" % d]['atomD'])   
        except KeyError:
            continue
    
    general_atoms = []
    for n_D in number_Ds:
        if len(pattern['Product'][n_D]['D']) != 0 and\
                len(pattern['Product'][n_D]['atom']) != 0:
            A, B, C, D = [(0, 0)] * 4
            O = (iStruct['Atoms'][RMList[0] - 1][3], iStruct['Atoms'][RMList[0] - 1][4]) # (x, y)
            # finds O2x position
            if pattern['Product'][n_D]['D'] == 'O2x':
                if 'O2x' in pattern['Product']["M"+str(number_R)+'_M1N'] or RMList[2] == 0:
                    A = (iStruct['Atoms'][RMList[1] - 1][3], iStruct['Atoms'][RMList[1] - 1][4])
                else:
                    A = (iStruct['Atoms'][RMList[2] - 1][3], iStruct['Atoms'][RMList[2] - 1][4])
            else:
                A = (iStruct['Atoms'][RMList[1] - 1][3], iStruct['Atoms'][RMList[1] - 1][4])
                if RMList[2] != 0:
                    B = (iStruct['Atoms'][RMList[2] - 1][3], iStruct['Atoms'][RMList[2] - 1][4])
                if RMList[3] != 0:
                    C = (iStruct['Atoms'][RMList[3] - 1][3], iStruct['Atoms'][RMList[3] - 1][4])
                if RMList[4] != 0:
                    D = (iStruct['Atoms'][RMList[4] - 1][3], iStruct['Atoms'][RMList[4] - 1][4])
    
            lenX = 0.8
            # Finding the position of D
            X = getThirdVector(O, lenX, A, B, C, D)
            new_x = []
            for x in X:
                if math.isnan(x):
                    new_x.append(0.0)
                else:
                    new_x.append(x)
            X = tuple(new_x)
           
            atomList = []
            for index, type, element, x, y in pattern['Product'][n_D]['atom']:
                atomList.append(index)
            atomD = atomList.index(pattern['Product'][n_D]['atomD'])
    
            # adding new atoms (D)
            atomMapping = [(pattern['Product'][n_D]['atomR'], RMList[0])]
            # adding also in the general mapping variable
            if not pattern['Product'][n_D]['atomR'] in [x1 for x1, x2 in generalMapping]:
                generalMapping.append((pattern['Product'][n_D]['atomR'], RMList[0]))
                
            for atom in pattern['Product'][n_D]['atom']:
                first_general_atoms = [x[0] for x in general_atoms]
                if atom[0] in first_general_atoms:
                    pos_general_atom = first_general_atoms.index(atom[0])
                    atomMapping.append(general_atoms[pos_general_atom])
                    generalMapping.append(general_atoms[pos_general_atom])
                elif atom in already_insert: # if already added the atom, not add again
                    if atom[0] in [x[0] for x in generalMapping]:
                        continue
                    tmp_pos_atom = [x[0] for x in pattern['KCF']['Align1'] if x[1] == atom[0]]
                    generalMapping.append((atom[0],tmp_pos_atom[0]))
                else:
                    index = max([x[0] for x in oStruct['Atoms']]) + 1
                    general_atoms.append((atom[0], index))
                    atomMapping.append((atom[0], index))
                    oStruct['Atoms'].append([index, atom[1], atom[2], # index, type, element
                        X[0] + (atom[3] - pattern['Product'][n_D]['atom'][atomD][3]), # x
                        X[1] + (atom[4] - pattern['Product'][n_D]['atom'][atomD][4])]) # y
                    already_insert.append(atom)
                    generalMapping.append((atom[0], index))
            
            # When there are multiple Rs, can happen to add an atom into another loop
            # even though they are the D of the R. Therefore it would be excluded from
            # the new atomMapping, so add it manually otherwise will lose the addition of the
            # bond.
            if not pattern['Product'][n_D]['atomD'] in [x[0] for x in atomMapping]:
                atomMapping.append((pattern['Product'][n_D]['atomD'],
                   [x[1] for x in generalMapping if pattern['Product'][n_D]['atomD'] == x[0]][0]))
            
            # adding the bond between R and D
            RDRowinBond = FindBondBetween2atoms(pattern['Product'][n_D]['bond'],
                                                pattern['Product'][n_D]['atomR'],
                                                pattern['Product'][n_D]['atomD'])
            
            # Add the bond between R and D
            for atom1,atom2 in atomMapping:
                if atom1 == pattern['Product'][n_D]['atomD']:
                    tmp_bond_check = pattern['Product'][n_D]['bond'][RDRowinBond - 1][1:]
                    tmp_used_bonds = [[x[1],x[2],x[3]] for x in used_bonds]
                    if not tmp_bond_check in tmp_used_bonds:
                        if len(oStruct['Bonds']) == 0:
                            index_bonds = 1
                        else:
                            index_bonds = max([x[0] for x in oStruct['Bonds']]) + 1
                        oStruct['Bonds'].append([index_bonds, 
                                     RMList[0], # index, atom1, atom2 (below)
                                     atom2, pattern['Product'][n_D]['bond'][RDRowinBond - 1][3] # order
                        ])
                        
                        used_bonds.append(pattern['Product'][n_D]['bond'][RDRowinBond - 1])
                else:
                    continue
    
            # adding other bonds in the functional group D
            for bond in pattern['Product'][n_D]['bond']:
                # print(bond)
                index, atom1, atom2, order = bond
                bRemoveBond = False
                if (atom1 in pos_d and not atom1 == pattern['Product'][n_D]['atomD']) \
                        or (atom2 in pos_d and not atom2 == pattern['Product'][n_D]['atomD']):
                    continue
                # skip bond if they are already added (a bond can be define as forward or reverse
                # so check both senses)
                if atom1 in [x[1] for x in used_bonds] and atom2 in [x[2] for x in used_bonds] \
                    and order in [x[3] for x in used_bonds if atom1 == x[1] and atom2 == x[2]]:
                    continue
                if atom1 in [x[2] for x in used_bonds] and atom2 in [x[1] for x in used_bonds] \
                    and order in [x[3] for x in used_bonds if atom1 == x[2] and atom2 == x[1]]:
                    continue
                if index != RDRowinBond:
                    bRemoveBond = False
                    newBond = [max([x[0] for x in oStruct['Bonds']]) + 1, 0, 0, order]
                    for j in [1, 2]: # atoms 1 and 2
                        try:
                            # use the generalMapping becuase it contains all the mapping
                            newBond[j] = [a2 for a1, a2 in generalMapping if a1 == bond[j]][0]
                        except IndexError:
                            if pattern['Product'][n_D]['D'] == 'O2x':
                                neighborsM1, _ = FindNeighbors(iStruct['Bonds'], RMList[1])
                                neighborsM2, _ = FindNeighbors(iStruct['Bonds'], RMList[2])
                                connectToM = 0
                                if len(pattern['Reactant']["M"+str(number_R)+'_M1N']) <      \
                                            len(pattern['Product']["M"+str(number_R)+'_M1N']):
                                    # M1 is connected to O2x
                                    if sorted(iStruct['Atoms'][x - 1][1] for x in neighborsM1) == \
                                       sorted(pattern['Reactant']["M"+str(number_R)+'_M1N']):
                                        connectToM = RMList[1]
                                    elif sorted(iStruct['Atoms'][x - 1][1] for x in neighborsM2) == \
                                         sorted(pattern['Reactant']["M"+str(number_R)+'_M1N']):
                                            connectToM = RMList[2]
                                elif len(pattern['Reactant']["M"+str(number_R)+'_M2N']) <    \
                                            len(pattern['Product']["M"+str(number_R)+'_M2N']):
                                    # M2 is connected to O2x
                                    if sorted(iStruct['Atoms'][x - 1][1] for x in neighborsM1) == \
                                       sorted(pattern['Reactant']["M"+str(number_R)+'_M2N']):
                                        connectToM = RMList[1]
                                    elif sorted(iStruct['Atoms'][x - 1][1] for x in neighborsM2) == \
                                         sorted(pattern['Reactant']["M"+str(number_R)+'_M2N']):
                                        connectToM = RMList[2]
                                else:
                                    bRemoveBond = True
                                newBond[j] = connectToM
                            else:
                                bRemoveBond = True
                    if not bRemoveBond: # comment C2
                        tmp_newBond = [newBond[1],newBond[2],newBond[3]]
                        tmp_BondPresent = [[x[1],x[2],x[3]] for x in oStruct['Bonds']]
                        if not tmp_newBond in tmp_BondPresent:
                            oStruct['Bonds'].append(newBond)
                            used_bonds.append(bond)
    
    # Keep track of the charge in the template
    if charge_Template[pattern['KCF']['compound2']]:
        position = []
        for key, value in charge_Template[pattern['KCF']['compound2']].items():
            try:
                if not int(key) in [a1 for a1,a2 in generalMapping]:
                    continue
                position.append([key,[a2 for a1, a2 in generalMapping if a1 == int(key)][0]])
            except:
                position.append([key,[a1 for a1, a2 in pattern['KCF']['Align%d' % number_R] \
                                      if a2 == int(key)][0]])
                
        tmp_charge_Template = deepcopy(charge_Template)
        charge_Template.clear()
        for p in range(len(position)):
            charge_Template[str(position[p][1])] = \
                tmp_charge_Template[pattern['KCF']['compound2']].pop(position[p][0])
    else:
        charge_Template.clear()
        
def XMirrorGroupofPoints(xpos, x_center):
    return [x - 2 * (x - x_center) for x in xpos]

def getThirdVector(O, lenX, A, B, C, D):
    def vecAdd(v1, v2):
        return (v1[0] + v2[0], v1[1] + v2[1])
    def vecSubtract(v1, v2):
        return (v1[0] - v2[0], v1[1] - v2[1])
    def vecScale(v, s):
        return (v[0] * s, v[1] * s)
    def vecRotate(v, theta):
        c, s = cos(theta), sin(theta)
        return (v[0] * c - v[1] * s, v[0] * s + v[1] * c)
    def vecNormalize(v):
        norm = sqrt(v[0] ** 2 + v[1] ** 2)
        if norm == 0:
            return (float('nan'), float('nan'))
        return (v[0] / norm, v[1] / norm)

    if A != (0, 0) and B == (0, 0) and C == (0, 0) and D == (0, 0):
        X = vecAdd(vecRotate(vecSubtract(A, O), 2 * m_pi / 3), O)
    else:
        nA = vecNormalize(vecSubtract(A, O))
        nB = vecNormalize(vecSubtract(B, O))
        nC = vecNormalize(vecSubtract(C, O))
        nD = vecNormalize(vecSubtract(D, O))
        if A == (0, 0):
            nA = (0, 0)
        if B == (0, 0):
            nB = (0, 0)
        if C == (0, 0):
            nC = (0, 0)
        if D == (0, 0):
            nD = (0, 0)
        nX = vecNormalize(vecScale(vecAdd(vecAdd(vecAdd(nA, nB), nC), nD), -1))
        X = vecAdd(O, vecScale(nX, lenX))
    return X

def FindBondBetween2atoms(bond, R, M):
    # finding the index (column1) in bond where elements in col2 and 3 are R and M
    for index, atom1, atom2, order in bond:
        if atom1 == M and atom2 == R:
            return index
    for index, atom1, atom2, order in bond:
        if atom1 == R and atom2 == M:
            return index
    return None

def FindNeighborsP(bonds, atom):
    neighbors = set()
    for index, atom1, atom2, order in bonds:
        if atom1 == atom or atom2 == atom:
            neighbors.add(index)
    return sorted(list(neighbors))

def ConstructMol(inputStructure, charge, charge_Template, substrate_ID, rxnID, comp1, comp2):
    
    def pairwise(iterable):
        a = iter(iterable)
        return zip(a, a)
    
    # constructing a mol file with atoms and bonds information
    mol = ' \n \n \n %2d %2d  0  0  0  0  0  0  0  0999 V2000\n' % \
          (len(inputStructure['Atoms']), len(inputStructure['Bonds']))
    for index, type, element, x, y in inputStructure['Atoms']:
        mol += ('%9s' % ('%.4f' % x)) + ('%10s' % ('%.4f' % y)) + \
               '    0.0000  %s  0 0 0 0 0 0 0 0 0 0 0 0\n' % element
    for index, atom1, atom2, order in inputStructure['Bonds']:
        mol += ' %2d %2d  %d  0     0  0\n' % (atom1, atom2, order)
        
    if charge_Template:
        mol += 'M  CHG  %d' % len(charge_Template)
        for key,value in charge_Template.items():
            mol += '  %2s  %2s' % (key, value)
        mol += '\n'
    
    if not len(charge) == 0:
        if re.findall('M  CHG', mol):
            old = re.findall("M  CHG[\s\S]+",mol)[0].strip()
            old_split = re.split("\s+", old)
            new = "M  CHG  %d" % (int(old_split[2])+ len(charge))
            for x, y in pairwise(old_split[3:]): 
                new += "  %2s  %2s" % (x,y)
                
            for key,value in charge.items():
                new += '  %2s  %2s' % (key, value)
            new += "\n"
            mol.replace(old,new)
        else:
            
            mol += 'M  CHG  %d' % len(charge)
            for key,value in charge.items():
                mol += '  %2s  %2s' % (key, value)
            mol += '\n'
    mol += 'M  END\n'        
    return mol

def ConstructMolFile(inputStructure, charge, charge_Template, filename, product, query_info):    
    # Instead of saving the mol just as text, generate a json with a more complete information within it
    # Store the information as Smiles
    try:
        smiles = Chem.MolToSmiles(Chem.MolFromMolBlock(ConstructMol(inputStructure, 
                                                        charge, charge_Template, query_info['ID'], 
                                                        product['ID'][0], product['KCF']['compound1'], 
                                                        product['KCF']['compound2'])))
        mol = Chem.MolToMolBlock(Chem.MolFromSmiles(smiles))
    except:
        smiles = ""
        mol = ConstructMol(inputStructure, charge, charge_Template, query_info['ID'], product['ID'][0], 
                                   product['KCF']['compound1'], product['KCF']['compound2'])
        
    ec = product['Enzyme']
    rxnID = product['Reaction']
    TemplateSubstrate = product['KCF']['compound1']
    TemplateProduct = product['KCF']['compound2']
    QueryName = query_info['name']
    QuerySmiles = query_info['smiles']
    QueryID = query_info['ID']
    file_to_save = {'GeneratedProduct':[
                        {'smiles':smiles,
                         'mol': mol}
                        ],
                    'TemplateReaction':[
                        {'ec':ec,
                         'ID':rxnID,
                         'Substrate':TemplateSubstrate,
                         'Product':TemplateProduct}
                        ],
                    'QueryInformation': [
                        {'name':QueryName,
                         'ID':QueryID,
                         'smiles':QuerySmiles}
                        ]}
    with open(filename, 'w') as f:
        json.dump(file_to_save, f, indent=2)

def ConstructKcfFile(inputStructure, filename):
    f = open(filename, 'w')
    f.write('ENTRY       C00000                      Compound\n')
    f.write('ATOM        %d\n' % len(inputStructure['Atoms']))
    for index, type, element, x, y in inputStructure['Atoms']:
        f.write('            %-3d %-3s %-2s  %8.4f  %8.4f\n' % (index, type, element, x, y))
    f.write('BOND        %d\n' % len(inputStructure['Bonds']))
    for index, atom1, atom2, order in inputStructure['Bonds']:
        f.write('            %-3d %3d %3d %d\n' % (index, atom1, atom2, order))
    f.write('///\n')
    f.close()

def UpdateNumbering(modifiedInputStructure, charge, charge_Template):
    
    relation = {}
    
    for idx in range(len(modifiedInputStructure['Atoms'])):
        relation[modifiedInputStructure['Atoms'][idx][0]] = idx+1
        
    for atom in modifiedInputStructure['Atoms']:
        atom[0] = relation[atom[0]]
    
    for bond in modifiedInputStructure['Bonds']:
        bond[1] = relation[bond[1]]
        bond[2] = relation[bond[2]]
    
    to_add = []
    to_rem = []
    for each_at in charge.keys():
        to_add.append([relation[int(each_at)],charge[each_at]])
        to_rem.append(each_at)
    
    for r in to_rem:
        charge.pop(r)
    
    for a in to_add:
        charge[a[0]] = a[1]
    
    to_add = []
    to_rem = []
    for each_at in charge_Template.keys():
        to_add.append([relation[int(each_at)],charge_Template[each_at]])
        to_rem.append(each_at)
    
    for r in to_rem:
        charge_Template.pop(r)
    
    for a in to_add:
        charge_Template[a[0]] = a[1]
     
    return modifiedInputStructure, charge, charge_Template

def CheckBondPresence(bondToAdd, oStruct):
    for bond in oStruct['Bonds']:
        if bondToAdd[0] == bond[1] and bondToAdd[1] == bond[2]:
            return True
        elif bondToAdd[1] == bond[1] and bondToAdd[0] == bond[2]:
            return True
    return False
