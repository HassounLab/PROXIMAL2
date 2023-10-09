# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 12:40:41 2022

@author: fbalzerani
"""
import pickle
from copy import deepcopy

from .Common2 import Mutate, FindNeighbors, ProximalException, modifyAtomList, \
                        modifyBondList, ExtractCharge
                        
from rdkit.Chem import rdFMCS
from rdkit import Chem
try:
    from kcfconvoy import converter
except ImportError:
    from kcfconvoy import KCFvec
    class converter:
        @staticmethod
        def rdkmol_to_kcf(rdkmol, cpd_name='NoName'):
            vec = KCFvec()
            vec.input_rdkmol(rdkmol, cpd_name=cpd_name)
            return vec.kcf
import re

##############################################################################

def GenerateOperators(reactionInfo, operatorsFileName, metabolites):
    
    # getting reaction information
    ecData, charge_Template = ExtractReactionData(reactionInfo, metabolites)
    if len(ecData) == 0:
        raise ProximalException('No transformations found for the substrate-product pair')
    else:
        ecData = ChangeSomeKeggTypeName(ecData)

        # sorts the operators
        operators, operatorsMainAtom, operatorsMainAtomPos = sortRs(ecData)

    f = open(operatorsFileName, 'wb')
    pickle.dump({'operators': operators,
                 'operatorsMainAtom': operatorsMainAtom,
                 'operatorsMainAtomPos': operatorsMainAtomPos,
                 'charge_Template': charge_Template}, 
                 f)
    f.close()

    return operators, operatorsMainAtom, operatorsMainAtomPos, charge_Template

def ExtractReactionData(reactions, metabolites):
    
    #Extract the reaction ID
    rxn = reactions.id
    
    rpairs_mine = []
    charge_Template = {}
        
    id_a = reactions.Pair[0]
    id_b = reactions.Pair[1]
    
    #Extract the Smiles
    smiles_a = metabolites.loc[metabolites['name'].isin([id_a]),'smiles'].values[0]
    smiles_b = metabolites.loc[metabolites['name'].isin([id_b]),'smiles'].values[0]
    
    if id_a == id_b and smiles_a == smiles_b:
        return [],[]
    
    #Generate the mol from the smiles
    a = Chem.MolFromSmiles(smiles_a)
    b = Chem.MolFromSmiles(smiles_b)
    
    #Extract charge information
    charge_a = ExtractCharge(a)
    charge_b = ExtractCharge(b)
    charge_Template[id_a] = charge_a
    charge_Template[id_b] = charge_b
    
    #Calculate Maximum Common Substructure
    match1, target_atm1, match2, target_atm2 = ExtractMCS(a, b)    
    
    if len(match1) == 0 or len(match2) == 0:
        return [],[]
    
    #Obtain the KCF of the molecules
    a_kcf = converter.rdkmol_to_kcf(a)
    b_kcf = converter.rdkmol_to_kcf(b)
    
    #Extract the atoms and bonds information of the 1st molecule
    a_atom = re.findall("ATOM\s+\d+\s+([\S\s]+)BOND",a_kcf)[0].strip().split("\n")
    a_atom = [x.strip() for x in a_atom]
    
    a_molblock = Chem.MolToMolBlock(a)
    tmp_a = re.split("\n", a_molblock)
    min_pos = min([x for x in range(len(tmp_a)) if len(re.split("\s+",tmp_a[x].strip())) == 4])
    max_pos = min([x for x in range(len(tmp_a)) if re.findall("M  ",tmp_a[x])])
    tmp_a = tmp_a[min_pos:max_pos]
    
    #When the enumeration is bigger than 100, the three digits can generate errors.
    #The extraction need to take into account this possibility.
    a_bond = []
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
        a_bond.append(' '.join(tmp_2))
    
    #Extract the atoms and bonds information of the 2nd molecule
    b_atom = re.findall("ATOM\s+\d+\s+([\S\s]+)BOND",b_kcf)[0].strip().split("\n")
    b_atom = [x.strip() for x in b_atom]
    
    b_molblock = Chem.MolToMolBlock(b)
    tmp_b = re.split("\n", b_molblock)
    min_pos = min([x for x in range(len(tmp_b)) if len(re.split("\s+",tmp_b[x].strip())) == 4])
    max_pos = min([x for x in range(len(tmp_b)) if re.findall("M  ",tmp_b[x])])
    tmp_b = tmp_b[min_pos:max_pos]
    
    #Same as before.
    b_bond = []
    count = 0
    for x in tmp_b:
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
        b_bond.append(' '.join(tmp_2))
    
    if len(match1) == len(a_atom) and len(match2) == len(b_atom):
        return [],[]

    #Calculate alignment and reaction centres
    a_atom_proc = [":".join(re.split("\s+",x)[:2]) for x in a_atom]
    b_atom_proc = [":".join(re.split("\s+",x)[:2]) for x in b_atom]
    try:
        a_bond = modifyBondList(a_bond)
        b_bond = modifyBondList(b_bond)
    except ValueError:
        return [],[]
    
    try:
        alignment = RDMatomsByBond(match1,target_atm1,a_bond,
                                            a_atom_proc,match2,target_atm2,b_bond,b_atom_proc)
    except ValueError:
        return [],[]
    
    if len(alignment) == 0:
        return [],[]
    
    New_arrangement = ExtractNewArrangement(match1,a_bond,a_atom_proc,id_a,
                                            match2,b_bond,b_atom_proc,id_b)
    
    
    #Define atoms and bonds of the 2 molecules
    atoms = [modifyAtomList(a_atom), modifyAtomList(b_atom)] 
    bonds = [a_bond, b_bond]
    #Define the ids of the molecules
    compoundIds = [id_a, id_b]
    #Define the ID of pair and reaction
    id_pair = rxn 
    reactionId = reactions.rxn
    ec = reactions.EC
    
    #Define the dictionary related to the rpair
    rpair = {'atoms': atoms,
             'bonds': bonds,
             'compoundIds': compoundIds,
             'id': id_pair,
             'reactionId': reactionId,
             'New_arrangement': New_arrangement,
             'EC':ec 
             }
    
    for idx in range(len(alignment)):
        rpair['alignment%d' % (idx+1)] = alignment[idx+1]
        
    #Add the rpair to the list of transformations
    rpairs_mine.append(rpair)

    rPairKEGG, _ = ReadRDMpattern(rpairs_mine)

    if len(rPairKEGG) == 0:
        return [],[]

    AugmentFields(rPairKEGG)
    ChangeSomeKeggTypeName(rPairKEGG)
    SortMs(rPairKEGG)
    rPairKEGG = RemoveRepeatedRPairs(rPairKEGG)
    return rPairKEGG, charge_Template

def ReadRDMpattern(rpairs):
    # VP: This section implements the file parsing logic in ReadRDMpattern, ReadEntries,
    #     ReadAtom, and ReadBond

    # Read RDM pattern from KEGG

    rPairKEGG = []

    for rpair in rpairs:
        #extract all the alignment present in the reaction (1 or more)
        r_atoms = list(filter(lambda x: "alignment" in x, list(rpair.keys())))
        # VP: This section implements tempEntries construction found in ReadEntries
        tempEntries = {
            'compound1': rpair['compoundIds'][0],
            'compound2': rpair['compoundIds'][1],
            'atom1': rpair['atoms'][0],
            'atom2': rpair['atoms'][1],
            'bond1': rpair['bonds'][0],
            'bond2': rpair['bonds'][1],
            'New_arrangement': rpair['New_arrangement']
        }
        # define RDM for each alignment and add to tempEntries every alignment
        for each_al in r_atoms:
            number = int(re.findall("\d",each_al)[0])
            tempEntries['R%d' % number] =  None
            tempEntries['D%d' % number] = []
            tempEntries['M%d' % number] = []
            atoms = []
            for alignLine in rpair[each_al]:
                atom1Index = int(alignLine[1].split(':')[0]) if alignLine[1] != '*' else -1
                atom2Index = int(alignLine[2].split(':')[0]) if alignLine[2] != '*' else -1
                atom = (atom1Index, atom2Index)
                atoms.append(atom)
    
                if len(alignLine) >= 4:
                    markers = alignLine[3].upper()
                    if ('#R%d' % number) in markers:
                        tempEntries['R%d' % number] = atom
                    elif ('#M%d' % number) in markers:
                        tempEntries['M%d' % number].append(atom)
                    elif ('#D%d' % number) in markers:
                        tempEntries['D%d' % number].append(atom)
    
            tempEntries['Align%d' % number] = list(set(atoms))

        # VP: ReadRDMpattern logic continues below

        rPairKEGG += GenerateRPairKEGG(tempEntries, rpair['id'],
                                       rpair['reactionId'], rpair['EC'])

    return rPairKEGG, None

def GenerateRPairKEGG(tempEntries, rpairId, reactionId, reactionEC):

    rPairKEGG = []
    
    #extract the fields with R atom info and check that there is at least one not None
    tmp_r = list(filter(lambda x: "R" in x, list(tempEntries.keys())))
    
    if len(tempEntries) != 0 and any(tempEntries[x] is not None for x in tmp_r):
        
        #Add the forward reaction
        dataRDM = ReadRDM(tempEntries)
        rPairKEGG.append({
            'ID': [rpairId],
            'Reaction': [reactionId],
            'Enzyme': [reactionEC],
            'Reactant': dataRDM['Reactant'],
            'Product': dataRDM['Product'],
            'KCF': dataRDM['KCF']
        })
        
        #Add the reversible rection
        
        tempEntries2 = SwapReactantProduct(tempEntries)
        dataRDM = ReadRDM(tempEntries2)
        rPairKEGG.append({
            'ID': [rpairId+"_r"],
            'Reaction': [[x+"_r" for x in reactionId]],
            'Enzyme': [reactionEC],
            'Reactant': dataRDM['Reactant'],
            'Product': dataRDM['Product'],
            'KCF': dataRDM['KCF']
        })

    return rPairKEGG

def ReadRDM(entries):
    smallerMolecule=0
    largerMolecule=1
    RPairKegg = {
        'KCF': deepcopy(entries),
        'Reactant':{},
        'Product':{}}
    
    #Calculate the Rs for the molecules
    tmp_r = list(filter(lambda x: "R" in x, list(entries.keys())))
    for each_r in tmp_r:
        RPairKegg['Reactant'][each_r] = entries['atom1'][entries[each_r][smallerMolecule] - 1][1]
        RPairKegg['Product'][each_r] = entries['atom2'][entries[each_r][largerMolecule] - 1][1]
    
    
    #Calculate the Ms for the molecules, coherently to the multi-R
    tmp_m = list(filter(lambda x: "M" in x, list(entries.keys())))
    for each_m in tmp_m:
        number = int(re.findall("\d",each_m)[0])
        for i in range(len(entries[each_m])):
            RPairKegg['Reactant'][each_m+'_M%d' % (i + 1)] = \
                entries['atom1'][entries[each_m][i][smallerMolecule] - 1][1]
            NNr, tempVar = FindNeighbors(entries['bond1'], entries[each_m][i][0])
            neighborsR = [entries['atom1'][n - 1][1] for n in NNr]
    
            neighborsROrdered = list(enumerate(neighborsR))
            neighborsROrdered.sort(key=lambda x: x[1])
            neighborsR = [x[1] for x in neighborsROrdered]
            order = [x[0] for x in neighborsROrdered]
    
            RPairKegg['Reactant'][each_m+'_M%dN' % (i + 1)] = neighborsR
            RPairKegg['Product'][each_m+'_M%d' % (i + 1)] = \
                entries['atom2'][entries[each_m][i][largerMolecule] - 1][1]
    
            # find the corresponding neighbors to the Reactant
            NNr = [NNr[x] for x in order]
            indexP = [None] * len(NNr)
            neighborsP = [None] * len(NNr)
            for z in range(len(NNr)):
                indexR = [x for x in range(len(entries['Align%d' % (number)])) 
                          if entries['Align%d' % (number)][x][0] == NNr[z]]
                if len(indexR) > 1:
                    indexP[z] = None
                    for x in indexR:
                        if entries['Align%d' % (number)][x][1] != -1:
                            assert(indexP[z] is None)
                            indexP[z] = entries['Align%d' % (number)][x][1] - 1
                elif len(indexR) == 1:
                    # note: if atom in alignment is -1, this will be -2, so check for error must be < 0
                    indexP[z] = entries['Align%d' % (number)][indexR[0]][1] - 1
                else:
                    indexP[z] = -1
    
                if indexP[z] < 0:
                    neighborsP[z] = '*'
                else:
                    neighborsP[z] = entries['atom2'][indexP[z]][1]
    
            # find all neighbors for the Product
            allNeighborsP, tempVar = FindNeighbors(entries['bond2'], entries[each_m][i][1])
            addedNeighbor = set(allNeighborsP) - set([x + 1 for x in indexP])
            if len(addedNeighbor) > 0:
                neighborsP += [entries['atom2'][x - 1][1] for x in addedNeighbor]
            RPairKegg['Product'][each_m+'_M%dN' % (i + 1)] = neighborsP
    
    #Calculate the Ds for the molecules, coherently to the multi-R
    tmp_d = list(filter(lambda x: "D" in x, list(entries.keys())))
    for each_d in tmp_d:
        number = int(re.findall("\d",each_d)[0])
        counterAdd = 1
        counterN = len(entries["M%d" % number]) + 1
        for z in range(len(entries[each_d])):
            tempStructure = {}
            atomD = entries[each_d][z][largerMolecule]
    
            if atomD == -1:
                RPairKegg['Reactant']['M'+str(number)+'_M%d' % counterN] = \
                    entries['atom1'][entries[each_d][z][smallerMolecule] - 1][1]
                NN, tempVar = FindNeighbors(entries['bond1'], entries[each_d][z][0])
    
                neighbors = [entries['atom1'][x - 1][1] for x in NN]
                neighbors.sort()
    
                RPairKegg['Reactant']["M"+str(number)+'_M%dN' % counterN] = neighbors
                RPairKegg['Product']["M"+str(number)+'_M%d' % counterN] = '*'
                RPairKegg['Product']["M"+str(number)+'_M%dN' % counterN] = '*'
                RPairKegg['KCF']['M'+str(number)].append(entries[each_d][z])
                RPairKegg['KCF'][each_d] = [x for x in RPairKegg['KCF'][each_d] \
                                            if x != entries[each_d][z]]
                counterN += 1
            else:
                atomR = entries['R%d' % number][largerMolecule]
                bondDetails = entries['bond%d' % (largerMolecule + 1)]
                neighbors, _ = FindNeighbors(bondDetails, atomR)
                if len(entries['M%d' % number]) == 0:
                    otherNeighbors = []
                else:
                    otherNeighbors = [a2 for a1, a2 in entries['M%d' % number]]
    
                RPairKegg['Product'][each_d+'_D%d' % counterAdd] = {
                    'D': entries['atom2'][atomD - 1][1],
                    # 'D': [entries['atom2'][n - 1][1] for n in neighbors],
                    'atomR': atomR,
                    'atomD': atomD,
                    'RDRowinBond': [],
                    'atom': [],
                    'bond': []
                }
                newAtoms = []
                removedAtoms = [atomR] + otherNeighbors
                for i, a in enumerate(entries['atom2']):
                    if a[0] not in removedAtoms:
                        newAtoms.append(deepcopy(a))
                if len(newAtoms) != 0:
                    connectedGraphAtoms, connectedGraphAtomsInfo = \
                        FindConnectedGraphAtoms(atomD, newAtoms, entries['bond2'])
                    if len(connectedGraphAtomsInfo) != 0:
                        connectedGraphBonds, linkRow = \
                            FindConnectedGraphBonds(connectedGraphAtoms, entries['bond2'], 
                                                    atomD, atomR)
                        tempStructure = {
                            'D': entries['atom2'][atomD - 1][1],
                            'atomR': atomR,
                            'atomD': atomD,
                            'RDRowinBond': linkRow,
                            'atom': connectedGraphAtomsInfo,
                            'bond': connectedGraphBonds
                        }
                        RPairKegg['Product'][each_d+'_D%d' % counterAdd] = tempStructure
                        RPairKegg['Reactant'][each_d+'_D%d' % counterAdd] = '*'
                        counterAdd += 1

    return RPairKegg

def SwapReactantProduct(entries):
    swapped = {
        'compound1': deepcopy(entries['compound2']),
        'compound2': deepcopy(entries['compound1']),
        'atom1': deepcopy(entries['atom2']),
        'atom2': deepcopy(entries['atom1']),
        'bond1': deepcopy(entries['bond2']),
        'bond2': deepcopy(entries['bond1']),
        'New_arrangement': deepcopy(entries['New_arrangement'])
        }
    
    tmp_al = list(filter(lambda x: 'Align' in x, list(entries.keys())))
    for each_al in tmp_al:
        number = int(re.findall("\d",each_al)[0])
        swapped['R%d' % number] = (entries['R%d' % number][1], entries['R%d' % number][0])
        swapped['D%d' % number] = [(a2, a1) for a1, a2 in entries['D%d' % number]]
        swapped['M%d' % number] = [(a2, a1) for a1, a2 in entries['M%d' % number]]
        swapped['Align%d' % number] = [(a2, a1) for a1, a2 in entries['Align%d' % number]]
    return swapped

def FindConnectedGraphAtoms(atomD, newAtoms, newBonds):
    # finding the structure connected to atomD, removing other atoms that are
    # not connected to atomD through new atoms bonds

    allAddedAtoms = [a[0] for a in newAtoms]
    connectedGraphAtoms = [atomD]
    j = 0
    connectedGraphAtomsInfo = []
    if atomD in allAddedAtoms:
        for a in newAtoms:
            if a[0] == atomD:
                connectedGraphAtomsInfo.append(a)
                break
        while j < len(connectedGraphAtoms):
            atom = connectedGraphAtoms[j]
            for bond in newBonds:
                if bond[1] == atom:
                    connectedAtom = bond[2]
                elif bond[2] == atom:
                    connectedAtom = bond[1]
                else:
                    continue
                if connectedAtom not in connectedGraphAtoms and connectedAtom in allAddedAtoms:
                    connectedGraphAtoms.append(connectedAtom)
                    for a in newAtoms:
                        if a[0] == connectedAtom:
                            connectedGraphAtomsInfo.append(a)
                            break
            j += 1
    return connectedGraphAtoms, connectedGraphAtomsInfo

def FindConnectedGraphBonds(connectedGraphAtoms, newBonds, atomD, atomR):
    # saves the bonds between the connected structure
    linkRow = None
    connectedGraphBonds = []
    for bond in newBonds:
        atom1 = bond[1]
        atom2 = bond[2]
        if (atom1 == atomD and atom2 == atomR) or (atom2 == atomD and atom1 == atomR):
            linkRow = bond[0]
        if atom1 in connectedGraphAtoms or atom2 in connectedGraphAtoms:
            connectedGraphBonds.append([len(connectedGraphBonds) + 1] + bond[1:4])
    return connectedGraphBonds, linkRow

def AugmentFields(rPairKEGG):
    
    for rpair in rPairKEGG:
        tmp_al = list(filter(lambda x: 'Align' in x, list(rpair['KCF'].keys())))
        for each_al in tmp_al:
            number = int(re.findall("\d",each_al)[0])
            for j in range(1, 5):
                if ("M"+str(number)+'_M%d' % j) not in rpair['Reactant']:
                    rpair['Reactant']["M"+str(number)+'_M%d' % j] = []
                    rpair['Reactant']["M"+str(number)+'_M%dN' % j] = []
                    rpair['Product']["M"+str(number)+'_M%d' % j] = []
                    rpair['Product']["M"+str(number)+'_M%dN' % j] = []
            check_d = list(filter(lambda x: re.findall("D%d" % number + "_D",x), \
                                  list(rpair['Reactant'].keys())))
            if check_d == []:
                rpair['Reactant']["D"+str(number)+'_D1'] = []
                rpair['Product']["D"+str(number)+'_D1'] = {'D': []}
    return rPairKEGG

def ChangeSomeKeggTypeName(rPairKEGG):
    # some kegg atom types have the same discriptions but have different names
    # C8x = C2x
    # C8y = C2y
    # N1x = N4x
    # N5x = N2x
    # N4y = N1y

    for rpair in rPairKEGG:
        rpair['Reactant'] = Mutate(rpair['Reactant'])
        rpair['Product'] = Mutate(rpair['Product'])
        rpair['KCF']['atom1'] = Mutate(rpair['KCF']['atom1'])
        rpair['KCF']['atom2'] = Mutate(rpair['KCF']['atom2'])
    return rPairKEGG

def SortMs(rPairKEGG):
    
    for rpair in rPairKEGG:
        tmp_r = list(filter(lambda x: 'R' in x, rpair['Reactant']))
        for each_r in tmp_r:
            number = int(re.findall("\d",each_r)[0])
            sortKeyPopulated = False
            tempData = {}
            if len(rpair['KCF']['M%d' % number]) > 4:
                tmp_repetition = range(len(rpair['KCF']['M%d' % number]))
            else:
                tmp_repetition = range(4)
            for j in tmp_repetition:
                if len(rpair['Reactant']["M"+str(number)+'_M%d' % (j + 1)]) == 0 or \
                    len(rpair['Reactant']["M"+str(number)+'_M%dN' % (j + 1)]) == 0:
                    break
                sortKeyPopulated = True
                tempData[j] = rpair['Reactant']["M"+str(number)+'_M%d' % (j + 1)] + \
                                ''.join(rpair['Reactant']["M"+str(number)+'_M%dN' % (j + 1)])
            if sortKeyPopulated:
                uniqueData = set(tempData.values())
                if len(uniqueData) != len(tempData):
                    for j in range(len(rpair['KCF']['M%d' % number])):
                        for bond in rpair['KCF']['bond1']:
                            if bond[1] == rpair['KCF']['M%d' % number][j][0] and \
                                bond[2] == rpair['KCF']['R%d' % number][0]:
                                # print(tempData)
                                tempData[j] += str(bond[3])
                        for bond in rpair['KCF']['bond1']:
                            if bond[2] == rpair['KCF']['M%d' % number][j][0] and \
                                bond[1] == rpair['KCF']['R%d' % number][0]:
                                tempData[j] += str(bond[3])
    
                sortedMs = list(tempData.items())
                sortedMs.sort(key=lambda x: x[1])
                order = [x[0] for x in sortedMs]
    
                oldRpair = deepcopy(rpair)
                rpair['KCF']['M%d' % number] = [oldRpair['KCF']['M%d' % number][i] for i in order]
                for i, j in enumerate(order):
                    rpair['Reactant']["M"+str(number)+'_M%d' % (i + 1)] = \
                        oldRpair['Reactant']["M"+str(number)+'_M%d' % (j + 1)]
                    rpair['Reactant']["M"+str(number)+'_M%dN' % (i + 1)] = \
                        oldRpair['Reactant']["M"+str(number)+'_M%dN' % (j + 1)]
                    rpair['Product']["M"+str(number)+'_M%d' % (i + 1)] = \
                        oldRpair['Product']["M"+str(number)+'_M%d' % (j + 1)]
                    rpair['Product']["M"+str(number)+'_M%dN' % (i + 1)] = \
                        oldRpair['Product']["M"+str(number)+'_M%dN' % (j + 1)]

    return rPairKEGG

def RemoveRepeatedRPairs(rPairKEGG):
    # if two RPairs have the same RDM, but different ID, reaction number or enzyme,
    # one is removed and its iformation is added to the other one

    output = [rPairKEGG[0]]
    for i in range(1, len(rPairKEGG)):
        indexSimilarity = IsSimilar(output, rPairKEGG[i])
        if indexSimilarity is not None:
            output[indexSimilarity] = AddRPairsReactionEnzymeInfo(rPairKEGG[i], output[indexSimilarity])
        else:
            output.append(rPairKEGG[i])
    return output

def IsSimilar(input1, input2):
    
    # Finding the similar struct to input2 in input1
    for i, in1 in enumerate(input1):
        flagSimilar = True
        flagSimilar2 = False
        
        
        tmp_r1 = list(filter(lambda x: 'R' in x, list(in1['KCF'].keys())))
        tmp_r2 = list(filter(lambda x: 'R' in x, list(input2['KCF'].keys())))
        
        if len(tmp_r1) != len(tmp_r2):
            flagSimilar = False

        # checking if the reactants are the same
        if flagSimilar:
            for each_r in tmp_r1:
                number = int(re.findall("\d", each_r)[0])
                if in1['Reactant'][each_r] != input2['Reactant'][each_r]:
                    flagSimilar = False
                elif (in1['Reactant']["M"+str(number)+'_M1'] != input2['Reactant']["M"+str(number)+'_M1']) or \
                     (in1['Reactant']["M"+str(number)+'_M2'] != input2['Reactant']["M"+str(number)+'_M2']) or \
                     (in1['Reactant']["M"+str(number)+'_M3'] != input2['Reactant']["M"+str(number)+'_M3']) or \
                     (in1['Reactant']["M"+str(number)+'_M4'] != input2['Reactant']["M"+str(number)+'_M4']):
                    flagSimilar = False
                elif (in1['Reactant']["M"+str(number)+'_M1N'] != input2['Reactant']["M"+str(number)+'_M1N']) or \
                     (in1['Reactant']["M"+str(number)+'_M2N'] != input2['Reactant']["M"+str(number)+'_M2N']) or \
                     (in1['Reactant']["M"+str(number)+'_M3N'] != input2['Reactant']["M"+str(number)+'_M3N']) or \
                     (in1['Reactant']["M"+str(number)+'_M4N'] != input2['Reactant']["M"+str(number)+'_M4N']):
                    flagSimilar = False

        # checking if the products are the same
                if flagSimilar:
                    if in1['Product'][each_r] != input2['Product'][each_r]:
                        flagSimilar = False
                    elif (in1['Product']["M"+str(number)+'_M1'] != input2['Product']["M"+str(number)+'_M1']) or \
                         (in1['Product']["M"+str(number)+'_M2'] != input2['Product']["M"+str(number)+'_M2']) or \
                         (in1['Product']["M"+str(number)+'_M3'] != input2['Product']["M"+str(number)+'_M3']) or \
                         (in1['Product']["M"+str(number)+'_M4'] != input2['Product']["M"+str(number)+'_M4']):
                        flagSimilar = False
                    elif in1['Product']["D"+str(number)+'_D1']['D'] != input2['Product']["D"+str(number)+'_D1']['D']:
                        flagSimilar = False
                    elif len(in1['Product']["D"+str(number)+'_D1']['D']) > 0:
                        if (len(in1['Product']["D"+str(number)+'_D1']['atom']) != \
                            len(input2['Product']["D"+str(number)+'_D1']['atom'])) or\
                           (len(in1['Product']["D"+str(number)+'_D1']['bond']) != \
                            len(input2['Product']["D"+str(number)+'_D1']['bond'])):
                            flagSimilar = False

        # checks for the case of the same neighbors but flipped products
        # R: 'C8x'  =>  R: 'C1y'    and  R: 'C1y'
        # M1: 'C8x'     M1: 'C8x'        M1: 'C1y'
        # M2: 'C8x'     M2: 'C1y'        M2: 'C8x'
                if not flagSimilar:
                    if in1['Reactant'][each_r] == input2['Reactant'][each_r] and \
                        len(input2['Reactant']["M"+str(number)+'_M2']) != 0 and len(in1['Reactant']["M"+str(number)+'_M2']) != 0 and \
                        len(input2['Reactant']["M"+str(number)+'_M3']) == 0 and len(in1['Reactant']["M"+str(number)+'_M3']) == 0 and \
                        len(input2['Reactant']["M"+str(number)+'_M4']) == 0 and len(in1['Reactant']["M"+str(number)+'_M4']) == 0 and \
                        input2['Reactant']["M"+str(number)+'_M1'] == input2['Reactant']["M"+str(number)+'_M2'] and \
                        in1['Reactant']["M"+str(number)+'_M1N'] == input2['Reactant']["M"+str(number)+'_M2N'] and \
                        in1['Reactant']["M"+str(number)+'_M2N'] == input2['Reactant']["M"+str(number)+'_M1N'] and \
                        in1['Product']["M"+str(number)+'_M1'] == input2['Product']["M"+str(number)+'_M2'] and \
                        in1['Product']["M"+str(number)+'_M2'] == input2['Product']["M"+str(number)+'_M1'] and \
                        in1['Product']["D"+str(number)+'_D1']['D'] == input2['Product']["D"+str(number)+'_D1']['D']:
                        flagSimilar2 = True

        if flagSimilar or flagSimilar2:
            return i

    return None

def AddRPairsReactionEnzymeInfo(input1, input2):
    # add RPair, Reaction and Enzyme number of input1 to input2 if not repeated

    output = input2
    # adding RPairs
    for id in input1['ID']:
        if id not in output['ID']:
            output['ID'].append(id)

    # Adding Reaction numbers
    for id in input1['Reaction']:
        if id not in output['Reaction']:
            output['Reaction'].append(id)

    # Adding Enzyme numbers, organism name, organism weight
    for ec in input1['Enzyme']:
        if ec not in output['Enzyme']:
            output['Enzyme'].append(ec)
            if 'Organism' in input1:
                if 'Organism' not in output:
                    output['Organism'] = []
                    output['OrganismWeight'] = []
                output['Organism'].append(input1['Organism'])
                output['OrganismWeight'].append(input1['OrganismWeight'])
                if 'CYPFamilies' in input1:
                    if 'CYPFamilies' not in output:
                        output['CYPFamilies'] = []
                    output['CYPFamilies'].append(input1['CYPFamilies'])
    return output

def sortRs(cypData):
    # This function sorts the sturuct based on phase.Reactant.R
    main = []
    for idx in range(len(cypData)):
        # print(idx)
        # print("###################")
        # print(cypData[idx]['Reactant'])
        # print("###################")
        tmp_r = list(filter(lambda x: 'R' in x, list(cypData[idx]['KCF'].keys())))
        # print([cypData[idx]['Reactant'][each_r] for each_r in tmp_r])
        main.append([cypData[idx]['Reactant'][each_r] for each_r in tmp_r])
        # print("###################")
        
    main = list(enumerate(main))
    main.sort(key=lambda x: x[1])
    sortedR = [x[1] for x in main]
    order = [x[0] for x in main]

    phaseI = [cypData[i] for i in order]
    # phaseImainAtom = sorted(list(set(sortedR)))
    phaseImainAtom = sortedR.copy()
    phaseImainAtomPos = order.copy()

    return phaseI, phaseImainAtom, phaseImainAtomPos

def ExtractMCS(mol1, mol2):
    #Extract the maximum common substructure (MCS)
    
    # 1: CompareOrderExact, in order to avoid misleading between single and double bond (or aromatic)   
    # 2: ringMatchesRingOnly is to assure association between atoms belonging to rings 
    # 3: matchValences is to have differences in a matter of atomic valence
    mcs = rdFMCS.FindMCS([mol1,mol2], bondCompare = rdFMCS.BondCompare.CompareOrderExact,
                          ringMatchesRingOnly = True,
                          matchValences = True)
    
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    #Find the MCS within the molecule 1 and extract the index of the different atoms
    match1 = mol1.GetSubstructMatch(mcs_mol)
    target_atm1 = []
    for atom in mol1.GetAtoms():
        if atom.GetIdx() not in match1:
            target_atm1.append(atom.GetIdx())
    
    #Find the MCS within the molecule 2 and extract the index of the different atoms
    match2 = mol2.GetSubstructMatch(mcs_mol)
    target_atm2 = []
    for atom in mol2.GetAtoms():
        if atom.GetIdx() not in match2:
            target_atm2.append(atom.GetIdx())
    
    return match1, target_atm1, match2, target_atm2

def ExtractRatoms(bond,match):
    #Extract the reaction centres
    r_atom = []
    for each in bond:
        if each[1]-1 in match and each[2]-1 in match:
            continue
        elif not each[1]-1 in match and not each[2]-1 in match:
            continue
        else:
            if each[1]-1 in match:
                r_atom.append(each[1]-1)
            else:
                r_atom.append(each[2]-1)
    return list(set(r_atom))

def ExtractDMatoms(bond,r_atom, match):
    #Extract the adjacent and distant neighbours
    m_atom = []
    d_atom = []
    label_m = []
    label_d = []
    
    for idx in range(len(r_atom)):
        neigh = FindNeighbors(bond,r_atom[idx]+1)[0]
        for n in neigh:
            if n - 1 in match:
                m_atom.append(n)
                label_m.append("#M%d" % (idx+1))
            else:
                d_atom.append(n)
                label_d.append("#D%d" % (idx+1))
    return m_atom, d_atom, label_m, label_d


def RDMatomsByBond(match1,different1,bond1,atom1,match2,different2,bond2,atom2):
    
    """
    Extract the r atoms depending on the modification of bond
    
    From the alignment extract the bonds that contain just one atom present in
    the common substructrue and one that doesn't.
    
    The atom present in the common substructure is the r atom, because that bond
    goes through a modification. (change of FG, or addition of FG)
    
    The losing of FG is taken into account in a second step.
    """
    
    r_atom_1 = ExtractRatoms(bond1,match1)
    r_atom_2 = ExtractRatoms(bond2,match2)
        
    # if there are no r_atom, add the correspondent from the other molecule
    # it can happen when the molecule loses the functional group
    
    for idx in range(len(r_atom_1)):
        if not match2[match1.index(r_atom_1[idx])] in r_atom_2:
            r_atom_2.insert(idx,match2[match1.index(r_atom_1[idx])])
    for idx in range(len(r_atom_2)):
        if not match1[match2.index(r_atom_2[idx])] in r_atom_1:
            r_atom_1.insert(idx,match1[match2.index(r_atom_2[idx])])
    
    # Check the proper corresponding position of the atoms
    sol = [None] * len(r_atom_1)
    for idx in range(len(r_atom_1)):
        try:
            if not r_atom_1[idx] == match1[match2.index(r_atom_2[idx])]:
                sol[idx] = match2[match1.index(r_atom_1[idx])]
            else:
                sol[idx] = r_atom_2[idx]
        except IndexError:
            return []
            
    r_atom_2 = sol
            
    label_r_1 = []
    for idx in range(len(r_atom_1)):
        label_r_1.append("#R%d" % (idx+1))
    
    label_r_2 = []
    for idx in range(len(r_atom_2)):
        label_r_2.append("#R%d" % (idx+1))
    
    
    #Define the m and d atoms based on the neighbours of the r atom
            
    m_atom_1, d_atom_1, label_m_1, label_d_1 = ExtractDMatoms(bond1,r_atom_1,match1)
    m_atom_2, d_atom_2, label_m_2, label_d_2 = ExtractDMatoms(bond2,r_atom_2,match2)
    
    #Define the alignmnet with the label
    # r atoms in position of mol
    # m and d atoms in position of kcf -> +1
    alignment = {}
    for i in range(1,len(label_r_1)+1):
        tmp_al = []
        tmp_m = [m_atom_1[x] for x in range(len(m_atom_1)) 
                     if label_m_1[x] == ("#M%d" % i)]
        tmp_r = [r_atom_1[x] for x in range(len(r_atom_1)) 
                     if label_r_1[x] == ("#R%d" % i)]
        tmp_d_1 = [d_atom_1[x] for x in range(len(d_atom_1)) 
                     if label_d_1[x] == ("#D%d" % i)]
        tmp_d_2 = [d_atom_2[x] for x in range(len(d_atom_2)) 
                     if label_d_2[x] == ("#D%d" % i)]
        for idx in range(len(match1)):
            if match1[idx]+1 in tmp_m:
                tmp_al.append([idx+1,
                                atom1[match1[idx]], 
                                atom2[match2[idx]],
                                "#M%d" % i])
            elif match1[idx] in tmp_r:
                tmp_al.append([idx+1,
                                atom1[match1[idx]], 
                                atom2[match2[idx]],
                                "#R%d" % i])
            else:
                tmp_al.append([idx+1,
                                atom1[match1[idx]], 
                                atom2[match2[idx]]])
        
        for idx in range(len(tmp_d_1)):
            tmp_al.append(["-",
                            atom1[tmp_d_1[idx]-1], 
                            "*",
                            "#D%d" % i])
        
        for idx in range(len(tmp_d_2)):
            tmp_al.append(["-", 
                            "*",
                            atom2[tmp_d_2[idx]-1],
                            "#D%d" % i])
        alignment[i] = tmp_al
                    
    return alignment

def ExtractNewArrangement(match1,bond1,atom1,id_1,match2,bond2,atom2,id_2):
    #Extract the atoms involved in a change of the maximum common substructure (example: from exagon to pentagon)
    bond_present_in_2_not_in_1 = []
    for b2 in bond2:
        if b2[1]-1 in match2 and b2[2]-1 in match2:
            tmp_1 = match1[match2.index(b2[1]-1)]+1
            tmp_2 = match1[match2.index(b2[2]-1)]+1
            flag = True
            for b1 in bond1:
                if tmp_1 == b1[1] and tmp_2 == b1[2]:
                    flag = False
                    break
                if tmp_2 == b1[1] and tmp_1 == b1[2]:
                    flag = False
                    break
                
            if flag:
                bond_present_in_2_not_in_1.append(b2)
    
    bond_present_in_1_not_in_2 = []
    for b1 in bond1:
        if b1[1]-1 in match1 and b1[2]-1 in match1:
            tmp_1 = match2[match1.index(b1[1]-1)]+1
            tmp_2 = match2[match1.index(b1[2]-1)]+1
            flag = True
            for b2 in bond2:
                if tmp_1 == b2[1] and tmp_2 == b2[2]:
                    flag = False
                    break
                if tmp_2 == b2[1] and tmp_1 == b2[2]:
                    flag = False
                    break
                
            if flag:
                bond_present_in_1_not_in_2.append(b1)
                
    bond_to_change_in_1 = []
    for b2 in bond2:
        if b2[1]-1 in match2 and b2[2]-1 in match2:
            tmp_1 = match1[match2.index(b2[1]-1)]+1
            tmp_2 = match1[match2.index(b2[2]-1)]+1
            flag = False
            for b1 in bond1:
                if tmp_1 == b1[1] and tmp_2 == b1[2] and not b1[3] == b2[3]:
                    flag = True
                    break
                if tmp_2 == b1[1] and tmp_1 == b1[2] and not b1[3] == b2[3]:
                    flag = True
                    break
                
            if flag:
                bond_to_change_in_1.append(b2)
    
    bond_to_change_in_2 = []
    for b1 in bond1:
        if b1[1]-1 in match1 and b1[2]-1 in match1:
            tmp_1 = match2[match1.index(b1[1]-1)]+1
            tmp_2 = match2[match1.index(b1[2]-1)]+1
            flag = False
            for b2 in bond2:
                if tmp_1 == b2[1] and tmp_2 == b2[2] and not b1[3] == b2[3]:
                    flag = True
                    break
                if tmp_2 == b2[1] and tmp_1 == b2[2] and not b1[3] == b2[3]:
                    flag = True
                    break
                
            if flag:
                bond_to_change_in_2.append(b1)
                
    New_arrangement = dict({id_1:{},id_2:{}})
    
    # in mol 1 to add
    to_add = []
    for each_new in bond_present_in_2_not_in_1:
        at1 = match1[match2.index(each_new[1]-1)]
        
        at2 = match1[match2.index(each_new[2]-1)]
        
        m_atom_1, _ = FindNeighbors(bond1,at1+1)
        m_atom_2, _ = FindNeighbors(bond1,at2+1)
        
        m_n_1 = []
        for m in m_atom_1:
            m_n_1.append(FindNeighbors(bond1, m)[0])
            
        m_n_2 = []
        for m in m_atom_2:
            m_n_2.append(FindNeighbors(bond1, m)[0])
        
        to_add.append({'R1': at1+1, 'R2': at2+1, 
                       'M1': m_atom_1, 'M1N': m_n_1,
                       'M2': m_atom_2, 'M2N': m_n_2,
                       'order': each_new[3]})
    
    # in mol 1 to remove
    to_remove = []    
    for each_new in bond_present_in_1_not_in_2:
        at1 = each_new[1]-1
        
        at2 = each_new[2]-1
        
        m_atom_1, _ = FindNeighbors(bond1,at1+1)
        m_atom_2, _ = FindNeighbors(bond1,at2+1)
        
        m_n_1 = []
        for m in m_atom_1:
            m_n_1.append(FindNeighbors(bond1, m)[0])
            
        m_n_2 = []
        for m in m_atom_2:
            m_n_2.append(FindNeighbors(bond1, m)[0])
        
        to_remove.append({'R1': at1+1, 'R2': at2+1, 
                          'M1': m_atom_1, 'M1N': m_n_1,
                          'M2': m_atom_2, 'M2N': m_n_2,})
        
    # in mol1 to change
    check = {}
    for each_new in bond_to_change_in_1:
        at1 = match1[match2.index(each_new[1]-1)]
        if at1 in check.keys():
            check[at1] += 1
        else:
            check[at1] = 1
        
        at2 = match1[match2.index(each_new[2]-1)]
        if at2 in check.keys():
            check[at2] += 1
        else:
            check[at2] = 1
    
    to_change = []
    for each_new in bond_to_change_in_1:
        at1 = match1[match2.index(each_new[1]-1)]
        
        at2 = match1[match2.index(each_new[2]-1)]
        
        if check[at1] > 1 and check[at2] > 1:
            continue
        
        m_atom_1, _ = FindNeighbors(bond1,at1+1)
        m_atom_2, _ = FindNeighbors(bond1,at2+1)
        
        m_n_1 = []
        for m in m_atom_1:
            m_n_1.append(FindNeighbors(bond1, m)[0])
            
        m_n_2 = []
        for m in m_atom_2:
            m_n_2.append(FindNeighbors(bond1, m)[0])
        
        to_change.append({'R1': at1+1, 'R2': at2+1, 
                       'M1': m_atom_1, 'M1N': m_n_1,
                       'M2': m_atom_2, 'M2N': m_n_2,
                       'order': each_new[3]})
        
    New_arrangement[id_1]['ToRemove'] = to_remove
    New_arrangement[id_1]['ToAdd'] = to_add
    New_arrangement[id_1]['ToChange'] = to_change
    
    # in mol 2 to add
    to_add = []
    for each_new in bond_present_in_1_not_in_2:
        at1 = match2[match1.index(each_new[1]-1)]
        
        at2 = match2[match1.index(each_new[2]-1)]
        
        m_atom_1, _ = FindNeighbors(bond2,at1+1)
        m_atom_2, _ = FindNeighbors(bond2,at2+1)
        
        m_n_1 = []
        for m in m_atom_1:
            m_n_1.append(FindNeighbors(bond2,m)[0])
            
        m_n_2 = []
        for m in m_atom_2:
            m_n_2.append(FindNeighbors(bond2,m)[0])
        
        to_add.append({'R1': at1+1, 'R2': at2+1, 
                       'M1': m_atom_1, 'M1N': m_n_1,
                       'M2': m_atom_2, 'M2N': m_n_2,
                       'order': each_new[3]})
    
    # in mol 2 to remove
    to_remove = []    
    for each_new in bond_present_in_2_not_in_1:
        at1 = each_new[1]-1
        
        at2 = each_new[2]-1
        
        m_atom_1, _ = FindNeighbors(bond2,at1+1)
        m_atom_2, _ = FindNeighbors(bond2,at2+1)
        
        m_n_1 = []
        for m in m_atom_1:
            m_n_1.append(FindNeighbors(bond2,m)[0])
            
        m_n_2 = []
        for m in m_atom_2:
            m_n_2.append(FindNeighbors(bond2,m)[0])
        
        to_remove.append({'R1': at1+1, 'R2': at2+1, 
                          'M1': m_atom_1, 'M1N': m_n_1,
                          'M2': m_atom_2, 'M2N': m_n_2,})
   
    # in mol2 to change
    check = {}
    for each_new in bond_to_change_in_2:
        at1 = match2[match1.index(each_new[1]-1)]
        if at1 in check.keys():
            check[at1] += 1
        else:
            check[at1] = 1
        
        at2 = match2[match1.index(each_new[2]-1)]
        if at2 in check.keys():
            check[at2] += 1
        else:
            check[at2] = 1
            
    to_change = []
    for each_new in bond_to_change_in_2:
        at1 = match2[match1.index(each_new[1]-1)]
        
        at2 = match2[match1.index(each_new[2]-1)]
        
        if check[at1] > 1 and check[at2] > 1:
            continue
        
        m_atom_1, _ = FindNeighbors(bond2,at1+1)
        m_atom_2, _ = FindNeighbors(bond2,at2+1)
        
        m_n_1 = []
        for m in m_atom_1:
            m_n_1.append(FindNeighbors(bond2,m)[0])
            
        m_n_2 = []
        for m in m_atom_2:
            m_n_2.append(FindNeighbors(bond2,m)[0])
        
        to_change.append({'R1': at1+1, 'R2': at2+1, 
                       'M1': m_atom_1, 'M1N': m_n_1,
                       'M2': m_atom_2, 'M2N': m_n_2,
                       'order': each_new[3]})
        
    New_arrangement[id_2]['ToRemove'] = to_remove
    New_arrangement[id_2]['ToAdd'] = to_add
    New_arrangement[id_2]['ToChange'] = to_change
    
    return New_arrangement
