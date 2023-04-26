import re
from rdkit import Chem

class ProximalException(Exception):
    pass

def Mutate(s):
    if isinstance(s, dict):
        for key, val in s.items():
            s[key] = Mutate(val)
    elif isinstance(s, list):
        for i, val in enumerate(s):
            s[i] = Mutate(val)
    else:
        keggatomTypeOld = ['C2x','C2y','N4x','N2x','N1y']
        keggatomTypeNew = ['C8x','C8y','N1x','N5x','N4y']
        for i in range(len(keggatomTypeOld)):
            if s == keggatomTypeOld[i]:
                return keggatomTypeNew[i]
    return s

def FindNeighbors(bonds, atom):
    # Finds the neighbors of an atom
    neighborsPos = set()
    neighboringAtom = []
    for index, atom1, atom2, order in bonds:
        if atom1 == atom or atom2 == atom:
            neighborsPos.add(index)
            neighboringAtom.append(atom2 if atom1 == atom else atom1)
    return neighboringAtom, sorted(list(neighborsPos))

def modifyAtomList(list_atoms):
    #Since the atoms are extracted through re, the numbers are not int or float
    #so modify them
    result = []
    for x in list_atoms:
        tmp = re.split("\s+", x)
        for y in range(len(tmp)):
            try:
                tmp[y] = int(tmp[y])
            except:
                try:
                    tmp[y] = float(tmp[y])            
                except:
                    continue
        result.append(tmp)
    return result  

def modifyBondList(list_bonds):
    #Since the bonds are extracted through re, the numbers are not int or float
    #so modify them
    result = []
    for x in list_bonds:
        tmp = re.split("\s+", x)
        for y in range(len(tmp)):
            tmp[y] = int(tmp[y])
        result.append(tmp)
    return result

def ExtractCharge(mol_file):
    #Extract charge
    molblock = Chem.MolToMolBlock(mol_file).split("\n")
    
    # Extract from the M  CHG field the charges. the 3rd field represents the number of
    # charged atoms. Then there are paired values (1st the position, 2nd the value of the charge)
    # Then extract all the pairs and store in the charge variable
    def pairwise(iterable):
        a = iter(iterable)
        return zip(a, a)

    charge = {}
    for x in molblock:
        if re.findall('M  CHG', x):
            tmp = re.split("\s+", x)
            for x, y in pairwise(tmp[3:]): 
                charge[x] = y
    return charge

def ExtractPairs(reactions, metabolites):
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
    import pandas as pd
    
    MetabolicSpace = pd.DataFrame(columns = ['id','rxn','Pair','EC'])
    n_rxn = 0
    for idx in reactions.index:
        print("Removing Redundancy among Pairs. Reaction: %d/%d" % (idx, len(reactions)-1))
        substrates = reactions.formula[idx].split(' -> ')[0].strip().split(' + ')
        products = reactions.formula[idx].split(' -> ')[1].strip().split(' + ')
        compoundPair = []
        if len(substrates) > 1 and len(products) > 1:
            # Version to avoid association when all the substrate and products have a low
            # similarity
            mol_sub = {}
            for s in substrates:
                a = Chem.MolFromSmiles(metabolites.loc[metabolites['name'].isin(\
                                                                        [s]),'smiles'].values[0])
                mol_sub[s] = a
            mol_prod = {}
            for p in products:
                a = Chem.MolFromSmiles(metabolites.loc[metabolites['name'].isin(\
                                                                        [p]),'smiles'].values[0])
                mol_prod[p] = a
            
            for s_key,s_value in mol_sub.items():
                f_a = AllChem.GetMorganFingerprint(s_value, 2)
                sim = {}
                for p_key,p_value in mol_prod.items():
                    tmp_fp = AllChem.GetMorganFingerprint(p_value, 2)
                    sim[p_key] = DataStructs.DiceSimilarity(f_a, tmp_fp)
                    if sim[p_key] < 0.5:
                        sim.pop(p_key)
                
                if all([x==0 for x in sim.values()]):
                    continue
                max_sim = max([x for x in sim.values()])
                compoundPair.append([s_key,[x for x,y in sim.items() if y == max_sim][0]])
        else:
            for sub in substrates:
                for prod in products:
                    compoundPair.append([sub,prod])
                    
        for each_p in compoundPair:
            tmp = MetabolicSpace.loc[MetabolicSpace.Pair.isin([each_p]),]
            if len(tmp) > 0:
                MetabolicSpace.loc[tmp.index[0],'rxn'] += [reactions.id[idx]]
                try:
                    tmp_ec = re.findall("[\d-]+\.[\d-]+\.[\d-]+\.[\d-]+",reactions.EC[idx])
                    tmp_ec = MetabolicSpace.loc[tmp.index[0],'EC'] + tmp_ec
                except TypeError:
                    tmp_ec = MetabolicSpace.loc[tmp.index[0],'EC']
                MetabolicSpace.loc[tmp.index[0],'EC'] = list(set(tmp_ec))
            else:
                MetabolicSpace.loc[n_rxn,'id'] = "R%d" % (n_rxn)
                MetabolicSpace.loc[n_rxn,'rxn'] = [reactions.id[idx]]
                MetabolicSpace.loc[n_rxn,'Pair'] = each_p
                try:
                    tmp_ec = re.findall("[\d-]+\.[\d-]+\.[\d-]+\.[\d-]+",reactions.EC[idx])
                except TypeError:
                    tmp_ec = []
                MetabolicSpace.loc[n_rxn,'EC'] = tmp_ec
                n_rxn += 1
    return MetabolicSpace
                