from proximal_functions.Common2 import ProximalException, ExtractPairs
from proximal_functions.Operators2 import GenerateOperators
from proximal_functions.Products2 import GenerateProducts
from proximal_functions.GenerateMolFiles2 import GenerateMolFiles

from os.path import isfile, isdir
from os import mkdir, makedirs, listdir
import pickle
import pandas as pd
import time
import shutil
from rdkit import Chem

################################# TEST FILES #################################

molecules_of_interest = pd.read_csv("./test/TESTinput/test_molecules_interest.csv")
metabolites = pd.read_csv("./test/TESTinput/test_metabolites.csv")
reaction_list = pd.read_csv("./test/TESTinput/test_reactions.csv")

OP_CACHE_DIRECTORY = 'test/TESToutput/operators'
OUTPUT_DIRECTORY = 'test/TESToutput/products'
path_finalReactions = "./test/TESTinput/"

##############################################################################

########################## APPLICATION FILES #################################
# The molecule of interest must have the structure expressed in smiles.
# It is implemented to accept a csv with tabulator as separator.
# molecules_of_interest = pd.read_csv("./input/XXXX.csv", sep = "\t")
# molecules_of_interest['smiles'] = [Chem.MolToSmiles(Chem.MolFromInchi(x)) for x in molecules_of_interest['Inchi']]
# molecules_of_interest['ID'] = ['Met'+str(x) for x in range(len(molecules_of_interest['name']))]

# The following two are created separately
# metabolites = pd.read_csv("./input/reachableMolecules.csv")
# reaction_list = pd.read_csv("./input/templateReactions.csv")

# OP_CACHE_DIRECTORY = 'output/operators'
# OUTPUT_DIRECTORY = 'output/products'
# path_finalReactions = "./input/"

##############################################################################
#Operator calculation (look-up tables)

if not isdir(OP_CACHE_DIRECTORY):
    makedirs(OP_CACHE_DIRECTORY)
if not isdir(OUTPUT_DIRECTORY):
    makedirs(OUTPUT_DIRECTORY)

# Remove redundancy among pairs
if not isfile(path_finalReactions+"FinalMetabolicSpace.pkl"):
    reactions = ExtractPairs(reaction_list,metabolites)
    with open(path_finalReactions+"FinalMetabolicSpace.pkl", "wb") as f:
        pickle.dump(reactions, f)
else:
    with open(path_finalReactions+"FinalMetabolicSpace.pkl", 'rb') as f: 
        reactions = pickle.load(f)

start_time = time.time()
for i in reactions.index:
    opFilename = OP_CACHE_DIRECTORY + '/' + reactions.id[i] + '.dat'
    print('Building operator: '+str(i)+'/'+str(len(reactions.index)-1))
    if not isfile(opFilename):
        print('Building operator for %s...' % reactions.id[i])
        try:
            operators, operatorsMainAtom, operatorsMainAtomPos, charge_Template = \
                    GenerateOperators(reactions.loc[i,], opFilename, metabolites)
        except ProximalException as e:
            print(str(e))
            continue
print("--- %s seconds = %s hours ---" % ((time.time() - start_time),(time.time() - start_time)/60/60))

##############################################################################
#Product prediction

start_time = time.time()
files = listdir(OP_CACHE_DIRECTORY)
for idx in molecules_of_interest.index:
    targetMoleculeDir = OUTPUT_DIRECTORY + '/' +  \
                                molecules_of_interest['ID'][idx]
    if isdir(targetMoleculeDir):
        shutil.rmtree(targetMoleculeDir)
    mkdir(targetMoleculeDir)
    
    for pos_file in range(len(files)):
        print("Product: " + str(idx) + '/' + str(len(molecules_of_interest.index)-1) + "\n" +
              "Operator: " + str(pos_file) + '/' + str(len(files)-1))
        opFilename = OP_CACHE_DIRECTORY + '/' + files[pos_file]
        f = open(opFilename, 'rb')
        operators = pickle.load(f)['operators']
        f.close()
        f = open(opFilename, 'rb')
        charge_Template = pickle.load(f)['charge_Template']
        f.close()
        
        rxn_id = files[pos_file].split(".")[0]
        try:
            inputList, inputStructure, inputListNumbering, products, charge = \
                    GenerateProducts(molecules_of_interest['name'][idx], 
                                      molecules_of_interest, opFilename, metabolites)
        except ProximalException as e:
            print('%s -> %s (%d ops): %s' % (molecules_of_interest['name'][idx], 
                                              rxn_id, len(operators), str(e)))
            continue
        print('%s -> %s (%d ops) => %d products' %
              (molecules_of_interest['name'][idx], rxn_id, 
                len(operators), len(products)))

        if len(products) >= 1:
            
            targetReactionDir = targetMoleculeDir + '/' + rxn_id

            mkdir(targetReactionDir)
            
            try:
                GenerateMolFiles(targetReactionDir, inputStructure, inputListNumbering, products,
                                  charge, charge_Template, molecules_of_interest.loc[idx,])
                
            except ProximalException as e:
                print(str(e))
                continue
            
            if len(listdir(targetReactionDir)) == 0:
                shutil.rmtree(targetReactionDir)
            
    if len(listdir(targetMoleculeDir)) == 0:
        shutil.rmtree(targetMoleculeDir)
            
print("--- %s seconds = %s hours ---" % ((time.time() - start_time),(time.time() - start_time)/60/60))
