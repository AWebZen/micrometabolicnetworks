#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: adele
"""
import os
import logging
import cPickle as cpk

import numpy as np
import pandas as pd

from utils import *

"""
Finds new species from Human Pan-Microbiote Communities Database, so mesophiles
from human microbiota (assumed 37ºC), and builds metabolic graphs.

(Adds 76 bacteria as of May 2018 to the ones we already had)
"""

#Human Pan-Microbiote Communities Database
#Roseburia faecis, Bacteroides faecis, Bacteroides finegoldii, Lachnospira pectinoschiza removed for convenience

#Setting logging preferences
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def parse_db_bacteriaclassification(fname):
    """
    Get species names and taxonomy IDs out of a tab separated file, with first
    column "name" and second column "taxonomyid", such as BacteriaClassification.txt
    from Human Pan-Microbiote Communities Database. Handles names with brackets
    (removes brackets from names).
    
    INPUT:
    fname - file name/path
    
    OUTPUT:
    species_names - numpy array of species names, without brackets
    list of taxonomy ids for each species
    """
    data = pd.read_csv(fname, sep = '\t')
    species_names = []
    for name in data["name"]:
        i_close = name.find("]")
        if i_close != -1 :
            if name[0] == "[" :
                new_name = name[1:i_close] + name[i_close+1:]
            else:
                i_open = name.find("[")
                if i_open != -1:
                    new_name = name[:i_open] + name[i_open+1:i_close] + name[i_close+1:]
                else:
                    logger.error("No opening bracket matched to closing bracket!")
                    raise Exception()
            name = new_name
        species_names.append(name)
    return np.array(species_names), list(data['taxonomyid'])


if __name__ == '__main__':
    org_dir = "graph_species/"
    backup_dir = "backup_cpkl/"
    
    if not os.path.exists(org_dir) or not os.path.exists(backup_dir):
        os.makedirs(org_dir)
        os.makedirs(backup_dir)
        
    microbiota_names, microbiota_taxIDs = parse_db_bacteriaclassification('bacteriaClassification_hpmcd_microbiota.txt') #Downloaded from HPMCD data base.
    valid_names = []
    #Downloading species
    for i, name in enumerate(microbiota_names):
        try :
            print "\n%d Species %s" %(i, name)
            code, tax_id = find_get_organism(microbiota_names[i], [microbiota_taxIDs[i]], [], val=False)
            logger.info("code : %s, %s", code, tax_id)
            obj = find_fasta_with_taxID(code, microbiota_names[i], work_dir=org_dir)
#            obj.directory = microbiota_names[i].replace(" ", "_") + "_" + code
            obj.get_reaction_graph(gname = "metabolites_reaction_" + code + ".graphml", pklname = "metabolites_reactions_graph_" + code + ".cpkl", dir_ec="EC_global_takemoto/")
            obj.build_reaction_graph(False, gname = "metabolites_reaction_" + code + ".graphml", pklname = "metabolites_reactions_graph_" + code + ".cpkl")
            obj.build_substrate_product_graph(gname = "metabolites_substrate_product_" + code + ".graphml", pklname = "metabolites_substrate_product_graph_" + code + ".cpkl")
            obj.build_substrate_product_graph(False, gname = "metabolites_substrate_product_" + code + ".graphml", pklname = "metabolites_substrate_product_graph_" + code + ".cpkl")
            #Pathway enzymes graphs
            obj.build_reaction_graph(gname = "metabolites_reaction_" + code + "_pathway.graphml", pklname = "metabolites_reactions_graph_" + code + "_pathway.cpkl", pathways=True)
            obj.build_reaction_graph(False, gname = "metabolites_reaction_" + code + "_pathway.graphml", pklname = "metabolites_reactions_graph_" + code + "_pathway.cpkl", pathways=True)
            obj.build_substrate_product_graph(gname = "metabolites_substrate_product_" + code + "_pathway.graphml", pklname = "metabolites_substrate_product_graph_" + code + "_pathway.cpkl", pathways=True)
            obj.build_substrate_product_graph(False, gname = "metabolites_substrate_product_" + code + "_pathway.graphml", pklname = "metabolites_substrate_product_graph_" + code + "_pathway.cpkl", pathways=True)
            valid_names.append((microbiota_names[i], code))
        except (SystemExit, IOError):
            pass
        
    cpk.dump(valid_names, open(backup_dir + "valid_microbiota.cpkl", "wb"))
