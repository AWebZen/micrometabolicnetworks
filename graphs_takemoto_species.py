#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: adele
"""

from __future__ import unicode_literals

import os
import logging
import cPickle as cpk

import pandas as pd
from bioservices import KEGG
import numpy as np

from utils import *

#Setting logging preferences
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

"""
This code builds graphs for the species from Takemoto et al's 2007 article.
100/113 species were found as of March 2018.
"""



if __name__ == '__main__':
    k = KEGG()
    data = pd.read_excel("12859_2006_1675_MOESM1_ESM.xls") #Takemoto et al suppl file
    species_codes = list(data.iloc[1:114,0]) #113 species
    org_dir = "graph_species/"
    backup_dir = "backup_cpkl/"
    
    if not os.path.exists(org_dir) or not os.path.exists(backup_dir):
        os.makedirs(org_dir)
        os.makedirs(backup_dir)
    
    species_names, dupli, dupli_index = get_species(species_codes) #Species atc not found! (index 96)
    species_names[-1] += " SS9" #Cheat to get right strain for Photobacterium profundum
    optimal_temps = np.array(list(data.iloc[1:114,1]))
    temp_class = np.array(list(data.iloc[1:114,3]))
    
    species_names_all = []
    
    # DOWNLOADING FILES, BUILDING GRAPHS
    for i, code in enumerate(species_codes) :
        print "\n", i, species_codes[i]
        try:
            spec = get_organism_name(code)
            print spec
        except SystemExit:
            species_names_all.append("None")
            logger.error("Organism name not found for code %s" %code)
            continue
            
        try : #Find through assembly ID
            obj = get_fasta(species_codes[i], spec, work_dir=org_dir, alternative_name=species_names[i])
        except SystemExit:
            species_names_all.append("None")
            continue
            
        try:
            obj.directory = org_dir + spec.replace(" ", "_")
            if species_names[i] in dupli:
                obj.directory += "_" + species_codes[i]
            obj.get_reaction_graph(gname="metabolites_reaction_" + species_codes[i] + ".graphml", 
                                   pklname="metabolites_reactions_graph_" + species_codes[i] + ".cpkl",
                                   dir_ec="EC_global_takemoto/")
            obj.build_reaction_graph(filtr=False, 
                                     gname="metabolites_reaction_" + species_codes[i] + ".graphml", 
                                     pklname="metabolites_reactions_graph_" + species_codes[i] + ".cpkl")
            obj.build_substrate_product_graph(gname="metabolites_substrate_product_" + species_codes[i] + ".graphml", 
                                              pklname="metabolites_substrate_product_graph_" + species_codes[i] + ".cpkl")
            obj.build_substrate_product_graph(filtr=False, 
                                              gname="metabolites_substrate_product_" + species_codes[i] + ".graphml", 
                                              pklname="metabolites_substrate_product_graph_" + species_codes[i] + ".cpkl")
            
            #Pathway enzymes graphs
            obj.build_reaction_graph(gname="metabolites_reaction_" + species_codes[i] + "_pathway.graphml", 
                                     pklname="metabolites_reactions_graph_" + species_codes[i] + "_pathway.cpkl", 
                                     pathways=True)
            obj.build_reaction_graph(filtr=False, 
                                     gname="metabolites_reaction_" + species_codes[i] + "_pathway.graphml", 
                                     pklname="metabolites_reactions_graph_" + species_codes[i] + "_pathway.cpkl", 
                                     pathways=True)
            obj.build_substrate_product_graph(gname="metabolites_substrate_product_" + species_codes[i] + "_pathway.graphml", 
                                              pklname="metabolites_substrate_product_graph_" + species_codes[i] + "_pathway.cpkl", 
                                              pathways=True)
            obj.build_substrate_product_graph(filtr=False, 
                                              gname="metabolites_substrate_product_" + species_codes[i] + "_pathway.graphml", 
                                              pklname="metabolites_substrate_product_graph_" + species_codes[i] + "_pathway.cpkl",
                                              pathways=True)
        except (SystemExit, IOError, TypeError): #Problem when constructing
            os.system("rm -f genomes_cdna/" + spec.replace(" ", "_") + "_cds_from_genomic.fna")
            os.system("rm -rf " + org_dir + spec.replace(" ", "_")) #remove directory
            species_names_all.append("None")
            logger.error("Species %s %s will not be handled" %(spec, code))
            continue
            
        species_names_all.append(spec)
        
    cpk.dump(species_names_all, open(backup_dir+"takemoto_species_names.cpkl", "wb"))
    cpk.dump(species_codes, open(backup_dir+"takemoto_species_codes.cpkl", "wb"))
