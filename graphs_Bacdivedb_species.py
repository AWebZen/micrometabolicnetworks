#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: adele
"""
import copy
import cPickle as cpk
import os
import logging

import urllib
from bs4 import BeautifulSoup
import numpy as np

from utils import *

"""
Finds new species from BacDive database, especially psychrophiles, thermophiles
and hyperthermophiles. Builds (or re-builds) their graphs.

(Adds 53 species to the 100 from Takemoto et al 2007 article, as of April 2018)
"""

#Setting logging preferences
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_strain_properties(url):
    """
    Gets from an url (of an advanced search in BacDive) the info we want, building
    a series of dictionaries : a global one and one per temperature class, with
    BacDive ID as key and a dictionary as value, with all kinds of attributes as
    keys.
    
    INPUT:
    url - BacDive url (advanced search, as info still appears in page)
    
    OUTPUT:
    db_dict - global dict
    thermo - thermophile dict
    hyper - hyperthermophile dict
    psychro - psychrophile dict
    meso - mesophile dict
    """
    soup = BeautifulSoup(urllib.urlopen(url).read(), "html.parser")
    hits = str(soup.body).split('href="/strain/')[1:] #List of information regarding a species
#    print hits
    if len(hits) == 0:
        logger.error("Uh oh... no species for url!")
        raise SystemExit()
        
    db_dict = dict()
    thermo = {}
    psychro = {}
    hyper = {}
    meso = {}
    for spec in hits :
        hit = BeautifulSoup(spec, "html.parser")
        start = True
        fields = ["strain designation: ", "strain number: ", "synonym: ", 
                  "Species: ", "Temperature: ", "Kind of pH: ", #Field "Species: " became Species name\n(e.g. Escherichia coli)\n:<name>
                  "Temperature range: ", "pH: ", "Kind of pH: ", 
                  "Associated NCBI tax ID: "]
        fields_seen = np.zeros((len(fields)))
        for string in hit.stripped_strings:
            if start:
                start = False
                if len(db_dict.keys()) > 0:
                    if len(db_dict[name]["Temperature range"]) > 0:
                        if db_dict[name]["Temperature range"][0] == "thermophilic":
                            thermo[name] = db_dict[name]
                        elif db_dict[name]["Temperature range"][0] == "hyperthermophilic" :
                            hyper[name] = db_dict[name]
                        elif db_dict[name]["Temperature range"][0] == "psychrophilic" :
                            psychro[name] = db_dict[name]
                        elif db_dict[name]["Temperature range"][0] == "mesophilic" :
                            meso[name] = db_dict[name]
                        else:
                            print "oh no"
                name = string[:string.find('"')]
                if not name.isdigit():
                    logger.warning("Could not find Bacdive ID : %s" %name)
                    continue
                else :
                    db_dict[name] = {'strain designation':[], "synonym":[],
                           "Associated NCBI tax ID":[], 'strain number':[], 
                           'Species':[], "pH":[], "Temperature":[], 
                           "Temperature range":[]}
            else:
                if string.startswith("/* <![CDATA[ */"): #Last species in page
                    break
                if string.startswith("type strain"): #don't care about field
                    continue
                if string.startswith("Kind of pH: "): #don't care about field
                    if fields_seen[4] == 1:
                        fields_seen[7] = 1
                    else:
                        fields_seen[4] = 1
                    continue
#                print string
                for i, field in enumerate(fields):
                    if string.startswith(field):
                        key = field.split(":")[0]
                        if field == "Temperature: " or field == "pH: ":
                            if string == "pH: 4.4.-4.5": #typos in database
                                string = "pH: 4.4-4.5"
                            elif string == "Temperature: 35 .0":
                                string = "Temperature: 35.0"
#                            vals = np.array(string.split(field)[1].split("-")).astype(float)
#                            db_dict[name][key].append(np.mean(vals))  #We take the mean!
                            db_dict[name][key] += string.split(field)[1].split(", ")
                        else:
                            db_dict[name][key] += string.split(field)[1].split(", ")
                        fields_seen[i] = 1
                        found = True
                        break
                    else:
                        found = False
                if not found :
                    try :
                        idx = np.where(fields_seen == 1)[0][-1]
                    except IndexError :
                        logger.warning("Weird line for species %s: %s", spec, string)
                    #Add line that we think has been splitted to last field seen
                    key = fields[idx].split(":")[0]
                    db_dict[name][key][-1] += string
    
    
    if len(db_dict[name]["Temperature range"]) > 0: #handling last element of dictionary
        if db_dict[name]["Temperature range"][0] == "thermophilic":
            thermo[name] = db_dict[name]
        elif db_dict[name]["Temperature range"][0] == "hyperthermophilic" :
            hyper[name] = db_dict[name]
        elif db_dict[name]["Temperature range"][0] == "psychrophilic" :
            psychro[name] = db_dict[name]
        elif db_dict[name]["Temperature range"][0] == "mesophilic" :
            meso[name] = db_dict[name]
        else:
            print "oh no"
    return db_dict, thermo, hyper, psychro, meso      


def correspondance_name_id(d):
    """
    Build a dictionary that makes the correspondance between name of species
    and BacDive id (which is the key to BacDive dictionary)
    """
    new_d = {}
    for key in d.keys():
        new_d[d[key]["Species"][0]] = key
    return new_d


def build_graphs(d, species_temp, validate=False, valid_fname=""):
    # Downloading and building graphs
    valid_list = []
    logger.info("%s ..." %species_temp)
    for name in d.keys():
        try :
            logger.info("Species %s (database ID) %s", name, d[name]["Species"][0])
            code, tax_id = find_get_organism(d[name]["Species"][0], 
                                             d[name]['Associated NCBI tax ID'], 
                                             d[name]['synonym'])
            logger.info("code KEGG, TAXid : %s, %s", code, tax_id)
            obj = find_fasta_with_taxID(code, d[name]["Species"][0])
            obj.directory = work_dir \
                            + d[name]["Species"][0].replace(" ", "_") \
                            + "_" + code
            obj.get_reaction_graph(gname="metabolites_reaction_" + code + ".graphml", 
                                   pklname="metabolites_reactions_graph_" + code + ".cpkl",
                                   dir_ec="EC_global_takemoto/")
            obj.build_reaction_graph(filtr=False, 
                                     gname="metabolites_reaction_" + code + ".graphml", 
                                     pklname="metabolites_reactions_graph_" + code + ".cpkl")
            obj.build_substrate_product_graph(gname="metabolites_substrate_product_" + code + ".graphml", 
                                              pklname="metabolites_substrate_product_graph_" + code + ".cpkl")
            obj.build_substrate_product_graph(filtr=False, 
                                              gname="metabolites_substrate_product_" + code + ".graphml", 
                                              pklname="metabolites_substrate_product_graph_" + code + ".cpkl")
            #Graphs with only enzymes in pathways
            obj.build_reaction_graph(gname="metabolites_reaction_" + code + "_pathway.graphml",
                                     pklname="metabolites_reactions_graph_" + code + "_pathway.cpkl", 
                                     pathways=True)
            obj.build_reaction_graph(filtr=False, 
                                     gname="metabolites_reaction_" + code + "_pathway.graphml", 
                                     pklname="metabolites_reactions_graph_" + code + "_pathway.cpkl", 
                                     pathways=True)
            obj.build_substrate_product_graph(gname="metabolites_substrate_product_" + code + "_pathway.graphml", 
                                              pklname="metabolites_substrate_product_graph_" + code + "_pathway.cpkl", 
                                              pathways=True)
            obj.build_substrate_product_graph(filtr=False, 
                                              gname="metabolites_substrate_product_" + code + "_pathway.graphml", 
                                              pklname="metabolites_substrate_product_graph_" + code + "_pathway.cpkl", 
                                              pathways=True)
            if validate:
                valid_list.append((d[name]["Species"][0], code, name))
        except (SystemExit, IOError):
            pass
    
    if validate:
        cpk.dump(valid_list, open(valid_fname, "wb"))
        return valid_list
        

if __name__ == '__main__':       
    backup_dir = "backup_cpkl/"
    work_dir = "graph_species/" #working directory : directory where we can find our graphs
    
    if not os.path.exists(work_dir) or not os.path.exists(backup_dir):
        os.makedirs(work_dir)
        os.makedirs(backup_dir)
        
    if os.path.exists(backup_dir + "valid_additional_hyper.cpkl") \
        and os.path.exists(backup_dir + "valid_additional_thermo.cpkl") \
        and os.path.exists(backup_dir + "valid_additional_psychro.cpkl"):
        
        #Get names of valid species: lists of tuples (species name, KEGG code, BacDive ID)
        valid_hyper = np.array(cpk.load(open(backup_dir + "valid_additional_hyper.cpkl", "rb"))) 
        valid_thermo = np.array(cpk.load(open(backup_dir + "valid_additional_thermo.cpkl", "rb")))
        valid_psychro = np.array(cpk.load(open(backup_dir + "valid_additional_psychro.cpkl", "rb")))
        additional_species = np.vstack((valid_hyper, valid_thermo, valid_psychro)) #concatenate all species
        
        for i, spec in enumerate(additional_species):
            try :
                org_name, code, name = spec
                code = str(code)
                print "\n%d, Species %s (database ID) %s" %(i, org_name, name)
                if os.path.exists("genomes_cdna2.0/" + org_name.replace(" ", "_") + "_cds_from_genomic.fna"):
                    obj = MetabolicGraph(org_name, "genomes_cdna2.0/" + org_name.replace(" ", "_") + "_cds_from_genomic.fna", code)
                else:
                    logger.warning("Fasta {} not found for species {}! Downloading...".format("genomes_cdna2.0/" + org_name.replace(" ", "_") + "_cds_from_genomic.fna", org_name))
                    obj = find_fasta_with_taxID(code, org_name)
                    logger.info("Fasta found for species %s" %org_name)
                obj.directory = work_dir + org_name.replace(" ", "_") + "_" + code
                obj.get_reaction_graph(gname = "metabolites_reaction_" + code + ".graphml", pklname = "metabolites_reactions_graph_" + code + ".cpkl", dir_ec="EC_global_takemoto/")
                obj.build_reaction_graph(False, gname = "metabolites_reaction_" + code + ".graphml", pklname = "metabolites_reactions_graph_" + code + ".cpkl")
                obj.build_substrate_product_graph(gname = "metabolites_substrate_product_" + code + ".graphml", pklname = "metabolites_substrate_product_graph_" + code + ".cpkl")
                obj.build_substrate_product_graph(False, gname = "metabolites_substrate_product_" + code + ".graphml", pklname = "metabolites_substrate_product_graph_" + code + ".cpkl")
                
                #Pathway enzymes graphs
                obj.build_reaction_graph(gname = "metabolites_reaction_" + code + "_pathway.graphml", pklname = "metabolites_reactions_graph_" + code + "_pathway.cpkl", pathways=True)
                obj.build_reaction_graph(False, gname = "metabolites_reaction_" + code + "_pathway.graphml", pklname = "metabolites_reactions_graph_" + code + "_pathway.cpkl", pathways=True)
                obj.build_substrate_product_graph(gname = "metabolites_substrate_product_" + code + "_pathway.graphml", pklname = "metabolites_substrate_product_graph_" + code + "_pathway.cpkl", pathways=True)
                obj.build_substrate_product_graph(False, gname = "metabolites_substrate_product_" + code + "_pathway.graphml", pklname = "metabolites_substrate_product_graph_" + code + "_pathway.cpkl", pathways=True)
            except (SystemExit, IOError):
                pass
    else:
        #BUILDING GRAPHS
        # Finding valid species and downloading and building their graphs
        
        
        if not os.path.exists(work_dir): #Checking existence
            os.makedirs(work_dir)
        
        global_dict = {} #2060
        thermo_glob = {} #165
        hyper_glob = {} #13
        psychro_glob = {} #106
        meso_glob = {} #1776
        for page in range(21):
            url = "https://bacdive.dsmz.de/advsearch?site=advsearch&searchparams[73][contenttype]=text&searchparams[73][typecontent]=contains&searchparams[73][searchterm]=&searchparams[73][not]=1&searchparams[117][contenttype]=text&searchparams[117][typecontent]=contains&searchparams[117][searchterm]=&searchparams[920][contenttype]=integer&searchparams[920][typecontent]=equal&searchparams[920][searchterm]=&searchparams[920][not]=1&searchparams[919][searchterm]=optimum&searchparams[909][searchterm]=*&searchparams[917][contenttype]=integer&searchparams[917][typecontent]=equal&searchparams[917][searchterm]=&searchparams[917][not]=1&searchparams[916][searchterm]=optimum&searchparams[330][contenttype]=integer&searchparams[330][typecontent]=equal&searchparams[330][searchterm]=&searchparams[330][not]=1&advsearch=search=&pfc="+str(page)
            try:
                db_dict, thermo, hyper, psychro, meso = get_strain_properties(url)
            except SystemExit:
                logger.error("No species for url")
            if page == 12:
                db12 = copy.deepcopy(db_dict)
            global_dict.update(db_dict)
            thermo_glob.update(thermo)
            hyper_glob.update(hyper)
            psychro_glob.update(psychro)
            meso_glob.update(meso)
        
        cpk.dump(global_dict, open(backup_dir + "bacdive_dict.cpkl", "wb"))
        
        
        valid_psychro = build_graphs(psychro_glob, "Psychrophiles", 
                                     validate=True, 
                                     valid_fname=backup_dir + "valid_additional_psychro.cpkl")
                                    #Problem with Granulicella - not right Tax ID
        
        valid_hyper = build_graphs(hyper_glob, "Hyperthermophiles", 
                                 validate=True, 
                                 valid_fname=backup_dir + "valid_additional_hyper.cpkl")


        valid_thermo = build_graphs(thermo_glob, "Thermophiles", 
                                 validate=True, 
                                 valid_fname=backup_dir + "valid_additional_thermo.cpkl")
