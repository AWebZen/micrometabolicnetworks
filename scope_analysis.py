#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 15:20:08 2019

@author: adele
"""
import logging
import os
import cPickle as cpk
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score
import pandas as pd
import seaborn as sns

from utils import GraphClass
from utils_general import COMPLETE_BACT_MEDIUM, sim_matrix, get_domain,\
    from_y_to_dict, plot_scope_matrix_clusters


#Setting logging preferences
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_graph(name, code, work_dir):
    """
    Load graph for species as GraphClass object.

    Parameters
    ----------
    name : str
        Species name (as used in directory of species)
    code : str
        KEGG organism codes
    work_dir : str
        All metabolic network species directory

    Returns
    -------
    g : GraphClass instance
    """
    g = GraphClass("reaction graph")
                
    if os.path.exists(work_dir + name.replace(" ", "_") + "_" + code): #directory
        directory = work_dir + name.replace(" ", "_") + "_" + code
    else:
        directory = work_dir + name.replace(" ", "_")
        
        
    if os.path.exists(directory+"/metabolites_reaction_" + code + ".graphml"):
        g.load_graph_graphml(directory+"/metabolites_reaction_" + code + ".graphml")
    else:
        logger.error("Graph not found")
        
    if len(g.nodes()) == 0:
        logger.error("Uh oh, empty graph!")
    return g
        

def get_scope_species(g, medium):
    """
    Evaluate scope for each node in graph of a given species and return node names.

    Parameters
    ----------
    g : GraphClass instance of a species

    Returns
    -------
    acc_nodes : list of str
        List of node names in scope for species

    """
    scope_dict = g.scope(inputs=medium)
    acc_true = [n for n, s in scope_dict.items() if s == "Accessible"]
    logger.info("%d/%d nodes in scope" %(len(acc_true), len(scope_dict)))
    acc_nodes = g.get_full_names(acc_true)
    return acc_nodes


def get_scope_all(good_names, good_codes, work_dir, medium):
    """
    Retrieve nodes in scope for each species and build a matrix of presence/absence
    of each node in the scope of each species.

    Parameters
    ----------
    good_names : list/array
        Species names (as used in directories of species)
    good_codes : list/array
        KEGG organism codes
    work_dir : str
        Metabolic network directory

    Returns
    -------
    presence_acc_nodes : numpy array, union of nodes x species
        Matrix of presence (1) or absence (0) of node in scope per species, by
        decreasing order of prevalence in scope for all species
    all_nodes_acc : list of str
        Node names in presence_acc_nodes matrix

    """
    all_scope = []
    
    for i, spec in enumerate(good_names) :        
        print i, spec

        #Get graph     
        g = get_graph(spec, good_codes[i], work_dir)
        #Evaluate scope
        try:
            scope_sp = get_scope_species(g, medium)
        except SystemExit:
            logger.error("Species will have no scope")
            all_scope.append([])
            continue
        
        all_scope.append(scope_sp)
        
    c_nodes = Counter(np.concatenate(all_scope))
    #Orders nodes in decreasing order of number of species with node in scope
    all_nodes_acc = list(np.array(c_nodes.most_common(len(c_nodes)))[:,0])
    
    #Build scope matrix
    presence_acc_nodes = np.zeros((len(all_nodes_acc), len(all_scope)))
    for i, acc in enumerate(all_scope):
        for nod in acc:
            idx = all_nodes_acc.index(nod)
            presence_acc_nodes[idx, i] = 1
    return presence_acc_nodes, all_nodes_acc


def spectral_clustering_scope(simplified_matrix, good_names, good_codes, 
                              optimal_temps, temp_class, domain):
    """
    Do a spectral clustering on scope per species matrix (simplified), for
    2 to 8 clusters. Evaluate silhouette score and save a PDF with cluster
    and species info, a swarmplot of cluster vs optimal growth temperature
    and a scope heatmap with cluster information.

    Parameters
    ----------
    simplified_matrix : array
        Scope nodes x species

    Returns
    -------
    None.

    """
    #Species similarity matrix
    similarity_matrix = sim_matrix(simplified_matrix.T) #Jaccard index

    for clst_no in range(2,8):
        #Spectral clustering
        sc = SpectralClustering(clst_no, affinity='precomputed', n_init=100,
                                assign_labels='discretize')
        sc.fit_predict(similarity_matrix)  
        print list(sc.labels_)
        
        
        d_labels = from_y_to_dict(sc.labels_) #keys=clusters, values=species indexes      
        
        sil = silhouette_score(1 - similarity_matrix, 
                               sc.labels_, metric="precomputed")
        print "k =", clst_no, ":", sil
        
        
        #SAVE INFO AS CSV
        df = pd.DataFrame({"Species_names":good_names,
                           "sp_codes":good_codes, 
                           "opt_temp":optimal_temps, 
                           "temp_class":temp_class, 
                           "clusters": sc.labels_, 
                           "life_domain":domain,
                           })
        df.to_csv("k{}_scope_Takemoto.txt".format(clst_no), index=False)

            
        #PLOTS
        #Cluster x temperature
        plt.figure(figsize=(8, 6))
        col_set = ["orangered", "goldenrod", "lightgreen", "lightblue"]
        sns.swarmplot(x=df.clusters, y=df.opt_temp, hue=df.temp_class, dodge=True, 
                      hue_order=["HT", "T", "M", "P"], size=7, 
                      palette=sns.xkcd_palette(col_set))
        plt.legend(loc="best")
        plt.ylabel("Optimal growth temperature (degrees Celsius)")
        plt.xlabel("Clusters")
        plt.xticks(range(clst_no), [str(cl) for cl in range(clst_no)])
        plt.savefig("temp_vs_clusters_k{}.pdf".format(clst_no))
        

        #Cluster x Scope heatmap
        plot_scope_matrix_clusters(simplified_matrix.T, 
                                   d_labels, 
                                   parameters_txt_xlim=[-30,-35],)
        plt.savefig("scope_heatmap_k{}.pdf".format(clst_no))
        
        

if __name__ == '__main__':
    medium = COMPLETE_BACT_MEDIUM
    
    #LOADING IMPORTANT VARIABLES
    work_dir = "graph_species/" #working directory : directory where we can find our graphs
    backup_dir = "backup_cpkl/"
    if not os.path.exists(work_dir):
        logger.error("Directory with graphs %s could not be found" %work_dir)
        raise SystemExit()
        
    if not os.path.exists(backup_dir):
        logger.error("Backup directory %s could not be found" %backup_dir)
        raise SystemExit()
        
    if not os.path.exists(backup_dir + "/names_228.cpkl"):
        logger.error("Run analysis_graph_properties.py")
        raise SystemExit()
        
    optimal_temps = cpk.load(open(backup_dir + "/optimal_temp_228.cpkl", "rb"))
    good_names = cpk.load(open(backup_dir + "/names_228.cpkl", "rb"))
    good_codes = cpk.load(open(backup_dir + "/codes_228.cpkl", "rb"))
    if not os.path.exists(backup_dir + "/domain_228.cpkl"):
        domain = get_domain(good_codes)
        cpk.dump(np.array(domain), open(backup_dir + "/domain_228.cpkl", "wb"))
    else:
        domain = cpk.load(open(backup_dir + "/domain_228.cpkl", "rb"))
        
    temp_class = np.full((len(optimal_temps), ), "MA") #MA to set size of objects as 2, for HT
    temp_class[temp_class == "MA"] = "M"
    temp_class[optimal_temps <= 25] = "P"
    temp_class[optimal_temps >= 45] = "T"
    temp_class[optimal_temps >= 80] = "HT"
    
    #GET SCOPE FOR ALL SPECIES
    presence_acc_nodes, all_nodes_acc = get_scope_all(good_names, good_codes, 
                                                      work_dir, medium)
    
    
    #REMOVE UNINFORMATIVE NODES (IN/OUT SCOPE FOR ALL SPECIES)
    simplified_matrix = presence_acc_nodes[presence_acc_nodes.sum(axis=1) != presence_acc_nodes.shape[1],:]
    nodes_simplified = np.array(all_nodes_acc)[presence_acc_nodes.sum(axis=1) != presence_acc_nodes.shape[1]]
    
    
    #PLOT OF WHOLE SCOPE MATRIX
    plt.matshow(presence_acc_nodes.T)
    plt.xlabel("Nodes in scope")
    plt.ylabel("Species")
    plt.title("Nodes in scope per species, pathway networks")   

    
    #SPECTRAL CLUSTERING
    spectral_clustering_scope(simplified_matrix, good_names, good_codes, 
                              optimal_temps, temp_class, domain)