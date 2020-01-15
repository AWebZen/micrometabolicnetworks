#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: adele
"""
from __future__ import unicode_literals

import os
import cPickle as cpk
import logging

import numpy as np
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt

from utils import *


 
"""
Laplacian spectrum in decreasing order for our 228 species.
"""

#Setting logging preferences
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

if __name__ == '__main__':
    undirect = False #If we want to work on undirected graphs or not
    pathway = False
    
    #LOADING IMPORTANT VARIABLES
    work_dir = "graph_species/" #working directory : directory where we can find our graphs
    if not os.path.exists(work_dir):
        logger.error("Directory with graphs %s could not be found" %work_dir)
        raise SystemExit()
        
    if not os.path.exists("backup_cpkl"):
        logger.error("Backup directory backup_cpkl/ could not be found")
        raise SystemExit()
        
    if not os.path.exists("backup_cpkl/names_228.cpkl"):
        logger.error("Run analysis_graph_properties.py")
        raise SystemExit()
        
    optimal_temps = cpk.load(open("backup_cpkl/optimal_temp_228.cpkl", "rb"))
    good_names = cpk.load(open("backup_cpkl/names_228.cpkl", "rb"))
    good_codes = cpk.load(open("backup_cpkl/codes_228.cpkl", "rb"))
    temp_class = cpk.load(open("backup_cpkl/temp_class_228.cpkl", "rb"))
    if pathway:
        all_nodes = cpk.load(open("backup_cpkl/all_nodes_pathway.cpkl", "rb")) #union of all nodes (all of our species)
    else:
        all_nodes = cpk.load(open("backup_cpkl/all_nodes.cpkl", "rb")) #union of all nodes (all of our species)
    
    
    #BUILDING SPARSE ADJACENCY MATRICES AND VECTORS
    vectors = [] #list of vectorised nodes * nodes matrices
    class_vectors = [] #list of mean vectorised adjacency matrix for each temperature class
    spectrum = [] #list of sorted Laplacian eigenvalues per species
    class_spectrum = [] #list of mean sorted Laplacian eigenvalues per temperature class
    eigs_nodes = []
    
    for i, spec in enumerate(good_names) : #Our 228 species
        print "\n", i, good_codes[i], spec, temp_class[i]
        try:
            g = GraphClass("substrate-product graph")
            
            if os.path.exists(work_dir + spec.replace(" ", "_") + "_" + good_codes[i]): #directory
                directory = work_dir + spec.replace(" ", "_") + "_" + good_codes[i]
            else:
                directory = work_dir + spec.replace(" ", "_")
            if pathway:
                if os.path.exists(directory+"/metabolites_substrate_product_graph_" + good_codes[i] + "_pathway.cpkl") : #loading graph
                    g.load_graph_cpkl(directory+"/metabolites_substrate_product_graph_" + good_codes[i] + "_pathway.cpkl")
                else:
                    logger.error("Graph not found")
                
            else:
                if os.path.exists(directory+"/metabolites_substrate_product_graph_" + good_codes[i] + ".cpkl") : #loading graph
                    g.load_graph_cpkl(directory+"/metabolites_substrate_product_graph_" + good_codes[i] + ".cpkl")
                else:
                    logger.error("Graph not found")
                
            #BUILDING ADJACENCY MATRICES
            #for sparse matrices (nodes * nodes adjacency matrices)
            rows = []
            cols = []
            data = []
            
            if undirect:
                for edge in g.graph.to_undirected().edges():
                    s,t = [all_nodes.index(g.node_name_equivalence[x]) for x in edge]
                    rows.append(s)
                    cols.append(t)
                    data.append(1)
                    if s != t:
                        rows.append(t)
                        cols.append(s)
                        data.append(1)
            else:
                for edge in g.graph.edges():
                    s,t = [all_nodes.index(g.node_name_equivalence[x]) for x in edge]
                    rows.append(s)
                    cols.append(t)
                    data.append(1)
            
            #adjacency matrix
            mat = coo_matrix((data, (rows, cols)), shape=(len(all_nodes), len(all_nodes)))
         
            #LAPLACIAN MATRIX AND EIGENVALUES
            norm = True #Normalised laplacian or not
            lapl, spect = laplacian_spectrum(mat, normalised=norm, use_out_degree=False)
            #Uncomment to study eigenvectors
#            lapl, spect, vects = laplacian_spectrum_vector(mat, normalised=norm, use_out_degree=False)
#            eigs_nodes.append(get_nodes_eigenvalue(spect, vects, all_nodes))
            
            spectrum.append(np.sort(spect)[::-1]) #sort in decreasing order.
            

            
            #PER TEMPERATURE CLASS
            if i == len(temp_class)-1 or temp_class[i] != temp_class[i+1] :
                index_class = np.where(temp_class == temp_class[i])[0] #Index of species in class
                
                #CLASS SPECTRUM
                spect_class = np.mean(np.array(spectrum)[index_class], axis=0)
                class_spectrum.append(spect_class)
                
        except SystemExit :
            logger.error("Uh oh, something went wrong with this species...")
            pass
        
    eigen = {tmp_cl:class_spectrum[j] for j, tmp_cl in enumerate(np.unique(temp_class))}
    
    if pathway:
        if undirect:
            if norm:
                cpk.dump(np.array(spectrum), open("backup_cpkl/spectrum_Norm_undir_pathway.cpkl", "wb"))
            else:
                cpk.dump(np.array(spectrum), open("backup_cpkl/spectrum_NotNorm_undir_pathway.cpkl", "wb"))
        else:
            if norm:
                cpk.dump(np.array(spectrum), open("backup_cpkl/spectrum_Norm_dir_pathway.cpkl", "wb"))
            else:
                cpk.dump(np.array(spectrum), open("backup_cpkl/spectrum_NotNorm_dir_pathway.cpkl", "wb"))
    else:
        if undirect:
            if norm:
                cpk.dump(np.array(spectrum), open("backup_cpkl/spectrum_Norm_undir.cpkl", "wb"))
            else:
                cpk.dump(np.array(spectrum), open("backup_cpkl/spectrum_NotNorm_undir.cpkl", "wb"))
        else:
            if norm:
                cpk.dump(np.array(spectrum), open("backup_cpkl/spectrum_Norm_dir.cpkl", "wb"))
            else:
                cpk.dump(np.array(spectrum), open("backup_cpkl/spectrum_NotNorm_dir.cpkl", "wb"))



    #Plot sorted eigenvalues per class
    plt.close()
    plt.figure()
    colors = ["b", "g", "r", "y"]
    colors2 = ["steelblue", "lime", "salmon", "orange"]
    t_class = ['P', 'M', 'T', 'HT']
    for j, tmp_cl in enumerate(t_class):
        plt.scatter(range(len(spect_class)), eigen[tmp_cl], c=colors[j], lw=0, label = tmp_cl, s=20)
        plt.legend()
        plt.ylabel("Sorted eigenvalues")
        #    plt.show()
        plt.xlabel('Rank')
        plt.savefig("eigen.pdf")    
        plt.close()
        for j, tmp_cl in enumerate(t_class):
                plt.scatter(range(80),eigen[tmp_cl][:80], c=colors[j], lw=0, label = tmp_cl, s=20)
        #plt.legend()
        plt.ylabel("Sorted eigenvalues")
        plt.xlabel('Rank')
        plt.savefig('eigen_close_up.pdf')
        #    plt.show()
    