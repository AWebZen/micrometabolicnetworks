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
import numpy as np
import networkx as nx

from utils import get_temperature_pH_bacdive, GraphClass, triadic_census, article_density, \
    degree_exponent_ML
from utils_general import is_enz, plotting_pathway_non_pathway_graph, plotting_graph
    

#Setting logging preferences
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

"""
Main analysis code, generating lots of variables and back_ups needed for other
codes.

Analyses graph properties for prokaryotic species, plots against optimal growth
temperature, etc.
"""

    

if __name__ == '__main__':
    work_dir = "graph_species/" #working directory : directory where we can find our graphs
    if not os.path.exists(work_dir):
        logger.error("Working directory %s (where graphs are located) does not exist. Modify." %work_dir)
        raise SystemExit()
        
        
    #SPECIES FROM TAKEMOTO et al. 2007 ARTICLE
    data = pd.read_excel("12859_2006_1675_MOESM1_ESM.xls")
    species_codes = list(data.iloc[1:114,0]) #113 species
    optimal_temps = np.array(list(data.iloc[1:114,1]))
    temp_class = np.array(list(data.iloc[1:114,3]))
    species_names = np.array(cpk.load(open("backup_cpkl/takemoto_species_names.cpkl", "rb")))



    # FETCHING DATA FOR SPECIES FROM BACDIVE AND HPMCD DATABASES
    if not os.path.exists("backup_cpkl/valid_additional_hyper.cpkl")\
    or not os.path.exists("backup_cpkl/valid_additional_hyper.cpkl")\
    or not os.path.exists("backup_cpkl/valid_additional_thermo.cpkl")\
    or not os.path.exists("backup_cpkl/valid_additional_psychro.cpkl")\
    or not os.path.exists("backup_cpkl/valid_microbiota.cpkl"):
        logger.error("Has data from BacDive and HPMCD databases been generated? Backups not found.")
        raise SystemExit()
    #BacDive
    valid_hyper = np.array(cpk.load(open("backup_cpkl/valid_additional_hyper.cpkl", "rb"))) #List of tuples (species name, KEGG code, BacDive ID)
    valid_thermo = np.array(cpk.load(open("backup_cpkl/valid_additional_thermo.cpkl", "rb")))
    valid_psychro = np.array(cpk.load(open("backup_cpkl/valid_additional_psychro.cpkl", "rb")))
    #HPMCD 
    valid_meso = np.array(cpk.load(open("backup_cpkl/valid_microbiota.cpkl", "rb")))
    dupes = np.where(np.in1d(valid_meso[:,1], species_codes))[0] #Remove species common with article
    valid_meso = np.delete(valid_meso, dupes, 0)
    valid_meso = np.hstack((valid_meso, np.ones((len(valid_meso),1)))) #To have same format as BacDive species

    additional_species = np.vstack((valid_hyper, valid_thermo, valid_meso, valid_psychro))
    
    #Temperature for additional species
    temp_class_additions = np.concatenate((np.full((valid_hyper.shape[0],),"HT"), 
                                           np.full((valid_thermo.shape[0],),"T"), 
                                           np.full((valid_meso.shape[0],),"M"), 
                                           np.full((valid_psychro.shape[0],),"P"),
                                           ))
    bacdive_dict = cpk.load(open("backup_cpkl/bacdive_dict.cpkl", "rb")) #Database information : {Bacdive ID : {Species:<>, Temperature:<>, pH:<>, etc}}
    temp_additions = []
    for valid_spec in [valid_hyper, valid_thermo, valid_meso, valid_psychro] :
        if valid_spec is not valid_meso :
            temp, _ = get_temperature_pH_bacdive(bacdive_dict, valid_spec)
        else:
            temp = [37. for _ in range(len(valid_meso))]
        temp_additions += temp
        
    
    
    # MERGING OLD DATA AND NEW DATA
    tmp = np.concatenate((temp_class, temp_class_additions))
    sort_old_new = np.argsort(tmp) #sort per temperature class
    temp_class_global = tmp[sort_old_new]
    species_codes_global = np.concatenate((species_codes, additional_species[:,1]))[sort_old_new]
    species_names_global = np.concatenate((species_names, additional_species[:,0]))[sort_old_new]
    optimal_temps_global = np.concatenate((optimal_temps, np.array(temp_additions)))[sort_old_new]
    
    
    #==============================================================================
    # GRAPH ANALYSIS
    #==============================================================================
    density_undir = []
    strong_comp_g1 = []
    strong_comp_g2 = []
    strong_comp_g4 = []
    strong_comp_g5 = []
    weak_g2 = []
    length_comp = []
    assortativity = []
    transitivity = []
    flow_hier = []
    triads = []
    nodes = []
    exponent = []
    clustering_coeff = []
    closeness = []
    betweenness = []
    in_degr = []
    out_degr = []
    edges = []
    enzs = []
    
    valid_species = [] #valid species indices (species with existing graphs)
    
    for i, spec in enumerate(species_names_global) :
        print "\n", i, species_codes_global[i], spec
        try:
            g1 = GraphClass("reaction graph")
            g2 = GraphClass("substrate-product graph")
            g4 = GraphClass("pathway reactions graph")
            g5 = GraphClass("pathway substrate-product graph")
            
            if os.path.exists(work_dir + spec.replace(" ", "_") + "_" + species_codes_global[i]): 
                directory = work_dir + spec.replace(" ", "_") + "_" + species_codes_global[i] #BacDive species have KEGG code in the directory name
            elif os.path.exists(work_dir + spec.replace(" ", "_")):
                directory = work_dir + spec.replace(" ", "_")
            
            if os.path.exists(directory+"/metabolites_reaction_" + species_codes_global[i] + ".graphml") \
            and os.path.exists(directory+"/metabolites_reaction_" + species_codes_global[i] + "_pathway.graphml") \
            and os.path.exists(directory+"/metabolites_substrate_product_" + species_codes_global[i] + ".graphml") \
            and os.path.exists(directory+"/metabolites_substrate_product_" + species_codes_global[i] + "_pathway.graphml") :
                g1.load_graph_graphml(directory+"/metabolites_reaction_" + species_codes_global[i] + ".graphml")
                g4.load_graph_graphml(directory+"/metabolites_reaction_" + species_codes_global[i] + "_pathway.graphml")
                g2.load_graph_graphml(directory+"/metabolites_substrate_product_" + species_codes_global[i] + ".graphml")
                g5.load_graph_graphml(directory+"/metabolites_substrate_product_" + species_codes_global[i] + "_pathway.graphml")
                
                #Robustness testing
                perc = 0
                g1.remove_nodes(perc)
                g2.remove_nodes(perc)
                g4.remove_nodes(perc)
                g5.remove_nodes(perc)
                
                #Degree - edge density
                density_undir.append([article_density(g1.graph.to_undirected()), 
                                      article_density(g2.graph.to_undirected()), 
                                      article_density(g4.graph.to_undirected()), 
                                      article_density(g5.graph.to_undirected()),
                                      ])
                                      
                # 0 in in degrees or out degrees
                in_degr.append([np.sum(np.array(dict(g1.graph.in_degree).values()) == 0), 
                                np.sum(np.array(dict(g2.graph.in_degree).values()) == 0),
                                np.sum(np.array(dict(g4.graph.in_degree).values()) == 0),
                                np.sum(np.array(dict(g5.graph.in_degree).values()) == 0),
                                ])
                out_degr.append([np.sum(np.array(dict(g1.graph.out_degree).values()) == 0), 
                                 np.sum(np.array(dict(g2.graph.out_degree).values()) == 0),
                                 np.sum(np.array(dict(g4.graph.out_degree).values()) == 0),
                                 np.sum(np.array(dict(g5.graph.out_degree).values()) == 0),
                                 ])
    
                edges.append([g1.graph.number_of_edges(), 
                              g2.graph.number_of_edges(),
                              g4.graph.number_of_edges(),
                              g5.graph.number_of_edges(),
                              ])
    
                enzs.append([len([n for n in g1.nodes() if is_enz(n)]),
                             len([n for n in g4.nodes() if is_enz(n)]),
                             ])
                

                # Strongly connected components SCC
                comps_1 = g1.strongly_connected_comps()
                comps_2 = g2.strongly_connected_comps()
                comps_4 = g4.strongly_connected_comps()
                comps_5 = g5.strongly_connected_comps()
                weak_2 = g2.weakly_connected_comps()
                strong_comp_g1.append(comps_1[1])
                strong_comp_g2.append(comps_2[1])
                strong_comp_g4.append(comps_4[1])
                strong_comp_g5.append(comps_5[1])
                weak_g2.append(weak_2[1])
                length_comp.append([comps_1[0], 
                                    comps_2[0], 
                                    comps_4[0], 
                                    comps_5[0]
                                    ])
    
                #build/save strongly and weakly connected component subgraphs
#                g2.scc_wcc_subgraph(species_codes_global[i], directory)
                
                #Degree exponent
                exponent.append([degree_exponent_ML(g1.graph), 
                                 degree_exponent_ML(g1.graph.to_undirected()), 
                                 degree_exponent_ML(g2.graph), 
                                 degree_exponent_ML(g2.graph.to_undirected()), 
                                 degree_exponent_ML(g4.graph.to_undirected()), 
                                 degree_exponent_ML(g5.graph.to_undirected())
                                 ])
                
                #Average node clustering coefficient
                clustering_coeff.append([nx.average_clustering(g1.graph.to_undirected()), 
                                         nx.average_clustering(g2.graph.to_undirected()), 
                                         nx.average_clustering(g4.graph.to_undirected()), 
                                         nx.average_clustering(g5.graph.to_undirected())
                                         ])

                
                #Many NetworkX measures
                assortativity.append([nx.degree_assortativity_coefficient(g1.graph), 
                                      nx.degree_assortativity_coefficient(g2.graph), 
                                      ])
                transitivity.append(nx.transitivity(g2.graph)) #0 for the 2 other graphs
                triads.append([triadic_census(g1.graph), 
                               triadic_census(g2.graph)
                               ]) #32 plots!

                #Flow hierarchy                               
                flow_hier.append([g2.hierarchy_flow(),
                                 g5.hierarchy_flow(),
                                 ]) #our function
                
                #Number of nodes
                nodes.append([len(g1.graph.nodes), 
                              len(g2.graph.nodes),
                              len(g4.graph.nodes), 
                              len(g5.graph.nodes),
                              ])
    

                #Centralities - closeness and betweenness centralities
                closeness.append(g2.centrality(20, centrality = nx.closeness_centrality))
                betweenness.append(g2.centrality(20, centrality = nx.betweenness_centrality))
    
    
                if i == 0:
                    all_nodes_set = set(g2.node_name_equivalence.values()) #Union of all nodes in all graphs g2 (substrate-product filtered all enz)
                    all_nodes_set_path = set(g5.node_name_equivalence.values()) #Union of all nodes in all graphs (substrate-product filtered pathway enz)
                    
                else:
                    all_nodes_set = all_nodes_set.union(set(g2.node_name_equivalence.values()))
                    all_nodes_set_path = all_nodes_set_path.union(set(g5.node_name_equivalence.values()))
                
                valid_species.append(i)
                
        except SystemExit :
            pass
            
    #Filtering prokaryotes that did not work when creating graphs
    optimal_temps_ok = optimal_temps_global[valid_species] #Values for prokaryotes for whom we have graphs
    temp_class_ok = temp_class_global[valid_species] #Values for prokaryotes for whom we have graphs
    good_names = np.array(species_names_global)[valid_species] #updating lists with valid species
    good_codes = np.array(species_codes_global)[valid_species]
    
    #SAVING IMPORTANT VARIABLES
    save = True
    if save :
        #We will sort per temperature the species for further studies
        index_temp = np.argsort(optimal_temps_ok)
        good_names_sort = good_names[index_temp]
        optimal_temps_sort = optimal_temps_ok[index_temp]
        good_codes_sort = good_codes[index_temp]
        temp_class_sort = temp_class_ok[index_temp]
        
        number = str(len(valid_species))
        cpk.dump(optimal_temps_sort, open("backup_cpkl/optimal_temp_"+number+".cpkl", "wb"))
        cpk.dump(good_names_sort, open("backup_cpkl/names_"+number+".cpkl", "wb"))
        cpk.dump(good_codes_sort, open("backup_cpkl/codes_"+number+".cpkl", "wb"))
        cpk.dump(temp_class_sort, open("backup_cpkl/temp_class_"+number+".cpkl", "wb"))
        cpk.dump(list(all_nodes_set), open("backup_cpkl/all_nodes_"+number+".cpkl", "wb"))
        cpk.dump(list(all_nodes_set_path), open("backup_cpkl/all_nodes_pathway_"+number+".cpkl", "wb"))
        
    #Arraying to make our life easier
    density_undir = np.array(density_undir)
    length_comp = np.array(length_comp)
    assortativity = np.array(assortativity)
    transitivity = np.array(transitivity)
    flow_hier = np.array(flow_hier)
    triads = np.array(triads)
    nodes = np.array(nodes)
    exponent = np.array(exponent)
    in_degr = np.array(in_degr)
    out_degr = np.array(out_degr)
    clustering_coeff = np.array(clustering_coeff)
    edges = np.array(edges)
    enzs = np.array(enzs)
    
    #Largest strongly connected component
    max_1 = []
    max_2 = []
    max_4 = []
    max_5 = []
    for i in range(len(strong_comp_g1)) :
        max_1.append(np.max(strong_comp_g1[i]))
        max_2.append(np.max(strong_comp_g2[i]))
        max_4.append(np.max(strong_comp_g4[i]))
        max_5.append(np.max(strong_comp_g5[i]))

    #Number of enzymes per species
    nb_enzs_all = nodes[:,0]-nodes[:,1]
    nb_enzs_path = nodes[:,2]-nodes[:,3]
    

    #==========================================================================
    #     SOME PLOTS
    #==========================================================================
    
    #NODES
    plotting_pathway_non_pathway_graph(optimal_temps_ok, 
                                       nodes[:,1], 
                                       nodes[:,3], 
                                       "Temperature (°C)", "Number of nodes")
    
    #EDGE DENSITY   
    plotting_pathway_non_pathway_graph(optimal_temps_ok, 
                                       density_undir[:,1], 
                                       density_undir[:,3], 
                                       "Temperature (°C)", "Edge density")

    #CLUSTERING COEFFICIENT
    plotting_graph(optimal_temps_ok, 
                   clustering_coeff[:,1], 
                   "Temperature (°C)", 
                   "Average clustering coefficient\nfor undirected substrate-product graph", .005)
    plotting_pathway_non_pathway_graph(optimal_temps_ok, 
                                       clustering_coeff[:,1], 
                                       clustering_coeff[:,3], 
                                       "Temperature (°C)", "Average clustering coefficient")
    
    #DEGREE EXPONENT
    plotting_graph(optimal_temps_ok, 
                   exponent[:,0], 
                   "Temperature (°C)", 
                   "Degree exponent for directed reaction graph", .015)
    plotting_graph(optimal_temps_ok, 
                   exponent[:,1], 
                   "Temperature (°C)", 
                   "Degree exponent for undirected reaction graph", .015)
    plotting_graph(optimal_temps_ok, 
                   exponent[:,2], 
                   "Temperature (°C)", 
                   "Degree exponent for directed substrate-product graph", .02)
    plotting_graph(optimal_temps_ok, 
                   exponent[:,3], 
                   "Temperature (°C)", 
                   "Degree exponent for undirected substrate-product graph", .02)
    # for filtered undirected substrate-product graph
    plotting_pathway_non_pathway_graph(optimal_temps_ok, 
                                       exponent[:,3], 
                                       exponent[:,5], 
                                       "Temperature (°C)", "Degree exponent")

    #STRONGLY CONNECTED COMPONENTS
    plotting_graph(optimal_temps_ok, 
                   length_comp[:,0], 
                   "Temperature (°C)", 
                   "Number of strongly connected comps for reaction graph", 1)
    plotting_graph(optimal_temps_ok, 
                   length_comp[:,1], 
                   "Temperature (°C)", 
                   "Number of strongly connected comps\nfor substrate-product graph", 0)
    plotting_pathway_non_pathway_graph(optimal_temps_ok, 
                                       length_comp[:,1], 
                                       length_comp[:,3], 
                                       "Temperature", "Number of strongly connected comps\nfor substrate-product graph")
    #Largest SCC    
    plotting_graph(optimal_temps_ok, 
                   max_1/nodes[:,0].astype(float), 
                   "Temperature (°C)", 
                   "Normalised maximal length of strongly connected comps\nfor reaction graph", .001)
    plotting_graph(optimal_temps_ok, 
                   max_2/nodes[:,1].astype(float), 
                   "Temperature (°C)", 
                   "Normalised maximal length of strongly connected comps\nfor substrate-product graph", .001)
    #filtered directed substrate-product graph
    plotting_pathway_non_pathway_graph(optimal_temps_ok, 
                                       max_2/nodes[:,1].astype(float), 
                                       max_5/nodes[:,3].astype(float), 
                                       "Temperature (°C)", 
                                       "Node normalised size of\nlargest strongly connected component")
    
    #FLOW HIERARCHY
    plotting_graph(optimal_temps_ok, 
                   flow_hier[:,0], 
                   "Temperature (°C)", 
                   "Node flow hierarchy for substrate-product graph", .01)
    plotting_graph(optimal_temps_ok, 
                   flow_hier[:,0]/nodes[:,1], 
                   "Temperature (°C)", 
                   "Node normalised flow hierarchy\nfor substrate-product graph", .0005, ylim=[0, 0.0035])
    #filtered directed substrate-product graph
    plotting_pathway_non_pathway_graph(optimal_temps_ok, 
                                       flow_hier[:,0]/nodes[:,1].astype(float), 
                                       flow_hier[:,1]/nodes[:,3].astype(float), 
                                       "Temperature (°C)", "Node normalised flow hierarchy", ylim = [0, 0.004])
    
    #ZERO IN/OUT DEGREE
    plotting_pathway_non_pathway_graph(optimal_temps_ok, 
                                       in_degr[:,1]/(nodes[:,1]).astype(float), 
                                       in_degr[:,3]/(nodes[:,3]).astype(float), 
                                       "Temperature (°C)", "Fraction of 0 in-degrees")
    plotting_pathway_non_pathway_graph(optimal_temps_ok, 
                                       out_degr[:,1]/(nodes[:,1]).astype(float), 
                                       out_degr[:,3]/(nodes[:,3]).astype(float), 
                                       "Temperature (°C)", "Fraction of 0 out-degrees")
    
    #NETWORKX MEASURES
    plotting_graph(optimal_temps_ok, 
                   assortativity[:,0], 
                   "Temperature (°C)", 
                   "Degree assortativity for reaction graph", .01)
    plotting_graph(optimal_temps_ok, 
                   assortativity[:,1], 
                   "Temperature (°C)", 
                   "Degree assortativity for substrate-product graph", .01)
    
    plotting_graph(optimal_temps_ok, 
                   transitivity, 
                   "Temperature (°C)", 
                   "Transitivity for substrate-product graph", .01) 
                    #fraction of all possible triangles present in G.
    
    """
    #TRIAD FRACTIONS (like subgraphs in Takemoto et al's 2007 article), UNCOMMENT WHEN WANTING TO SEE THE FOLLOWING 32 PLOTS...
    for i, tr in enumerate(['201',  '021C', '021D', '210', 
                            '120U', '030C', '003', '300', 
                            '012', '021U', '120D', '102', 
                            '111U', '030T', '120C', '111D']):
        plotting_graph(optimal_temps_ok, 
                       triads[:,1][:,i], 
                       "Temperature (°C)",
                       "Triad fraction "+tr+" for substrate-product graph", 0.003)
        
    for i, tr in enumerate(['201',  '021C', '021D', '210', 
                            '120U', '030C', '003', '300', 
                            '012', '021U', '120D', '102', 
                            '111U', '030T', '120C', '111D']):
        plotting_graph(optimal_temps_ok, 
                       triads[:,0][:,i], 
                       "Temperature (°C)", 
                       "Triad fraction "+tr+" for reaction graph", 0.003)
    """

    # ENZS/NODES
    plotting_pathway_non_pathway_graph(optimal_temps_ok, 
                                       nb_enzs_all.astype(float)/nodes[:,1], 
                                       nb_enzs_path.astype(float)/nodes[:,3], 
                                       "Temperature (°C)", "#enzymes/#nodes for substrate-product graph", 
                                       ylim = [0.5, 1.2])
    plotting_pathway_non_pathway_graph(optimal_temps_ok, 
                                       nodes[:,1]/nb_enzs_all.astype(float), 
                                       nodes[:,3]/nb_enzs_path.astype(float), 
                                       "Temperature (°C)", 
                                       "#nodes/#enzymes for substrate-product graph", 
                                       ylim = [0.8, 2])
    
    

    

