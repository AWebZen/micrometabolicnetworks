#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 13:20:53 2019

@author: adele
"""

from utils import *

"""
This will create a Cronobacter_condimenti directory in the current directory where
species gene files and networks (in graphml and cpkl) will be found. 
Other directories will be created:
- genomes_cdna/ where the fasta file will be stored.
- EC_test/, where the KEGG EC files will be stored.
"""

# =============================================================================
#  Building metabolic network for species Cronobacter condimenti
# =============================================================================

"""
 To find a cDNA fasta file that works, we will try to find it through the Assembly ID,
 or Taxonomy ID if it doesn't work, or species name if the previous don't work
 You can otherwise choose one single method, functions : 
 find_fasta_with_assemblyID(), find_fasta_with_taxID(), find_fasta_with_species_name()
 """
dmu = get_fasta("dmu", "Desulfurococcus mucosus", 
                 work_dir="./") #create species directory in current directory


"""
get_substrate_product_graph will download species gene files and EC files, then
build the substrate-product graph. If after that you also want to build another
graph, it is recommended to use build_reaction_graph() or build_substrate_product_graph(),
as that will only build the graphs and avoid the time consuming steps of getting
the genes and EC files (though they will be quicker once already downloaded).
"""
dmu.get_substrate_product_graph(filtr = True, #filter ubiquitous metabolites
                                 save = True, #save graph as graphml and cpkl
                                 pathways = False, #all enzymes, not enzymes in pathways only
                                 gname="substrate_product_dmu.graphml", 
                                 pklname="substrate_product_graph_dmu.cpkl",
                                 dir_ec="EC_test/", #directory where EC files will be saved
                                 )


dmu.build_reaction_graph(gname="reaction_graph_dmu.graphml", 
                       pklname="reaction_graph_dmu.cpkl",
                       )


#Examples of network structural properties
dmu.substrate_product_graph.degree_distribution()
print "\nHierarchy flow value:", dmu.reaction_graph.hierarchy_flow()



