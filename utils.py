#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: adele
"""

from StringIO import StringIO
import logging
import os
import cPickle as cpk
import collections

from bioservices import KEGG
from Bio import SeqIO 
from Bio.KEGG import Enzyme
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress, kruskal, bartlett, f_oneway
import numpy as np
from scipy.sparse import csgraph

#Setting logging preferences
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GraphClass:
    """
    Class for the different attributes and properties of the graphs from the MetabolicGraph class
    """
    def __init__ (self, name) :
        """
        name - name/type of graph, to appear in plot titles, etc
        """
        self.graph = nx.DiGraph()
        self.node_name_equivalence = {} #global name equivalence dictionary, node ids as keys and names as values
        self.name = name
        
    def nodes(self):
        """Shortcut for node list"""
        return self.graph.nodes
    
    
    def degree_distribution(self) :
        """
        Plots the in/out/global degree distribution of our graph
        """
        in_deg = collections.Counter(dict(self.graph.in_degree).values())
        out_deg = collections.Counter(dict(self.graph.out_degree).values())
        deg = collections.Counter(dict(self.graph.degree).values())
        plt.figure()
        plt.plot(in_deg.keys(), in_deg.values(), label = "in degree", marker = "o", c = "orange")
        plt.plot(out_deg.keys(), out_deg.values(), label = "out degree", marker = "o", c = "blue")
        plt.plot(deg.keys(), deg.values(), label = "global degree", marker = "o", c = "green")
        plt.legend(loc = "best")
        plt.title("Degree distribution for the "+self.name)
        plt.xlabel("Degree")
        plt.ylabel("Frequence")
        plt.show()
        
        
    def place_compound_distribution(self, cpd):
        """
        Placing a given compound in degree distribution.
        
        cpd - node ID
        """
        in_deg = collections.Counter(dict(self.graph.in_degree).values())
        out_deg = collections.Counter(dict(self.graph.out_degree).values())
        deg = collections.Counter(dict(self.graph.degree).values())
        plt.figure()
        plt.plot(in_deg.keys(), in_deg.values(), label = "in degree", marker = "o", c = "orange", ms = 5)
        plt.plot(out_deg.keys(), out_deg.values(), label = "out degree", marker = "o", c = "blue", ms = 5)
        plt.plot(deg.keys(), deg.values(), label = "global degree", marker = "o", c = "green", ms = 5)
        plt.plot(self.graph.in_degree[cpd], in_deg[self.graph.in_degree[cpd]], label = "in "+cpd, marker = "o", c = "tomato", ms = 10)
        plt.plot(self.graph.out_degree[cpd], out_deg[self.graph.out_degree[cpd]], label = "out "+cpd, marker = "o", c = "red", ms = 10)
        plt.plot(self.graph.degree[cpd], deg[self.graph.degree[cpd]], label = "global "+cpd, marker = "o", c = "darkred", ms = 10)
        plt.legend(loc = "best")
        plt.title("Degree distribution for the "+self.name+"\n with "+ cpd+" in red")
        plt.xlabel("Degree")
        plt.ylabel("Frequence")
        plt.show()        
        
        
    def remove_nodes(self, perc) :
        """
        Removes random perc % (freq, so 0 <= perc <= 1) nodes from graph. (check
        robustness)
        
        perc - proportion of random nodes to remove from graph
        """
        assert perc <= 1 and perc >= 0, "Wrong value for percentage"
        rmv_nodes = np.random.choice(self.graph.nodes(), int(perc*len(self.graph.nodes())), replace = False)
        self.graph.remove_nodes_from(rmv_nodes)
        
        
    def remove_ecs(self, perc) :
        """
        Removes random perc % (freq, so 0 <= perc <= 1) ECs from graph. (check
        robustness), which is +/- equivalent to a random % of edges for substrate 
        graph.
        
        perc - proportion of random edges to remove from graph
        """
        assert perc <= 1 and perc >= 0, "Wrong value for percentage"
        enz_nodes = [nod for nod in self.nodes() if is_enz(nod)]
        assert len(enz_nodes) > 0, "No enzyme codes in nodes! Is this a reaction graph?"
        rmv_edges = np.random.choice(enz_nodes, 
                                     int(perc*len(enz_nodes)), 
                                     replace = False)
        self.graph.remove_nodes_from(rmv_edges)
        
        
    def mean_degree(self) :
        """
        Returns mean in degree, mean out degree and mean global degree
        """
        mean_in = np.mean(dict(self.graph.in_degree).values())
        mean_out = np.mean(dict(self.graph.out_degree).values())
        mean_glob = np.mean(dict(self.graph.degree).values())
        return mean_in, mean_out, mean_glob
    
    
    def article_density(self):
        """
        Edge density from Takemoto's 2007 and Tamura et al's 2009 article :
        D = E/N, (half of the mean degree since the sum of all degrees count twice
        every edge).
        """
        dens = self.graph.number_of_edges()/float(self.graph.number_of_nodes())
        return dens
    
    
    def eccentricity (self) :
        """
        Gives eccentricity dictionary. Does not consider the need of multiple substrates 
        to have a product. Not networkx function to ignore the fact that not all nodes are 
        connected.
        """
        lngth = nx.shortest_path_length(self.graph, weight='length')
        ecc = {}
        for l in lngth :
            ecc[l[0]] = np.max(l[1].values())
        return ecc
    
    
    def diameter(self) :
        """
        Gives graph longest shortest path (diameter) for graph. Does not consider
        the need of multiple substrates to have a product. Not networkx function 
        to ignore the fact that not all nodes are connected.
        """
        ecc = self.eccentricity()
        return np.max(ecc.values())
    
    
    def graph_density(self) :
        """
        Returns Networkx's density measure of the graph
        """
        return nx.density(self.graph)
        
    
    def weakly_connected_comps(self, plot=False):
        """
        Gives the number of weakly connected components found and the list of the
        number of components per subgraph.
        If plot == True, plots the boxplot of the list.
        """
        comp_length = []
        comps = []
        for comp in nx.weakly_connected_components(self.graph):
            if len(comp) > 1 :
                comp_length.append(len(comp))
                comps.append(comp)
        comps.sort(key = len, reverse = True)
        if plot:
            plt.figure()
            plt.boxplot(comp_length)
            plt.title("Boxplot of the number of components\nper the weakly connected components subgraph")
            plt.show()
        return (len(comp_length), sorted(comp_length, reverse = True), comps)
    
    
    def strongly_connected_comps(self, plot=False):
        """
        Gives the number of strongly connected components found and the list of the
        number of components per subgraph and the components.
        If plot == True, plots the boxplot of the list.
        """
        comp_length = []
        comps = []
        for comp in nx.strongly_connected_components(self.graph):
            if len(comp) > 1 :
                comp_length.append(len(comp))
                comps.append(comp)
        comps.sort(key = len, reverse = True)
        if plot:
            plt.figure()
            plt.boxplot(comp_length)
            plt.title("Boxplot of the number of components\nper strongly connected components subgraph")
            plt.show()
        return (len(comp_length), sorted(comp_length, reverse = True), comps)

    
    def scc_wcc_subgraph(self, suffix, directory, 
                         gname = "strongest_connected_comp", 
                         gname2 = "weak_connected_comp"):
        """
        Creates subgraph of largest strongly connected component (or 2 largests 
        if they are the same size) and the subgraph(s) of the weakly connected
        component(s) containing it (them).
        
        INPUT:
        suffix - suffix for graphml output name (organism code name for example)
        directory - directory where the graphml files will be created
        gname - name of graph for SCC, defaults to "strongest_connected_comp".
        gname2 - name of graph for WCC, defaults to "weak_connected_comp".
        """
        #Largest strongly connected component(s) subgraph(s)
        strong_graphs = sorted(nx.strongly_connected_component_subgraphs(self.graph), 
                               key=len, reverse=True)
        if len(strong_graphs[0]) == len(strong_graphs[1]): #In case first 2 SCC are of the same size, we save both
            if len(strong_graphs[0]) == len(strong_graphs[2]) :
                logger.error("At least 3 first strongly connected components of same size! - Not implemented/supported")
                raise SystemExit()
            logger.info("Saving 2 largest strongly connected components in "
                        + directory + " as " + gname + "_" + suffix 
                        + "_1.graphml and " + gname + "_" + suffix + "_2.graphml")
            nx.write_graphml(strong_graphs[0], directory + "/"+ gname + "_" + suffix + "_1.graphml")
            nx.write_graphml(strong_graphs[1], directory + "/"+ gname + "_" + suffix + "_2.graphml")
            strongest_comps = np.array([-1,-1])
        else:
            logger.info("Saving largest strongly connected components in "
                        + directory + " as "+gname + "_" + suffix + ".graphml")
            nx.write_graphml(strong_graphs[0], directory + "/"+ gname + "_" + suffix + "_1.graphml")
            strongest_comps = np.array([-1])
            
        #Weakly connected component subgraph(s) containing largest strongly connected component(s)
        weak_graphs = sorted(nx.weakly_connected_component_subgraphs(self.graph), key=len, reverse=True)
        for i, weak in enumerate(weak_graphs) :
            if len(set(strong_graphs[0].nodes()) & set(weak.nodes())) == len(strong_graphs[0].nodes()) :
                strongest_comps[0] = i #index of weakly connected component containing this strongly connected comp
            if len(strongest_comps) == 2 :
                if len(set(strong_graphs[1].nodes()) & set(weak.nodes())) == len(strong_graphs[1].nodes()) :
                    strongest_comps[1] = i
                    
            if np.any(strongest_comps == i) :
                if len(strongest_comps) == 2 and strongest_comps[0] == strongest_comps[1]:
                    logger.info("Saving weakly connected component containing 2 largest strongly connected component...")
                    nx.write_graphml(weak, directory + "/"+ gname2 + "_" + suffix + "_1and2.graphml")
                    break
                else:
                    ind = np.where(strongest_comps == i)[0]
                    logger.info("Saving weakly connected component containing largest strongly connected component...")
                    nx.write_graphml(weak, directory + "/"+ gname2 + "_" + suffix + "_" + str(ind[0]+1) + ".graphml")
                        
            if np.all(strongest_comps != -1):
                break
                        
               
    def get_full_names(self, ids):
        """
        Gives the list of names of a list of node labels from a GraphClass graph.
        Keeps KEGG compound code (CXXXXX) code as is.
        
        INPUT:
        ids - group of node names from g
        
        OUTPUT:
        names - list of names
        """
        if not np.all(np.in1d(list(ids), list(self.graph.nodes)) == True):
            logger.error("At least one wrong compound name for inputs")
            raise SystemExit()
        names = []
        for cpd in ids :
            if len(cpd) == 6 and cpd[1:6].isdigit():
                names.append(cpd)
            elif len(cpd.split(".")) == 4 and np.all(np.array([sp.isdigit() for sp in cpd.split(".")]) == True):
                names.append(cpd)
            else:
                names.append(self.node_name_equivalence[cpd])
        return names
    
    
    def centrality(self, topX, centrality = nx.closeness_centrality):
        """
        Gives topX most closeness central nodes.
        
        INPUT:
        topX - number of top nodes in output
        centrality - NetworkX's closeness centrality or betweenness... Defaults to closeness
        
        OUTPUT:
        topX_items - list of arrays node label - closeness value
        """
        cl = centrality(self.graph)
        items = np.vstack((np.array(self.get_full_names(cl.keys())), np.array(cl.values()))).T
        topX_items = sorted(items, key = lambda x: float(x[1]), reverse = True)[:topX]
        return topX_items

    
    def hierarchy_flow(self) :
        """
        Node fraction (instead of number of edges) not contained in cycle (SCC).
        We consider it to be the number of nodes in the WCC that contains SCC that 
        are not in the SCC.
        """
        weak_comps = self.weakly_connected_comps()[2]
        strong_comps = self.strongly_connected_comps()[2]
        for w in weak_comps :
            if len(strong_comps[0] & w) == len(strong_comps[0]) :
                hier_flow = (len(w) - len(strong_comps[0]))/float(len(w))
                break
        return hier_flow
    

    def high_low_in_degree(self, n) :
        """
        Gives a list of n highest in-degreed compounds and a list of n lowest in-degreed compounds 
        (with degree-value associated).
        
        n - number of components wanted.
        """
        high = sorted(list(self.graph.in_degree), key=lambda x: x[1], reverse = True)[:n]
        low = sorted(list(self.graph.in_degree), key=lambda x: x[1], reverse = True)[-n:]
        return high, low
        
    
    def high_low_out_degree(self, n) :
        """
        Gives a list of n highest out-degreed compounds and a list of n lowest out-degreed compounds 
        (with degree-value associated).
        
        n - number of components wanted.
        """
        high = sorted(list(self.graph.out_degree), key=lambda x: x[1], reverse = True)[:n]
        low = sorted(list(self.graph.out_degree), key=lambda x: x[1], reverse = True)[-n:]
        return high, low
    
    
    def high_low_degree(self, n) :
        """
        Gives a list of n highest degreed compounds and a list of n lowest degreed compounds 
        (with degree-value associated).
        
        n - number of components wanted.
        """
        high = sorted(list(self.graph.degree), key=lambda x: x[1], reverse = True)[:n]
        low = sorted(list(self.graph.degree), key=lambda x: x[1], reverse = True)[-n:]
        return high, low
            
    
    def load_graph_graphml(self, fname):
        """
        Loads graph from existing graphml file and node name equivalence file.
        
        INPUT:
        fname - graph file path
        fname_equiv - input file path (cpickle or csv) for node name equivalence file.
        cpkl -  boolean if cpickle file for node name equivalence file.
        """
        if os.path.exists(fname) :
            self.graph = nx.read_graphml(fname)
            self.node_name_equivalence = nx.get_node_attributes(self.graph, "id")
        else:
            logger.error("File {} does not exist!".format(fname))
        
        
    def load_graph_cpkl(self, fname):
        """
        Loads graph from existing cPickle binary file and node name equivalence file.
        
        INPUT:
        fname - file path
        fname_equiv - input file path (cpickle or csv) for node name equivalence file.
        cpkl -  boolean if cpickle file for node name equivalence file.
        """
        if os.path.exists(fname) :
            self.graph = cpk.load(open(fname, "rb"))
            self.node_name_equivalence = nx.get_node_attributes(self.graph, "id")
        else:
            logger.error("File {} does not exist!".format(fname))
            
    def merge_graphs(self, list_graphs) :
        """
        Brute forcing the merge of a list of graphs.
        
        Problems of label incongruencies between labels not handled (same nodes
        different labels, different nodes same labels)... Hope this doesn't happen.
        
        list_graphs - list of NetworkX graphs.
        """
        self.graph = nx.compose_all(list(list_graphs))
        self.node_name_equivalence = nx.get_node_attributes(self.graph, "id")
        #Handling same compounds with different node labels
        reverse_equiv = dict(zip(self.node_name_equivalence.values(), self.node_name_equivalence.keys()))
        if len(self.node_name_equivalence.keys()) != len(reverse_equiv.keys()) :
#            print len(self.node_name_equivalence.keys()), len(reverse_equiv.keys())
            logger.info("Merging duplicates...")
            c = collections.Counter(self.node_name_equivalence.values())
            duplicates = map(lambda x: x[0], filter(lambda x: x[1] >= 2, c.items()))
            items = np.array(self.node_name_equivalence.items())
            for dup in duplicates :
                indexes = np.where(items[:,1] == dup)[0]
                nodes = items[indexes, 0]
                for i in xrange(1,len(nodes)) :
                    self.graph = nx.contracted_nodes(self.graph, nodes[0], nodes[i], self_loops=False)
                self.node_name_equivalence = nx.get_node_attributes(self.graph, "id")
            if len(self.node_name_equivalence.keys()) != len(reverse_equiv.keys()) :
#                print len(self.node_name_equivalence.keys()), len(reverse_equiv.keys())
                logger.error("Uh oh... something wrong when handling label incongruencies...")
                raise SystemError()
            
                    
                    



class MetabolicGraph:
    """
    Builds metabolic graphs thanks to the KEGG database for a given species.
    """
    k = KEGG() #Bioservices' KEGG interface
    k.settings.TIMEOUT = 1000 #Changing timeout
    ubi_metab = ["C00001", "C00002", "C00008", "C00003", "C00004", "C00005", 
                     "C00006",  "C00011",  "C00014", "C00059", "C00342", "C00009", 
                     "C00013", "C00080"] #C00006 - NADP added
    
    def __init__(self, organism_name, fasta_name, 
                 code="", KO=False, merged=False, work_dir="./"):
        """
        To build hybrid species, we need a list of fasta files instead of just
        the name of one, and a list of codes instead of just one.
        
        organism_name - organism name given as input, space separated (or list)
        fasta_name - file name/path for fasta file. Needs to be a cDNA fasta file with a "gene:", "gene="
                or "locus_tag=" field in the description of each sequence, such as the ones from Ensembl DB.
                Can be a list of fasta names for artificial hybrid species.
        code - if KEGG organism 3-letter code is already known, it can be given.
                Else it will be deduced from organism name.
                If hybrid species, must be a list of already known codes.
        KO - boolean. If we build graphs through KO or not. Defaults to False.
            Instead of True, can put a KO dictionary with KOs as keys and ECs
            as values to quicken code.
        merged - boolean if we want to create an object of merged graphs.
        work_dir - working directory where the species directory will be located. Defaults to current dir.
        """
        if type(organism_name) == list:
            organism_name = "_".join(organism_name)

        self.organism_name = organism_name
        self.directory = work_dir + organism_name.replace(" ", "_") #directory - name of the directory with KEGG gene files
        self.code = code #KEGG organism code
        self.number_genes = 0 #Number of genes kept (valid)
        self.valid_ecs = [] #list of one or multiple EC codes as strings
        self.enzs_parsed = [] #List of KEGG enzyme entry file parsers
        self.genes = []
        self.multiple = [] #List of starting indexes of the genes in other species
        self.KO = KO
        
        self.reaction_graph = GraphClass("reaction graph") #reaction graph, filtered
        self.unfiltered_reaction_graph = GraphClass("unfiltered reaction graph") #unfiltered reaction graph (all metabolites)
        self.pathway_reaction_graph = GraphClass("pathway reaction graph") #reaction graph, filtered
        self.pathway_unfiltered_reaction_graph = GraphClass("pathway unfiltered reaction graph") #unfiltered reaction graph (all metabolites)
        
        self.substrate_product_graph = GraphClass("substrate-product graph") # substrate-product graph
        self.unfiltered_substrate_product_graph = GraphClass("unfiltered substrate-product graph") #unfiltered substrate-product graph (all metabolites)
        self.pathway_substrate_product_graph = GraphClass("pathway substrate-product graph") # substrate-product graph
        self.pathway_unfiltered_substrate_product_graph = GraphClass("pathway unfiltered substrate-product graph") #unfiltered substrate-product graph (all metabolites)
        
        self.in_out_graph = GraphClass("in-out graph") # in-out graph
        self.pathway_in_out_graph = GraphClass("pathway in-out graph")
        
        if not KO and not merged: #If not building graphs with KO, need organism code and fastas
            if type(fasta_name) == list: #Hybrid species
                assert type(self.code) == list and len(self.code) > 1, "Missing multiple codes as list"
                for i, fname in enumerate(fasta_name) :
                    self.multiple.append(len(self.genes))
                    self.load_data_file(fname) #Get gene names
            else:
                self.load_data_file(fasta_name) #Get gene names
                
            if self.code == "" : #Find organism code
                self.findandtest_organism()
            elif type(self.code) == str or type(self.code) == unicode or type(self.code) == np.unicode or type(self.code) == np.unicode_ :
                self.test_code() #Testing gene name - organism code correspondance (tests a single gene)
            else :
                self.test_multiple_codes()
    
    
    def merge_graphs(self, graphs_list, graph_type, filtr = True, pathway = False) :
        """
        Merge a list of metabolic graphs together.
        
        INPUT:
        graphs_list - list of networkx graphs
        graph_type - must be among "reaction", "substrate", "in-out"
        filtr - 13 ubiquitous metabolites filtered or not. Defaults to filtered.
        pathway - only enzymes in known pathways or not. Defaults to False.
        """
        logger.warning('Graph type must be among "reaction", "substrate", "in-out"')
        graph_type = graph_type.lower()
        assert graph_type in ["reaction", "substrate", "in-out"], "Wrong graph type"
        if filtr :
            if graph_type == "reaction":
                if pathway:
                    self.pathway_reaction_graph.merge_graphs(graphs_list)
                else:
                    self.reaction_graph.merge_graphs(graphs_list)
            elif graph_type == "substrate":
                if pathway:
                    self.pathway_substrate_product_graph.merge_graphs(graphs_list)
                else:
                    self.substrate_product_graph.merge_graphs(graphs_list)
            elif graph_type == "in-out":
                if pathway:
                    self.pathway_in_out_graph.merge_graphs(graphs_list)
                else:
                    self.in_out_graph.merge_graphs(graphs_list)
        else:
            if graph_type == "reaction":
                if pathway:
                    self.pathway_unfiltered_reaction_graph.merge_graphs(graphs_list)
                else:
                    self.unfiltered_reaction_graph.merge_graphs(graphs_list)
            elif graph_type == "substrate":
                if pathway:
                    self.pathway_unfiltered_substrate_product_graph.merge_graphs(graphs_list)
                else:
                    self.unfiltered_substrate_product_graph.merge_graphs(graphs_list)
        
    
    def load_data_file(self, fname):
        """ 
        Loads fasta file using Bio parser and returns genes names
        
        INPUT:
        fname - file name/path for fasta file. Needs to be a cDNA fasta file with a "gene:"
                field in the description of each sequence, such as the ones from Ensembl DB.
                Also supports cDNA fasta files with [locus_tag=XXXXXX] field, such as the 
                fastas from Genbank.
        """
        seqs = [s for s in SeqIO.parse(fname, "fasta")]
        self._genes = []
        for seq in seqs :
            descr = seq.description
            i_gen = descr.find("gene:") #Find gene name field
            i_gen2 = descr.find("locus_tag=") #For fasta from Genbank
            i_gen3 = descr.find("gene=") #For fasta from Genbank
            if i_gen != -1 :
                gene = descr[i_gen+5:].split()[0]
                self.genes.append(gene)
            elif i_gen2 != -1:
                gene = descr[i_gen2+10:].split("]")[0]
                self.genes.append(gene)
                if i_gen3 != -1: #Two possibilities from Genbank present
                    gene = descr[i_gen3+5:].split("]")[0]
                    self._genes.append(gene)
            elif i_gen3 != -1: #Last priority
                gene = descr[i_gen3+5:].split("]")[0]
                self.genes.append(gene)
        if len(self.genes) != len(seqs):
            if len(self.genes) <= int(.5*len(seqs)):
                logger.error("Could not find enough gene names. Is field 'gene:'/'gene='/'locus_tag=' present in your fasta descriptions?")
                raise SystemExit()
            else :
                logger.warning("Not all gene names found.")
    
    
    def test_code(self) :
        """
        Tests if 3-letter KEGG species code works with the set of gene names 
        (tests only one) from the fasta file (if there is a correspondance in KEGG).
        """
        if type(MetabolicGraph.k.get(self.code + ":" + self.genes[3])) == type(1) :
            if len(self._genes) == 0 or type(MetabolicGraph.k.get(self.code + ":" + self._genes[3])) == type(1) : #Priority to locus_tag= rather than gene=
                logger.error("Uh oh! 3-letter KEGG species code does not work with fasta file genes...")
                raise SystemExit()
            else:
                self.genes = self._genes
    
    
    def test_multiple_codes(self) :
        """
        Tests if 3-letter KEGG species code works with the set of gene names 
        (tests only one) from the fasta file (if there is a correspondance in KEGG).
        """
        for i, code in enumerate(self.code) :
            if type(MetabolicGraph.k.get(code + ":" + self.genes[self.multiple[i]])) == type(1) :
                logger.error("Uh oh! 3-letter KEGG species code does not work with fasta file genes for species %d!" %(i+1))
                raise SystemExit()
            
        
    def get_organism(self, org_name):
        """
        Finds the KEGG organism code name through the organism name. Tests hits found.
        
        INPUT:
        org_name - name of organism or parts of it, space separated
        
        OUTPUT:
        code - KEGG organism ID/code or None, if not found
        """
        org_list = MetabolicGraph.k.lookfor_organism(org_name)
        if len(org_list) > 0: #Found possible organism hits
            for org in org_list: #Test hits
                code = org.split()[1]
                txt = MetabolicGraph.k.get(code + ":" + self.genes[3])
                try: #If organism code works, keep code
                    int(txt)
                except ValueError:
                    self.code = code
                    return code
        return None
    
    
    def findandtest_organism(self, work_dir = "./"):
        """
        Finds the KEGG organism code name through the organism name. If not found, tests
        with parts of the name as query. If not found, asks the user for a new name.
        Raises an error if no code name found at the end.
        
        INPUT:
        work_dir - working directory where the species directory will be located. 
                   Defaults to current dir.
        
        OUTPUT:
        code - KEGG organism ID/code or Error, if not found
        """
        logger.info("Looking for organism code in KEGG...")
        code = self.get_organism(self.organism_name)
        if code == None:
            org_name_list = self.organism_name.split()
            org_name_list.append(org_name_list.pop(0)) #reshuffling : putting the first element of the name (usually the genus) as the last one to test, as it probably has a lot more hits
            logger.info("No hits for whole organism name, testing with parts of name...")
            for name in org_name_list: #Test parts of name
                code = self.get_organism(name)
                if code != None:
                    break
            if code == None:
                new_name = raw_input("Organism name " + self.organism_name + " was not found in KEGG, write another name for it (enter S to stop) : ")
                if new_name.lower() == "s" :
                    logger.error("Uh oh! Organism name not found in KEGG database!")
                    raise SystemExit()
                else:
                    self.organism_name = new_name #Updating attributes
                    self.directory = work_dir + self.organism_name.replace(" ", "_")
                    code = self.findandtest_organism()
        if code != None:
            self.code = code
            logger.info("Organism code found!")
            
    
    def get_kegg_genes(self) :
        """
        Downloads KEGG gene files into org_name directory.
        """
        logger.info("Fetching KEGG gene entries...")
        count = 0
        if type(self.code) == str or type(self.code) == unicode or type(self.code) == np.unicode or type(self.code) == np.unicode_ :
            code = [self.code]
        i_cod = 0
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        for i_gen, gene in enumerate(self.genes) :
            if not os.path.exists(self.directory + "/" + gene + "_gene.txt") : #Download missing gene files
                if i_gen in self.multiple :
                    i_cod = self.multiple.index(i_gen)
                txt = MetabolicGraph.k.get(code[i_cod] + ":" + gene)
                try:
                    int(txt)
                    count += 1
                    if count > 0.5 * len(self.genes) :
                        break
                except ValueError:
                    open(self.directory + "/" + gene + "_gene.txt", "w").write(txt)
        if count > 0.5 * len(self.genes) :
            logger.error("Not enough gene hits in KEGG database!")
            raise SystemExit()
        elif count != 0:
            logger.warning("No hits in the gene KEGG database for %d genes" %count)
        self.number_genes = len(self.genes) - count
        
    
    def get_ECs(self, dir_ec):
        """
        Extracts ECs for all gene files in our directory and downloads KEGG enzyme
        entries in directory/EC_files/ directory.
        """
        
        def extract_ec_number(fname):
            """
            Extracts EC and KO number(s) (if found) from the orthology field of a KEGG entry for a given gene.
        
            INPUT:
            fname - name/path of a KEGG gene database file (downloaded with get_kegg_genes() or otherwise)
            
            OUTPUT:
            KO - list of KO ids retrieved in ORTHOLOGY field
            ec_all - string of space-separated EC numbers, empty string otherwise
            """
            iOF = open(fname, "r")
            line = iOF.readline()
            KO = []
            ec_all = ""
            while line != "" and not line.startswith("ORTHO"):
                line = iOF.readline() #Skip lines until ORTHOLOGY field reached
            if line.startswith("ORTHO"): #ORTHOLOGY field exists
                while line[0] == " " or line[0:5] == "ORTHO":  #while still in ORTHOLOGY field
                    line = line.lstrip("ORTHOLOGY ") #Any of these characters are stripped from the beginning of str, order does not matter
                    ll = line.split()
                    if ll[0][1:].isdigit() and line[0] == "K" and len(ll[0]) == 6 : #Check presence of KO id
                         KO.append(ll[0])
                    i_ec = line.find("EC")
                    if i_ec != -1: #There should be at least one EC
                        ec = line[i_ec+3:].split("]")[0] #Cropping first 3 characters ("EC:") and last ] of [EC:XXXXXXX] field
                        ECs = ec.split() #List of ECs
                        for EC in ECs:
                            if EC.find(".") != -1 : #EC confirmed
                                if EC not in ec_all :
                                    ec_all += " " + EC
                    line = iOF.readline()
            iOF.close()
            return KO, ec_all
        
        
        logger.info("Fetching KEGG enzyme entries...")

        all_ECs = [] #List of ECs with hits in KEGG db
        gene_files = os.listdir(self.directory)
        if not os.path.exists(dir_ec):
            logger.info("Creating EC directory")
            os.makedirs(dir_ec)
        if not os.path.exists(self.directory) or not os.path.exists(self.directory + "/backups/"):
            os.makedirs(self.directory + "/backups/")
            
        #Check if shortcut exists (if user has already run function once, and EC list has been saved)
        if os.path.exists(self.directory + "/backups/valid_EC_list.cpkl") :
            all_ECs = cpk.load(open(self.directory + "/backups/valid_EC_list.cpkl", "rb"))
            for valid_ecs in all_ECs: #valid ECs taken from one of the gene files
                name = valid_ecs[1:].replace(" ", "_")
                if not os.path.exists(dir_ec + "/ec_" + name + ".txt") : #get missing EC files
                    txt = MetabolicGraph.k.get("ec:" + valid_ecs)
                    try:
                        int(txt)
                    except ValueError:
                        open(dir_ec + "/ec_" + name + ".txt", "w").write(txt)
        else: #Complete download. Possible improvement (?) : with bash check if number of ec_ files in EC_files/ is the same as 'grep -l "EC:" '+ self.directory + '/*|wc' ??
            for fi in gene_files :
                if fi.endswith("_gene.txt") :
                    fname = self.directory + "/" + fi
                    KO, ECs = extract_ec_number(fname)
                    if len(ECs) > 0 :
                        name = ECs[1:].replace(" ", "_") #We don't gain much time since we parse every gene file...
                        if not os.path.exists(dir_ec + "/ec_" + name + ".txt") : #If not first time dowloading, will check only non valid ECs 
                            txt = MetabolicGraph.k.get("ec:" + ECs)
                            try:
                                int(txt)
                            except ValueError:
                                all_ECs.append(ECs)
                                open(dir_ec + "/ec_" + name + ".txt", "w").write(txt)
                        else:
                            if ECs not in all_ECs :
                                all_ECs.append(ECs)
            cpk.dump(all_ECs, open(self.directory + "/backups/valid_EC_list.cpkl", "wb"))
        self.valid_ecs = all_ECs
        self.dir_ec = dir_ec
        
    
    def get_ecs_from_KOs(self, KO_list, dir_ec):
        def extract_ec_number_KO(ko, ko_dict):
            """
            Extracts EC (if found) from the definition field of a KEGG entry for a KO.
        
            INPUT:
            ko - Kegg Orthology (KO) code name, string
            ko_dict - boolean or dict of KO keys and their associated ECs as values
            
            OUTPUT:
            ec_all - string of space-separated EC numbers, empty string otherwise
            """
            try:
                if ko in ko_dict.keys():
                    return ko_dict[ko]
            except TypeError :
                pass
            txt = k.get("ko:"+ko)
            try :
                int(txt)
                return ""
            except ValueError:
                txt = txt.split("\n")
                ec_all = ""
                i = 0
                line = txt[i]
                while line != "" and not line.startswith("DEFINITION"):
                    i += 1
                    line = txt[i] #Skip lines until DEFINITION field reached
                if line.startswith("DEFINITION"): #DEFINITION field exists
                    while line[0] == " " or line[0:5] == "DEFIN":  #while still in DEFINITION field
                        line = line.lstrip("DEFINITION ") #Any of these characters are stripped from the beginning of str, order does not matter
                        i_ec = line.find("EC:")
                        if i_ec != -1: #There should be at least one EC
                            ec = line[i_ec+3:].split("]")[0] #Cropping first 3 characters ("EC:") and last ] of [EC:XXXXXXX] field
                            ECs = ec.split() #List of ECs
                            for EC in ECs:
                                if EC.find(".") != -1 : #EC confirmed
                                    if EC not in ec_all :
                                        ec_all += " " + EC
                        i += 1
                        line = txt[i]
            return ec_all
        
        
        logger.info("Fetching KEGG enzyme entries...")
        all_ECs = [] #List of ECs with hits in KEGG db

        if not os.path.exists(dir_ec):
            logger.error("{} directory given in command does not exist! Check path (current one: {})".format(dir_ec, os.getcwd()))
            raise SystemExit()
        if not os.path.exists(self.directory) or not os.path.exists(self.directory + "/backups/"):
            os.makedirs(self.directory + "/backups/")
        
        #Check if shortcut exists (if user has already run function once, and EC list has been saved)
        if os.path.exists(self.directory + "/backups/valid_EC_list.cpkl") :
            all_ECs = cpk.load(open(self.directory + "/backups/valid_EC_list.cpkl", "rb"))
            for valid_ecs in all_ECs: #valid ECs taken from one of the gene files
                name = valid_ecs[1:].replace(" ", "_")
                if not os.path.exists(dir_ec + "/ec_" + name + ".txt") : #get missing EC files in global EC directory
                    logger.info("Fetching undownloaded EC files: {}".format(valid_ecs))
                    txt = MetabolicGraph.k.get("ec:" + valid_ecs)
                    try:
                        int(txt)
                    except ValueError:
                        open(dir_ec + "/ec_" + name + ".txt", "w").write(txt)
        
        else: #Complete download
            for ko in KO_list :
                ECs = extract_ec_number_KO(ko, self.KO)
                if len(ECs) > 0 :
                    name = ECs[1:].replace(" ", "_") #We don't gain much time since we parse every gene file...
                    if not os.path.exists(dir_ec + "/ec_" + name + ".txt") : #If not first time dowloading, will check only non valid ECs 
                        txt = MetabolicGraph.k.get("ec:" + ECs)
                        try:
                            int(txt)
                        except ValueError:
                            all_ECs.append(ECs)
                            open(dir_ec + "/ec_" + name + ".txt", "w").write(txt)
                    else:
                        if ECs not in all_ECs :
                            all_ECs.append(ECs)
            cpk.dump(all_ECs, open(self.directory + "/backups/valid_EC_list.cpkl", "wb"))
        self.valid_ecs = all_ECs
        self.dir_ec = dir_ec
        
    
    def parse_enzymes(self) :
        """
        Retrieves all KEGG enzyme records with Biopython parser. Saves them as cpickle
        object for backup.
        
        OUTPUT :
        enzs - list of enzyme records
        """
        enzs = []
        logger.info("Parsing enzymes...")
        
        
        if os.path.exists(self.directory+"/EC_files"):
            if os.path.exists(self.directory + "/EC_files/enzs_parser_backup.cpkl"): #Gains only a few seconds...
                enzs = cpk.load(open(self.directory + "/EC_files/enzs_parser_backup.cpkl", "rb"))
            else:
                for fi in sorted(os.listdir(self.directory + "/EC_files")):
                    if fi.startswith("ec_"):
                        enzs += list(Enzyme.parse(open(self.directory+"/EC_files/"+fi)))
                cpk.dump(enzs, open(self.directory + "/EC_files/enzs_parser_backup.cpkl", "wb"))
                
        else:
            try:
                if not os.path.exists(self.dir_ec):
                    logger.error("<{}> global EC directory does not exist! Check path (current one: {})".format(self.dir_ec, os.getcwd()))
                    raise SystemExit()
            except (NameError, AttributeError):
                logger.error("self.dir_ec does not exist. Run get_ecs_from_KOs?")
                raise SystemExit()
            if not self.valid_ecs and not os.path.exists(self.directory + "/backups/valid_EC_list.cpkl"):
                logger.error("Run get_ecs_from_KOs")
                raise SystemExit()
            if os.path.exists(self.directory + "/backups/enzs_parser_backup.cpkl"): #Gains only a few seconds...
                enzs = cpk.load(open(self.directory + "/backups/enzs_parser_backup.cpkl", "rb"))
            else:
                for ecs in sorted(self.valid_ecs):
                    name = ecs[1:].replace(" ", "_")
                    fi = self.dir_ec + "/ec_" + name + ".txt"
                    try:
                        enzs += list(Enzyme.parse(open(fi)))
                    except IOError:
                        logger.error("<{}> file does not exist".format(fi))
                        raise SystemExit()
                cpk.dump(enzs, open(self.directory + "/backups/enzs_parser_backup.cpkl", "wb"))
        return enzs
        
        
    def get_substrates_products(self, e, filtr, graphe):
        """
        Finds unique substrate and products node ids and updates name equivalence dictionary.
        May filter following compounds : water, ATP, ADP, NAD, NADH, NADPH, carbon dioxide, 
        ammonia, sulfate, thioredoxin, (ortho) phosphate (P), pyrophosphate (PPi), H+ and NADP.
        
        Will consider as different compounds the metabolites that also appear in compounds that are
        actually a list, or slightly different name versions of same metabolite.
        
        INPUT: 
        e - KEGG enzyme/reaction entry parser (Biopython)
        filtr - boolean. If True, filters list of ubiquitous metabolites.
        graphe - determines to which graph these compounds need to be added
        
        OUPUT:
        subs - list of substrate node ids for given reaction, each being 10-char long
        prod - list of product node ids for given reaction, each being 10-char long
        """
        
        def extract_compound(comp) :
            """
            Extracts compound code or first 10 characters if code is not present.
            
            INPUT:
            comp - string of compound
            
            OUTPUT:
            compound code or 10 first compound characters
            i_cpd - -1 when no compound code
            """
            i_cpd = comp.find('CPD:')
            if i_cpd == -1:
                return comp[:10].upper(), i_cpd #+/- random 10-char code
            else:
                return comp[i_cpd+4:].split("]")[0], i_cpd #CPD code
    
        
        ubi_metab = ["C00001", "C00002", "C00008", "C00003", "C00004", "C00005", 
                     "C00006",  "C00011",  "C00014", "C00059", "C00342", "C00009", 
                     "C00013", "C00080"] #C00006 - NADP added
        subs = [] #Substrate node ids
        prod = [] #Product node ids
        for s in e.substrate :
            sub, i_cpd = extract_compound(s)
            if filtr :
                if sub in ubi_metab :
                    continue
            if s not in graphe.node_name_equivalence.values(): #Check if substrate exists in our equivalence dictionary
                i = 0
                while sub in graphe.node_name_equivalence.keys() and i_cpd == -1 : #Check if by bad luck our random compound node id exists in dictionary. Compound code should be unique.
                    if s[i*10+10:] != "" :
                        sub, i_cpd = extract_compound(s[i*10+10:]) #Find new compound node id in name
                    else :
                        sub += str(i) #add number if no unique compound node id can be found
                    i += 1
                graphe.node_name_equivalence[sub] = s
            else:
                sub = [k for k,name in graphe.node_name_equivalence.items() if name == s][0]
            subs.append(sub)
        for p in e.product :
            prd, i_cpd = extract_compound(p)
            if filtr :
                if prd in ubi_metab :
                    continue
            if p not in graphe.node_name_equivalence.values(): #Check if product exists in our equivalence dictionary
                i = 0
                while prd in graphe.node_name_equivalence.keys() and i_cpd == -1 : #Check if by bad luck our random compound node id exists
                    if p[i*10+10:] != "" :
                        prd, i_cpd = extract_compound(p[i*10+10:]) #Find new compound node id
                    else :
                        prd += str(i)
                    i += 1
                graphe.node_name_equivalence[prd] = p
            else:
                prd = [k for k,name in graphe.node_name_equivalence.items() if name == p][0]
            prod.append(prd)
        return subs, prod
    
    
    def build_reaction_graph (self, filtr = True, save = True, gname = "metabolites_reaction.graphml", pklname = 'metabolites_reactions_graph.cpkl', pathways = False):
        """
        Builds a directed reaction graph (substrates -> enzyme -> products).
        Skips enzymes without product and substrate entries.
        
        INPUT:
        filtr - boolean. If True, filters list of ubiquitous metabolites. Defaults to True.
        save - if True, saves graph as graphml. Defaults to True.
        gname - graph name if save = True. Defaults to "metabolites_reaction.graphml".
        pklname - cpickle graph name if save = True. Defaults to "metabolites_reactions_graph.cpkl".
        pathways - if we only want enzymes known to be in a pathway. Defaults to False.
        
        OUTPUT:
        graphe - reaction graph
        """
        if len(self.enzs_parsed) == 0 : 
            enzs = self.parse_enzymes()
            self.enzs_parsed = enzs
        else: #skips step if already built a reaction graph -> already parsed enzyme files
            enzs = self.enzs_parsed
        logger.info("Building graph...")
        count_skip = 0
        count_skip_paths = 0
        if filtr :
            if pathways :
                graphe = self.pathway_reaction_graph
            else:
                graphe = self.reaction_graph
        else:
            if pathways :
                graphe = self.pathway_unfiltered_reaction_graph
            else:
                graphe = self.unfiltered_reaction_graph
            
        for e in enzs:
            if len(e.substrate) == 0 and len(e.product) == 0: #Skipping enzymes without substrate and product entries
                count_skip += 1
                continue
            if pathways :
                if len(e.pathway) == 0:
                    count_skip_paths += 1
                    continue
            #Substrate and product names
            if filtr :
                subs, prod = self.get_substrates_products(e, True, graphe)
            else:
                subs, prod = self.get_substrates_products(e, False, graphe)
            er = e.entry #EC code
            graphe.node_name_equivalence[er] = e.name[0]
            
            #Building graph
            graphe.graph.add_node(er, id=e.name[0]) #takes the first out of all synonymous enzyme names
            for s in subs:
                graphe.graph.add_node(s, id=graphe.node_name_equivalence[s])
                graphe.graph.add_edge(s, er)
            
            for p in prod:
                graphe.graph.add_node(p, id=graphe.node_name_equivalence[p])
                graphe.graph.add_edge(er, p)
                
        #Writing output file with name equivalences
        if save:
            if filtr:
                logger.info("Saving graph as %s and as cpickle object %s in directory %s", 
                            gname, pklname, self.directory)
                nx.write_graphml(graphe.graph, self.directory + "/" + gname)
                cpk.dump(graphe.graph, open(self.directory + "/" + pklname,'wb'))
            else:
                logger.info("Saving graph as %s and as cpickle object %s in directory %s", 
                            "unfiltered_" + gname, "unfiltered_"  +pklname, self.directory)
                nx.write_graphml(graphe.graph, self.directory + "/unfiltered_"+gname)
                cpk.dump(graphe.graph, open(self.directory + "/unfiltered_"+pklname,'wb'))
        
        if count_skip > 0 :
            logger.warning("%d/%d enzyme entries without substrate nor product information have been skipped, and will not appear in graph.", count_skip, len(enzs))
        if count_skip_paths > 0:
            logger.warning("%d/%d enzyme entries without pathway have been skipped, and will not appear in graph", count_skip_paths, len(enzs))
        return graphe
    
    
    def build_substrate_product_graph (self, filtr = True, save = True,
                                       gname = "metabolites_substrate_product.graphml", 
                                       pklname = 'metabolites_substrate_product_graph.cpkl', 
                                       pathways = False):
        """
        Builds a directed substract-product graph (substrates -> products of same reaction).
        Skips enzymes without product and substrate entries.
        
        Similar graph to Takemoto et al.'s 2007 article graphs but directed.
        
        INPUT:
        filtr - boolean. If True, filters list of ubiquitous metabolites. Defaults to True.
        save - if True, saves graph as graphml. Defaults to True.
        gname - graph name if save = True. Defaults to "metabolites_substrate_product.graphml".
        pklname - cpickle graph name if save = True. Defaults to "metabolites_substrate_product_graph.cpkl".
        pathways - if we only want enzymes known to be in a pathway. Defaults to False.
        
        OUTPUT:
        graphe - substrate product graph
        """
        if len(self.enzs_parsed) == 0 : 
            enzs = self.parse_enzymes()
            self.enzs_parsed = enzs
        else: #skips step if already built a reaction graph -> already parsed enzyme files
            enzs = self.enzs_parsed
        logger.info("Building graph...")
        count_skip = 0
        count_skip_paths = 0
        if filtr :
            if pathways:
                graphe = self.pathway_substrate_product_graph
            else:
                graphe = self.substrate_product_graph
        else:
            if pathways:
                graphe = self.pathway_unfiltered_substrate_product_graph
            else:
                graphe = self.unfiltered_substrate_product_graph
            
        for e in enzs:
            if len(e.substrate) == 0 and len(e.product) == 0: #Skipping enzymes without substrate and product entries
                count_skip += 1
                continue
            if pathways :
                if len(e.pathway) == 0:
                    count_skip_paths += 1
                    continue
            #Substrate and product names
            if filtr :
                subs, prod = self.get_substrates_products(e, True, graphe)
            else:
                subs, prod = self.get_substrates_products(e, False, graphe)
            #Building graph
            if len(subs) != 0:
                for s in subs:
                    graphe.graph.add_node(s, id=graphe.node_name_equivalence[s])
                    for p in prod:
                        graphe.graph.add_node(p, id=graphe.node_name_equivalence[p])
                        graphe.graph.add_edge(s, p)
            else :
                for p in prod:
                    graphe.graph.add_node(p, id=graphe.node_name_equivalence[p])
                
        #Writing output file with name equivalences
        if save:
            if filtr:    
                logger.info("Saving graph as %s and as cpickle object %s in directory %s", 
                            gname, pklname, self.directory)
                nx.write_graphml(graphe.graph, self.directory + "/" + gname)
                cpk.dump(graphe.graph, open(self.directory + "/" + pklname,'wb'))
            else:
                logger.info("Saving graph as %s and as cpickle object %s in directory %s", 
                            "unfiltered_"+gname, "unfiltered_"+pklname, self.directory)
                nx.write_graphml(graphe.graph, self.directory + "/unfiltered_"+gname)
                cpk.dump(graphe.graph, open(self.directory + "/unfiltered_"+pklname,'wb'))
        
        if count_skip > 0 :
            logger.warning("%d/%d enzyme entries without substrate nor product information have been skipped, and will not appear in graph.", count_skip, len(enzs))
        if count_skip_paths > 0:
            logger.warning("%d/%d enzyme entries without pathway have been skipped, and will not appear in graph", count_skip_paths, len(enzs))
        return graphe


    def get_reaction_graph(self, filtr = True, save = True, 
                           gname="metabolites_reaction.graphml", 
                           pklname='metabolites_reactions_graph.cpkl', 
                           pathways=False, ko_list=[], 
                           dir_ec=""):
        """
        Global function for building a reaction graph, if you don't want to do every 
        step of it (fetching gene entries, fetching enzyme entries, building the graph).
        Once a graph has already been built (gene and enzyme entries already fetched), it 
        is recommended to use build_reaction_graph() or build_substrate_product_graph() directly.
        
        INPUT:
        filtr - boolean. If True, filters list of ubiquitous metabolites. Defaults to True.
        save - if True, saves graph as graphml. Defaults to True.
        gname - graph name if save = True. Defaults to "metabolites_reaction.graphml".
        pklname - cpickle graph name if save = True. Defaults to "metabolites_reactions_graph.cpkl".
        pathways - if we only want enzymes known to be in a pathway. Defaults to False.
        ko_list - list of KOs for when self.KO != False
        dir_ec - global directory where ECs have already been downloaded
        
        OUTPUT :
        graphe - self.reaction_graph or self.unfiltered_reaction_graph
        """
#        if len(self.reaction_graph.node_name_equivalence.keys()) == 0 and len(self.unfiltered_reaction_graph.node_name_equivalence.keys()) == 0 :
        
        if not self.KO :
            self.get_kegg_genes()
            self.get_ECs(dir_ec)
        else :
            assert len(ko_list) > 0, "ko_list argument must have at least one KO code"
            self.get_ecs_from_KOs(ko_list, dir_ec)
            

        graphe = self.build_reaction_graph (filtr, save, gname, pklname, pathways)
        return graphe
    
    
    def get_substrate_product_graph(self, filtr = True, save = True, 
                            gname = "metabolites_substrate_product.graphml", 
                            pklname = 'metabolites_substrate_product_graph.cpkl', 
                            pathways = False, ko_list = [], 
                            dir_ec = ""):
        """
        Global function for building a substrate-product graph, if you don't want to
        do every step of it (fetching gene entries, fetching enzyme entries, building the graph).
        Once a graph has already been built (gene and enzyme entries already fetched), it is
        recommended to use build_reaction_graph() or build_substrate_product_graph() directly.
        
        INPUT:
        filtr - boolean. If True, filters list of ubiquitous metabolites. Defaults to True.
        save - if True, saves graph as graphml. Defaults to True.
        gname - graph name if save = True. Defaults to "metabolites_substrate_product.graphml".
        pklname - cpickle graph name if save = True. Defaults to "metabolites_substrate_product_graph.cpkl".
        pathways - if we only want enzymes known to be in a pathway. Defaults to False.
        ko_list - list of KOs for when self.KO != False
        dir_ec - global directory where ECs have already been downloaded
        
        OUTPUT :
        graphe - self.substrate_product_graph or self.unfiltered_substrate_product_graph
        """
#        if len(self.substrate_product_graph.node_name_equivalence.keys()) == 0 and len(self.unfiltered_substrate_product_graph.node_name_equivalence.keys()) == 0 :
            
        if not self.KO :
            self.get_kegg_genes()
            self.get_ECs(dir_ec)
        else :
            assert len(ko_list) > 0, "ko_list argument must have at least one KO code"
            self.get_ecs_from_KOs(ko_list, dir_ec)

        graphe = self.build_substrate_product_graph (filtr, save, gname, pklname, pathways)
        return graphe


# =============================================================================
# General functions
# =============================================================================
def is_enz (node):
    """Tells if node is an enzyme (EC) or not"""
    split = node.split(".")
    if len(split) == 4 :#and np.all(np.array([sp.isdigit() for sp in split]) == True) :
        return True
    else:
        return False


def plotting_graph (x, y, xlabel, ylabel, yless, 
                    topredictX=False, ylim=[], title="") :
    """
    Plots with linear regression, giving Pearson's r and its associated p-value,
    and prints them.
    
    INPUT:
    x, y - data
    xlabel - label for x
    ylabel - label for y
    yless - value to deduct to minimal y value in data, to place the text giving r and p-value.
    topredictX - y value, to print expected X
    ylim - y min and y max in plot span, if needed (otherwise matplotlib's default).
    title - title of the plot, if needed
    """
    p = linregress(x, y)
    plt.figure()
    plt.subplots_adjust(left=0.16)
    plt.scatter(x, y)
    plt.plot(x, x*p.slope+p.intercept)
    plt.text(np.mean(x) - np.mean(x)/3. , np.min(y)-yless, "r = "+str(round(p.rvalue,4)))
    plt.text(np.mean(x) + np.mean(x)/3. , np.min(y)-yless, "p-val = "+str(round(p.pvalue,10)))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if len(ylim) != 0 :
        plt.ylim(ylim[0],ylim[1])
    if len(title) != 0:
        plt.title(title)
    plt.show()
    if p.pvalue < 0.05:
        print str(round(p.rvalue,4))+"\t"+str(round(p.pvalue,4)), ylabel#+"\t"+str(p.slope)
    if topredictX :
        print "Expected", xlabel, ":", (topredictX - p.intercept)/float(p.slope)
    
        
    
def plotting_pathway_non_pathway_graph (x, y, y_path, xlabel, ylabel, 
                                        ylim=[], title=""):
    """
    Plots with linear regression, giving Pearson's r and its associated p-value.
    Plots at the same time data for graph with all enzymes (y) and for graph with
    only enzymes in known pathways.
    
    INPUT:
    x, y - data
    y_path - data for graphs for enzymes in known pathways
    xlabel - label for x
    ylabel - label for y
    ylim - y min and y max in plot span, if needed (otherwise matplotlib's default).
    title - title of the plot, if needed
    """
    p = linregress(x, y)
    p2 = linregress(x, y_path)
    plt.figure()
    plt.scatter(x, y, label = "All enzymes", c = "b" )
    plt.scatter(x, y_path, label = "Pathway enzymes", c = "g")
    plt.plot(x, x*p.slope+p.intercept, c = "b")
    plt.plot(x, x*p2.slope+p2.intercept, c = "g")
    plt.legend(loc = "best")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if len(ylim) != 0 :
        plt.ylim(ylim[0], ylim[1])
    if len(title) != 0:
        plt.title(title)
    plt.show()
    print "All", str(p.rvalue)+"\t"+str(p.pvalue)+"\t"+str(p.slope)
    print "pathway", str(p2.rvalue)+"\t"+str(p2.pvalue)+"\t"+str(p2.slope)
    
    
def common_nodes (nodes1, nodes2) :
    """
    Returns list of the intersection of two sets/lists of nodes
    """
    nodes1 = set(nodes1)
    nodes2 = set(nodes2)
    return nodes1 & nodes2

# =============================================================================
# Additional reconstruction functions
# =============================================================================



def get_organism_name(code):
    """
    From an organism code get organism name
    """
    k = KEGG()
    entry = k.get("genome:" + code)
    if type(entry) != type(1):
        i_name = entry.find("DEFINITION")
        if i_name != -1:
            org_name = entry[i_name+10:].split("\n")[0].lstrip(" ").replace("(", "--").replace(")", "").replace("/", "-")
        else:
            logger.error("No DEFINITION field (organism name) found")
            raise SystemExit()
    else:
        logger.error("No hits in KEGG genome database!")
        raise SystemExit()
    return org_name



def get_species(code_list) :
    """
    From KEGG organism code name, get organism name.
    
    INPUT:
    code_list - list of KEGG organism codes
    
    OUTPUT:
    species - list of corresponding KEGG species names
    duplicates - list of duplicate species (same species for multiple codes)
    dupli_index - list of duplicate species index in code_list
    """
    k = KEGG()
    org_list = k.list("organism")
    org = pd.read_csv(StringIO(org_list), sep="\t")
    org.columns = ["Species_code", "Letter_code", "Species", "Taxonomy"]
    species = []
    duplicates = []
    dupli_index = []
    for code in code_list :
        spc = org.loc[org["Letter_code"] == code, "Species"].to_string().lstrip("0123456789 ")
        i = spc.find("(")
        if i != -1 :
            spc = spc[:i]
        if not spc.startswith("Series"):
            if spc.rstrip(". ") in species : #duplicate
                if spc.rstrip(". ") in duplicates : #existing duplicate (seen more than twice)
                    i = duplicates.index(spc.rstrip(". "))
                    dupli_index[i].append(len(species))
                else :
                    j = species.index(spc.rstrip(". "))
                    dupli_index.append([j, len(species)])
                    duplicates.append(spc.rstrip(". "))
            species.append(spc.rstrip(". "))
        else:
            species.append("None")
    return species, duplicates, dupli_index



# Modified function from MetabolicGraph class : does not test validity of code through
# testing code with genes, but tests validity by checking if the taxonomy id in KEGG
# entry is the same as the one from a list of taxonomy ids found in BacDive.
def get_organism(org_name, tax_id_list, val = True):
    """
    Finds the KEGG organism code name through the organism name. Tests hits found
    through validating that the tax id (potential species) for the KEGG species 
    is in tax_id_list (actual species).
    
    INPUT:
    org_name - name of organism or parts of it, space separated
    tax_id_list - list of tax ids (we get them from BacDive db)
    val - boolean. If True, evaluates and validates code found. Else takes taxonomy
            ID found. Defaults to True.
    
    OUTPUT:
    code - KEGG organism ID/code or None, if not found
    """
    k = KEGG()
    org_list = k.lookfor_organism(org_name)
    if len(org_list) == 1: #Found possible organism hits
        code = org_list[0].split()[1]
        entry = k.get("genome:" + code)
        if type(entry) != type(1):
            i_tax = entry.find("TAX:") 
            if i_tax != -1:
                tax_id = entry[i_tax+4:].split("\n")[0]
                if not val or tax_id in tax_id_list:
                    return code, tax_id
    elif len(org_list) > 1 :
        for org in org_list :
            code = org.split()[1]
            entry = k.get("genome:" + code)
            if type(entry) != type(1):
                i_tax = entry.find("TAX:") 
                if i_tax != -1:
                    tax_id = entry[i_tax+4:].split("\n")[0]
                    if not val or tax_id in tax_id_list :
                        return code, tax_id
    return None, None


def find_get_organism(org_name, tax_id_list, syn, val = True):
    """
    Tries to find a KEGG species code for a species through the species name.
    If it doesn't find it with whole name, tries parts of name, if still not found,
    tests synonym names, otherwise raises an error. As an option, can choose not to
    evaluate validity of tax ID (when you don't have to compare to some other
    data linked to this tax ID).
    
    INPUT:
    org_name - organism name
    tax_id_list - list of organism tax ids, for validation of KEGG code
    syn - list of synonym names for species
    
    OUTPUT:
    code - KEGG code
    tax_id - NCBI tax ID
    """
#    logger.info("Looking for organism code in KEGG...")
    code, tax_id = get_organism(org_name, tax_id_list, val)
    if code == None:
        org_name_list = org_name.split()
        if val :
            org_name_list.append(org_name_list.pop(0)) #reshuffling : putting the first element of the name (usually the genus) as the last one to test, as it probably has a lot more hits
        else :
            org_name_list.pop(0)
#        logger.info("No hits for whole organism name, testing with parts of name...")
        for name in org_name_list: #Test parts of name
            code, tax_id = get_organism(name, tax_id_list, val)
            if code != None:
                return code, tax_id
        if code == None:
            if len(syn) > 0 and syn[0] != org_name:
                new_name = syn[0]
                code, tax_id = find_get_organism(new_name, tax_id_list, syn, val)
                return code, tax_id
            else :
                logger.error("Uh oh! Organism name not found in KEGG database!")
                raise SystemExit()
                
    if code != None:
        logger.info("Organism code found!")
        return code, tax_id


def find_fasta_with_taxID(code, org_name, work_dir="./") :
    """
    From a KEGG organism code, download cDNA fasta file via its TAX ID, and 
    return a MetabolicGraph object.
    
    INPUT:
    code - KEGG organism code
    org_name - name of organism (in this case doesn't have to be accurate, as it is
               not use to parse or download)
    work_dir - directory where species folders will be stored. Defaults to current directory.
    OUTPUT:
    g - MetabolicGraph object
    """
    k = KEGG()
    entry = k.get("genome:" + code)
    if type(entry) != type(1):
        i_tax = entry.find("TAX:")
        if i_tax != -1:
            tax_id = entry[i_tax+4:].split("\n")[0]
            os.system("bash get_cDNA2.0.sh assembly_summary_prokaryotes.txt " 
                      + org_name.replace(" ", "_") + " " + tax_id)
            g = MetabolicGraph(org_name, "genomes_cdna/" 
                               + org_name.replace(" ", "_") 
                               + "_cds_from_genomic.fna", 
                               code=code, work_dir=work_dir)
            return g
        else:
            logger.error("No TAXONOMY ID found")
            raise SystemExit()
    else:
        logger.error("No hits in KEGG genome database!")
        raise SystemExit()
        
        
def find_fasta_with_assemblyID(code, org_name,
                               work_dir="./") :
    """
    From a KEGG organism code, download cDNA fasta file via its assembly ID, and 
    return a MetabolicGraph object.
    
    INPUT:
    code - KEGG organism code
    org_name - name of organism (in this case doesn't have to be accurate, as it is
               not use to parse or download)
    work_dir - directory where species folders will be stored. Defaults to current directory.
    OUTPUT:
    g - MetabolicGraph object
    """
    k = KEGG()
    entry = k.get("genome:" + code)
    if type(entry) != type(1):
        i_id = entry.find("Assembly:")
        if i_id != -1:
            assm_id = entry[i_id+9:].split("\n")[0].rstrip(")")
            os.system("bash get_cDNA3.0.sh assembly_summary_prokaryotes.txt " 
                      + org_name.replace(" ", "_") + " " + assm_id)
            g = MetabolicGraph(org_name, "genomes_cdna/" 
                               + org_name.replace(" ", "_") 
                               + "_cds_from_genomic.fna", 
                               code=code, work_dir=work_dir)
            return g
        else:
            logger.error("No assembly ID found")
            raise SystemExit()
    else:
        logger.error("No hits in KEGG genome database!")
        raise SystemExit()
        

def find_fasta_with_species_name(code, org_name, work_dir="./") :
    """
    From a KEGG organism code, download cDNA fasta file via its species name, and 
    return a MetabolicGraph object. Tries to find name in Ensembl assembly file.
    
    INPUT:
    code - KEGG organism code
    org_name - name of organism, needs to be accurate
    work_dir - directory where species folders will be stored. Defaults to current directory.
    OUTPUT:
    g - MetabolicGraph object
    """
    os.system("bash get_cDNA1.0.sh species.txt " + org_name.replace(" ", "_"))
    g = MetabolicGraph(org_name, 
                       "genomes_cdna/" + org_name.replace(" ", "_") + "_cds_from_genomic.fna", 
                       code=code, work_dir=work_dir)
    return g


def get_fasta(code, spec_name, work_dir="./", alternative_name=""):
    """
    Try to find cDNA fasta file among our 3 searching methods.
    
    Alternative name for the method using name to find fasta. Alternative name is then the only one tested.
    """
    try : #Find through assembly ID
        obj = find_fasta_with_assemblyID(code, spec_name, work_dir=work_dir)
    except (SystemExit, IOError, TypeError):
        try: #Find through taxonomy ID
            logger.info("Trying to find species with taxonomy ID...")
            obj = find_fasta_with_taxID(code, spec_name, work_dir=work_dir)
        except (SystemExit, IOError, TypeError):
            try: #Fing through name
                logger.info("Trying to find species with species name in Ensembl summary file")
                if alternative_name != "":
                    obj = find_fasta_with_species_name(code, alternative_name, work_dir=work_dir)
                else:
                    obj = find_fasta_with_species_name(code, spec_name, work_dir=work_dir)
            except (SystemExit, IOError, TypeError): #Not found fasta ok
                os.system("rm -f genomes_cdna" + spec_name.replace(" ", "_") + "_cds_from_genomic.fna")
                os.system("rm -rf " + work_dir + spec_name.replace(" ", "_")) #remove directory
                logger.error("Species %s %s will not be handled" %(spec_name, code))
                raise SystemExit()
    return obj


        
def _transform_temp_ph(values):
    """
    Extracts a single value from possible list of values from BacDive's pH and
    Temperature field (BacDive may give a single value, a range of values, or
    combinations).
    
    INPUT:
    values - list of optimum temperature/pH for a single species
    
    OUTPUT:
    values - mean of interval or single value of list of temperatures/pH
    """
    try :
        values = np.mean(np.array(values).astype(float))
    except ValueError: #Temperature/pH interval
        if len(values) == 1 :
            vals = np.array(values[0].split("-")).astype(float)
            values = np.mean(vals)
        elif len(values) == 2 : #a value and an interval - we take the value
            for val in values :
                if val.find("-") == -1:
                    values = float(val)
        else :
            logger.warning("Wrong value : %s"  %values)
            raise SystemExit()
    return values


def get_temperature_pH_bacdive(d, valid_spec):
    """
    Extracts a single value from possible list of values from BacDive's pH and
    Temperature field (BacDive may give a single value, a range of values, or
    combinations), for all species.
    
    INPUT:
    d - bacdive dictionnary (created with get_bacdive_db.py)
    valid_spec - list of species name, KEGG code and BacDive ID (key to the dictionary),
                 for a bunch of species
    
    OUTPUT:
    temp - list of single float temperatures per species
    ph - list of single float pH per species
    """
    temp = []
    ph = []
    for spec in valid_spec :
        temp.append(_transform_temp_ph(d[spec[2]]["Temperature"]))
        ph.append(_transform_temp_ph(d[spec[2]]["pH"]))
    return temp, ph

# =============================================================================
# Additional graph analysis functions
# =============================================================================

def article_density(g):
    """
    Edge density from Takemoto's 2007 and Tamura et al's 2009 article :
    D = E/N, (half of the mean degree since the sum of all degrees count twice
    every edge).
    
    Additional function for non GraphClass objects (e.g. here : undirected graph)
    """
    dens = g.number_of_edges()/float(g.number_of_nodes())
    return dens
   
    
def triadic_census(graph) :
    """
    Gives a list of each triad fraction for our graph, from the 16 possible 
    types of triads (see them  here: 
    http://www.stats.ox.ac.uk/~snijders/Trans_Triads_ha.pdf)
    """
    census = nx.triadic_census(graph)
    fractions = []
    triads = ['201',  '021C', '021D', '210', '120U', '030C', 
              '003', '300', '012', '021U', '120D', '102', 
              '111U', '030T', '120C', '111D']
    all_census = np.sum(census.values())
    for tr in triads :
        fractions.append(float(census[tr])/all_census)
    return fractions


def degree_exponent_ML(graph):
    """
    Takemoto's 2007 maximum likelihood estimate of degree exponent. Nodes with 
    degree 0 are ignored.
    
    INPUT :
    graph - networkx graph
    
    OUTPUT :
    gamma - degree exponent
    """
    deg = np.array(dict(graph.degree).values())
    deg = deg[np.where(deg != 0)[0]]
    kmin = np.min(deg)
    fraction_sum = np.sum(np.log(deg/kmin))
    gamma = 1 + len(deg) * 1./fraction_sum
    return gamma
    
    
def laplacian_spectrum(adj_mat, normalised, use_out_degree=False):
    """
    Calculates Laplacian matrix and deduce the module of the spectrum.
    
    INPUT:
    adj_mat - graph adjacency matrix (sparse or not)
    normalised - Boolean. Normalised Laplacian or not.
    use_out_degree - Boolean. Use out-degree to calculate Laplacian or not. Defaults to False.
    
    OUTPUT:
    lapl - Laplacian matrix (sparse if adj_mat spares, dense otherwise)
    spectr_mod - array of Laplacian eigenvalues (absolute value!)
    """
    #nx.directed_laplacian_matrix - but cannot use our adjacency matrices
    lapl = csgraph.laplacian(adj_mat, 
                             normed=normalised,
                             use_out_degree=use_out_degree)
    spectr_mod = np.linalg.eigvals(lapl.toarray()) #Laplacian eigenvalues
    return lapl, np.abs(spectr_mod)



def laplacian_spectrum_vector(adj_mat, normalised, use_out_degree):
    """
    Calculates Laplacian matrix and deduce eigenvectors and the module of the eigenvalues.
    
    
    INPUT:
    adj_mat - graph adjacency matrix (sparse)
    normalised - Boolean. Normalised Laplacian or not.
    use_out_degree - Boolean. Use out-degree to calculate Laplacian or not.
    
    OUTPUT:
    lapl - Laplacian matrix (sparse if adj_mat spares, dense otherwise)
    spectr_mod - array of Laplacian eigenvalues
    eig_vectors - array of eigenvectors associated to eigenvalues
    """
    #nx.directed_laplacian_matrix - but cannot use our adjacency matrices
    lapl = csgraph.laplacian(adj_mat, 
                             normed=normalised,
                             use_out_degree=use_out_degree)
    spectr_mod, eig_vects = np.linalg.eig(lapl.toarray()) #Laplacian eigenvalues and eigenvectors
    return lapl, np.abs(spectr_mod), eig_vects


def get_nodes_eigenvalue(eivals, eivects, all_nodes):
    """
    Look for nodes associated to null-sum eigenvectors, associated to eigenvalues = 1.
    
    INPUT:
    eivals - Laplacian eigenvalues, numpy vector
    eivects - vertical eigenvectors as a numpy array
    all_nodes - list of node names (union of all graph node names)
    
    OUTPUT:
    eigen_nodes - list of nodes associated to null-sum eigenvectors (for 
                  eigenvalues = 1) for a species
    """
    eigen_nodes = []    
    #Look for eigenvalues of value 1.
    ei_idx = np.where(np.round(eivals,6) == 1)[0]
    
    #Look for eigenvectors associated, of null sum (column eigenvectors)
    v_idx = np.where(np.sum(eivects[:, ei_idx], axis=0) == 0)[0]
    for v in v_idx:
        nodes = np.where(np.round(eivects[:, ei_idx[v]], 6) != 0)[0] #Look for non-null nodes associated to vector
        for n in nodes:
            eigen_nodes.append(all_nodes[n])
    return eigen_nodes





if __name__ == '__main__':
    pass