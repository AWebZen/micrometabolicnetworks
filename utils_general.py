#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 09:45:33 2019

@author: adele
"""
import logging
import collections

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.stats import linregress
from scipy.sparse import issparse
from bioservices import KEGG


#Setting logging preferences
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

k = KEGG()
k.settings.TIMEOUT = 1000 #Changing timeout


COMPLETE_BACT_MEDIUM = [ #letort c et al, 2001 : https://mediadb.systemsbiology.net/defined_media/media/322/
"C00568", #4-Aminobenzoate
"C00147", #Adenine
"C00041", #Alanine
"C01342", "C00158", #Ammonium citrate
"C00062", #Arginine
"C00072", #Ascorbate
"C00152", #Asparagine
"C00049", #Aspartate
"C00120", #Biotin
"C08130", "C00076", "C00698", #Calcium chloride anhydrous, Ca, Cl
"C00864", #Calcium pantothenate
"C00175", #Cobalt chloride: Cb
"C00070", "C00059", #Cupric sulfate: Cu, SO42-
"C00097", #Cysteine
"C00740", "C00065", #DL-Serine
"C06420", "C00082", #DL-Tyrosine
"C00855", "C00073", #DL-methionine
"C14818", "C14819",  #Ferrous chloride
"C00504", #Folate
"C00064", #Glutamine
"C00037", #Glycine
"C00242", #Guanine
"C00135", #Histidine
"C00294", "C00081", "C00104",#Inosine, Inosine TP, IDP
"C00407", #Isoleucine
"C00025", #L-Glutamate
"C00243", #Lactose
"C00123", #Leucine
"C00725", #Lipoate
"C00047", #Lysine
"C00305", #Magnesium chloride
"C00034", "C19610", "C19611", #Mn, Mn2+, Mn3+ :Manganese sulfate
"C00253", #Nicotinate
"C00295", #Orotate
"C00079", #Phenylalanine
"C13197", #Potassium dibasic phosphate
"C00238", "C00009", #Potassium dihydrogen phosphate
"C00148", #Proline
"C00534", #Pyridoxamine HCl
"C00314", #Pyridoxine HCl
"C00255", #Riboflavin
"C01330", "C00033", #Sodium acetate
"C00378", #Thiamine HCl
"C00188", #Threonine
"C00214", #Thymidine
"C00078", #Tryptophan
"C00106", #Uracil
"C00183", #Valine
"C05776", #Vitamin B12
"C00385", #Xanthine
"C00038", #Zinc sulfate
]




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
    
    
def get_domain(codes):
    """From list of KEGG codes, get organism domain"""
    domain = []
    orgs = k.list("organism")
    for code in codes:
        hit = orgs.find("\t"+code+"\t")
        if hit == -1:
            logger.error("Did not find species %s" %code)
            domain.append(None)
        else:
            dom = orgs[hit:].split("\n")[0].split("\t")[3].split(";")[1]
            domain.append(dom)
    return domain
    

def from_y_to_dict(y):
    """From vector of values, build dict where keys are the values and dict values
    are indexes of vector with such values"""
    return {k:np.where(y == k)[0] for k in np.unique(y)}


def jaccard_index(y1, y2):
    """Evaluate Jaccard Index between y1 and y2"""
    sum_array = np.array([y1, y2]).sum(axis=0)
    if np.sum(sum_array) == 0:
        return 0
    jacc_ind = np.sum(sum_array == 2)/float(np.sum(sum_array != 0))
    return jacc_ind


def sim_matrix(X):
    """
    For a species x other thing (nodes) matrix, build a species x species
    jaccard index similarity matrix of species vectors.
    """
    sim = np.zeros((len(X), len(X)))
    for i, spec in enumerate(X):
        print i
        for j, spec2 in enumerate(X):
            if j > i:
                break
            if issparse(spec):
                spec = spec.todense()
            if issparse(spec2):
                spec2 = spec2.todense()
            sim[i,j] = jaccard_index(spec, spec2)
    return sim + sim.T - np.diag(sim.diagonal()) #Symmetrizing triangular matrix


   
    

def get_taxID(code):
    """From KEGG species code, get Taxonomy ID"""
    global k
    entry = k.get("genome:" + code)
    if type(entry) != type(1):
        i_tax = entry.find("TAX:")
        if i_tax != -1:
            tax_id = entry[i_tax+4:].split("\n")[0]
    else:
        logger.error("Uh oh, organism entry not found for %s" %code)
        raise SystemExit()
    return tax_id




def reorder_matrix (m, d) :
    """
    Reorder matrix : put species in same cluster together.
    
    INPUT:
    m - matrix
    d - cluster dictionary : {cluster : [list of species index in cluster]}
    
    OUTPUT :
    m in new order
    new_order - order of species indexes in matrix
    """
    new_order = []
    for i, med_class in enumerate(d.values()):
        new_order.append(med_class)
    return m[np.concatenate(new_order), :], new_order
            
    

def plot_scope_matrix_clusters(scope, 
                               d,
                               parameters_txt_xlim=[-6,-7],
                               xbar=0.92,
                               xlab_acc="Nodes in scope",
                               ylab_acc="Species",
                               figsiz=(20,13),
                               subset_sp=[]):
    """
    Plots scope heatmap, with a colored line on the left to define clusters

    Parameters
    ----------
    scope : matrix
        Scope matrix, species x nodes. 
    d : dict
        {cluster : [list of species index in cluster]} dict.
    parameters_txt_xlim : list size 2, optional
        First number sets xlim min, so as to allow space for text. 
        Second number sets x for text (cluster names). The default is [-6,-7].
    xbar : float, optional
        Left position in figure of colour bar. The default is 0.92.
    xlab_acc : str, optional
        X label. The default is "Nodes in scope".
    ylab_acc : str, optional
        Y label. The default is "Species".
    figsiz : list size 2, optional
        Figure size. The default is (20,13).
    subset_sp : list, optional
        List of species indexes to plot, if only want a subset. The default is [].

    Returns
    -------
    None.

    """
    if subset_sp:
        d2 = {}
        for k in d.keys():
            d2[k] = np.array([sp for sp in d[k] if sp in subset_sp])
    else:
        d2 = d
    sim, new_order = reorder_matrix(scope, d2)
    colors = ["k", "mediumseagreen",  "dodgerblue", "crimson", "lightgrey", 
              "yellow", "peru", "magenta", "b", "darkorchid", "brown", 
              "olive", "wheat", "purple", "cadetblue", "pink", "red", "grey", 
              "turquoise", "lime", "orange", "salmon", "cyan", "g", "hotpink",
              "tan", "lavender", "teal", "darkorange", "seagreen"]

    fig = plt.figure(figsize=figsiz)
    
    gs = gridspec.GridSpec(nrows=1, ncols=1)
    
    #First subplot
    ax0 = plt.subplot(gs[0])
    im = ax0.imshow(sim, vmin=0, vmax=1)
    
    print im.axes.get_position(), im.axes.get_position(original=True)
    
    ax0.set_xlabel(xlab_acc)
    for tick in ax0.xaxis.get_major_ticks():
        tick.label1.set_horizontalalignment('right')
    ax0.set_ylabel(ylab_acc)
    ax0.set_ylim(-0.5, sim.shape[0])
    plt.subplots_adjust(left=0.07, right=0.9, top=0.98, bottom=.07)
    
    length = 0
    for i, med_class in enumerate(new_order):
        ax0.plot([-1, -1], [length, length+len(med_class)-1], lw = 4, c=colors[i])
        ax0.text(parameters_txt_xlim[0], (2*length+len(med_class)-1)/2., "Cl."+str(i), 
                 fontsize = 10, fontweight="bold", color=colors[i]) 
        length += len(med_class)
    ax0.set_xlim(parameters_txt_xlim[1], sim.shape[1]) 
    
    position = ax0.axes.get_position()
    print position, position.y0
    if figsiz[0] < 5:
        k = .05
    elif figsiz[0] < 10:
        k = .02
    else:
        k = .01
    cbax = fig.add_axes([xbar, position.y0, k, position.y1-position.y0]) 
    fig.colorbar(im, cax=cbax)
    fig.subplots_adjust(wspace=.0)
    
    