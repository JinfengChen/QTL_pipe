#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import re
import os
import argparse
from Bio import SeqIO
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.spatial.distance as distance
import matplotlib as mp
mp.use('pdf')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pandas as pd


def usage():
    test="name"
    message='''
perl cluster.py --input Problem_lulu.list.tranpose.dist

    '''
    print message

def set_ticks_XY(ax, ypos, ylim, chrs):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(True)
    ax.yaxis.set_ticks_position('right')
    ax.xaxis.set_ticks(chrs[2], minor=True)

    # put the major ticks at the middle of each cell
    ax.set_yticks(ypos, minor=False)
    ax.set_xticks(chrs[1], minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels(chrs[0], minor=False)
    ax.set_yticklabels(ylim, minor=False)

    # rotate the
    #plt.xticks(rotation=0)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()


    for t in ax.xaxis.get_major_ticks():
    #    t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
 

    return ax


def set_ticks_XY_Right(ax, ypos, ylim):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(True)
    ax.yaxis.set_ticks_position('right')

    # put the major ticks at the middle of each cell
    ax.set_yticks(ypos, minor=False)
    #ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels([], minor=False)
    ax.set_yticklabels(ylim, minor=False)

    # rotate the
    #plt.xticks(rotation=0)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax

def set_ticks_XY_empty(ax):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    #ax.set_yticks(np.arange(xlen) + 0.5, minor=False)
    #ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels([], minor=False)
    #ax.set_yticklabels([], minor=False)

    # Set y ticklabels font size
    for tck in ax.get_yticklabels():
        tck.set_fontsize(2)

    # rotate the
    plt.xticks(rotation=90)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax



def set_ticks_XY_empty0(ax):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    #ax.set_yticks(np.arange(xlen) + 0.5, minor=False)
    #ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels([], minor=False)
    ax.set_yticklabels([], minor=False)

    # rotate the
    plt.xticks(rotation=90)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax



def tree(mat, lab):
    # clustering
    dist_mat = mat
    linkage_matrix = linkage(dist_mat, "single")

    # Create a figure.
    figsize=(20,20)
    fig = plt.figure(figsize=figsize)
    plt.subplots_adjust(bottom=0.2)
  
    # split plot into 4 rows X 2 cols. gs[1, 0] is 2th row and 1th column, which is left bottom region.
    # gs[0:, 0]  is all rows of 1th column
    ncurve = 274.00
    gs = gridspec.GridSpec(int(ncurve), 2, wspace=0.01, hspace=0.15, width_ratios=[0.15,1], height_ratios=[1/ncurve]*int(ncurve))

    #dendrogram
    ax0 = fig.add_subplot(gs[0:, 0])
    ddata = dendrogram(linkage_matrix,
                   color_threshold=1,
                   orientation='right',
                   labels=lab)
    ax0 = set_ticks_XY_empty(ax0)

    fig.savefig('%s.pdf' %('tree'), bbox_inches='tight')


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    #lodfile = '../input/MPR.cross.uniq.QTL.mr.table.new'
    #lodcutfile = '../input/MPR.cross.uniq.QTL.mr.table.LOD_threshold.new'
    #chrfile = '../input/MSU7.Chr.midpoint'
    pdist = pd.read_table(args.input, index_col=0, header=None)
    #midpoint = pd.read_table(chrfile, header=None)
    #binmap   = pd.read_table(args.bins)
    #lodcut= pd.read_table(lodcutfile)   
 
    mat = np.array([[100,  80,  50],
                [1,  0.8, 0.5],
                [30, 80,  80],
                [0.5,  0.8, 1]])

    #pairwise_dists = distance.squareform(distance.pdist(mat))
    tree(pdist, pdist.index)

if __name__ == '__main__':
    main()

