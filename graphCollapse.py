#!/usr/bin/env python

import pandas as pd
import argparse
import sys

####
#
# Parses arguments from command line
#
####

parser = argparse.ArgumentParser(description = 'graphCollapse: drops specified \
    edges from a gfa, the connects newly uniquley connected edges and writes \
    the final collapsed sequence to a fasta')

parser.add_argument('--gfa', help = 'input graph gfa file')
parser.add_argument('--drop_list', help = 'file listing edges to drop from \
    graph. File should have one edge name per line')
parser.add_argument('--min_length', default = 0, type = int,
    help = 'minimum length for output contigs (default = 0)')

args = parser.parse_args()

####
#
# Defines accessory funtions for sequence manipulation
#
####

def reverse_compliment(sequence):
    pairdict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    upperseq = sequence.upper()
    newseq = ""
    for i in range(1,len(sequence) +1):
        newseq += pairdict[upperseq[-i]]
    return newseq

def revor(orient):
    if orient == "+":
        return "-"
    elif orient == "-":
        return "+"

def tig_extend(tigname,link_df,seq_dict,directions,orient):
    tig_links = link_df[link_df['edge1'] == tigname]
    used_tigs = [tigname]
    if orient == '+':
        newtig = seq_dict[tigname]
    elif orient == '-':
        newtig = reverse_compliment(seq_dict[tigname])
    if 'right' in directions:
        extend_df = tig_links[tig_links['orient1'] == orient]
        if len(extend_df.index) == 1:
            next_orient = extend_df['orient2'].iloc[0]
            next_edge = extend_df['edge2'].iloc[0]
            recip_df = link_df[link_df['edge2'] == next_edge]
            recip_df = recip_df[recip_df['orient2'] == next_orient]
            if len(recip_df.index) == 1 and not next_edge == edge:
                extender = tig_extend(next_edge,link_df,seq_dict,'right',next_orient)
                newtig += "NNNNNNNNNN" + extender[0]
                used_tigs += extender[1]
    if 'left' in directions:
        extend_df = tig_links[tig_links['orient1'] == revor(orient)]
        if len(extend_df.index) == 1:
            next_orient = revor(extend_df['orient2'].iloc[0])
            next_edge = extend_df['edge2'].iloc[0]
            recip_df = link_df[link_df['edge2'] == next_edge]
            recip_df = recip_df[recip_df['orient2'] == revor(next_orient)]
            if next_edge == 'edge_91':
                debug_return = recip_df
            if len(recip_df.index) == 1 and not next_edge == edge:
                extender = tig_extend(next_edge,link_df,seq_dict,'left',next_orient)
                newtig = extender[0] + "NNNNNNNNNN" + newtig
                used_tigs = extender[1] + used_tigs
    return [newtig,used_tigs]

####
#
# Reads in graph lines and drops edges or connections to edges in the drop list
#
####

edge_list = []
for line in open(args.drop_list):
    edge_list.append(line.replace('\r','').replace('\n',''))
while "" in edge_list:
    edge_list.remove("")

seq_dict = {}
link_dict = {}
count = 0
for line in open(args.gfa):
    fields = line.split()
    if line[0] == 'S':
        if not line.split()[1] in edge_list:
            seq_dict[fields[1]] = fields[2]
    elif line[0] == 'L':
        if not line.split()[1] in edge_list and not line.split()[3] in edge_list:
            link_dict[count] = fields[1:]
            count += 1

link_df = pd.DataFrame.from_dict(link_dict,orient='index')
link_df.columns = ['edge1','orient1','edge2','orient2','overlap']

####
#
# Connects uniquely connected edges into new contigs and writes to stdout
#
####


newtigs = {}
used_edges = []

count = 0
for edge in seq_dict.keys():
    if not edge in used_edges:
        tigname = "contig_" + str(count)
        count += 1
        try:
            extendee = tig_extend(edge,link_df,seq_dict,'rightleft','+')
        except:
            sys.stderr.write(edge)
            break
        newtigs[tigname] = extendee[0]
        used_edges.extend(extendee[1])

for tig in newtigs:
    if len(newtigs[tig]) > args.min_length:
        sys.stdout.write(">" + tig + "\n" + newtigs[tig] + "\n")
