import struct
import numpy as np
import os
import re
import sys
from os.path import join 
from anytree import *
from anytree.exporter import DictExporter
from anytree.importer import DictImporter
from collections import Counter

if len(sys.argv) != 2:
    sys.exit('USAGE: python improve_taxonomy.py kraken_db_dir')


# helper function to get a subset of data from the super string
def get_string_from_offset(super_string, location):
    read_string = ''
    read_loc = location
    while True:
        read_byte = super_string[read_loc:read_loc+1].decode()
        if read_byte == '\x00':
            break
        else:
            read_string += read_byte
            read_loc += 1 
    return(read_string)

# find the root of the list of nodes
def get_root(nodes, name_super_string, level_super_string):
    # construct tree structure from nodes
    new_name = get_string_from_offset(name_super_string, nodes[1,3])
    new_rank = get_string_from_offset(level_super_string, nodes[1,4])
    new_taxid = nodes[1, 5]
    root = Node(new_name, rank=new_rank, taxid=new_taxid)
    node_list = [root]
    for i in range(2, np.shape(nodes)[0]):
        new_name = get_string_from_offset(name_super_string, nodes[i,3])
        new_parent = node_list[nodes[i, 0]-1 ]
        new_rank = get_string_from_offset(level_super_string, nodes[i,4])
        new_taxid = nodes[i, 5]
        node_list.append(Node(new_name, parent=new_parent, rank=new_rank, taxid=new_taxid))
    return root

# helper function to fill taxonomy ranks
def fill_taxonomy(path, ranks, taxid, name, uhgg=False):
    if uhgg: 
        tax_dict = {'name': name, 
                    'taxid': taxid, 
                    'ro': 'root', 
                    'kingd': '',
                    'phyl': '',
                    'cla': '',
                    'ord': '',
                    'fami': '',
                    'gen': '',
                    'speci': '',
                    }
    else:
        tax_dict = {'name': name, 
                    'taxid': taxid, 
                    'root': 'root', 
                    'kingdom': '',
                    'phylum': '',
                    'class': '',
                    'order': '',
                    'family': '',
                    'genus': '',
                    'species': '',
                    'subspecies': ''
                    }
    orig_keys = list(tax_dict.keys())
    for a,b in zip(ranks, path):
        tax_dict[a]=b
    # check for domain
    if 'domain' in tax_dict.keys():
        tax_dict['kingdom'] = tax_dict['domain']
    # ensure no invalid ranks
    tax_dict = {key: value for key, value in tax_dict.items() if key in orig_keys}
    to_return = list(tax_dict.values())
    return(to_return)

def main():
    taxdir = sys.argv[1]

    data = open(join(taxdir, "taxo.k2d"), "rb").read()
    magic_validation = data[0:8].decode()
    assert magic_validation == 'K2TAXDAT',"Data structure not valid K2TAXDAT"

    # this assumes the data structure is the same
    # unpack the binary data
    read_nodes, read_names, read_levels = struct.unpack("6i", data[8:32])[::2]

    # node data begins at byte 32
    # read 14 4 byte integers for each node 
    node_offset = 32
    read_positions = 14
    read_bytes = read_positions * 4
    nodes = np.array([struct.unpack("14i", data[(node_offset+(i*read_bytes)):node_offset+(i*read_bytes)+read_bytes])[0::2] for i in range(read_nodes)])
    # position names start in the structure
    name_offset = node_offset + (read_nodes * read_bytes)

    # then read the rest of the data as text
    names = data[name_offset:name_offset + read_names -1].decode().split('\x00')
    levels = data[name_offset + read_names:len(data)-1].decode().split('\x00')
    # giant string of names and levels
    name_super_string = data[name_offset:name_offset + read_names]
    level_super_string = data[name_offset + read_names:]


    # save temporary data for later?
    savetemp = False
    if savetemp:
        np.savetxt(join(taxdir, 'nodes_array.txt'), np.array(nodes), delimiter='\t', fmt='%2d')
        with open(join(taxdir, 'names_from_binary.txt'),'w') as f:
            for i in names:
                print(i,file=f)
        with open(join(taxdir, 'levels_from_binary.txt'),'w') as f:
            for i in levels:
                print(i,file=f)

    #### MODIFICATIONS WE'RE MAKING TO THE TAXONOMY ####
    # 1) anything below a species gets set as rank subspecies
    # 2) Prune the taxonomy to a specified set of levels. 
    # 3) To fix cases of "unclassified XXX" or "environmental samples"
    #    being at weird places in the tree and not having a useful name, 
    #    use the heuristic:
    #        if the name is "unclassified X" or "environmental samples"
    #        and everything below it is a species, and parent is not a genus or species
    #        than mark that as a genus 
    # 4) a few custom changes
    #   4a) set crass-like viruses to genus 
    #   4b) set unclassified bacterial viruses to family
    #   4c) set all superkingdoms to kingdom 
    #   4d) set Fungi and metazoa to be direct descendants of root

    # detect if we're working with the UHGG database here, 
    # as the taxonomic levels are different. 
    if 'speci' in levels and 'gen' in levels:
        working_with_uhgg = True
        print ("Detected working with UHGG database")
        root_rank_name = 'ro'
        species_rank_name = 'speci'
        keep_ranks = ['cla', 'fami', 'gen', 'kingd', 'ord', 'phyl', 'ro', 'speci']
    else:
        working_with_uhgg = False
        root_rank_name = 'root'
        species_rank_name = 'species'
        keep_ranks = ['superkingdom', 'kingdom', 'domain', 'phylum', 'class', 'order', 'family',
                       'genus', 'subspecies', 'species']

    root = get_root(nodes, name_super_string, level_super_string)
    root.rank = root_rank_name
    all_levels = [node.rank for node in PreOrderIter(root)]
    

    # mark all levels below species as subspecies
    # find species, get all their children, mark them 
    species_nodes = [node for node in PreOrderIter(root, filter_=lambda n: n.rank==species_rank_name)]
    for n in species_nodes:
        for c in n.children:
            c.rank = 'subspecies'

    # Custom changes
    # set crass-like viruses to genus 
    try: 
        find(root, lambda node: node.name == "crAss-like viruses").rank='genus'
    except: 
        AttributeError 

    try:
        # set these viruses to family
        find(root, lambda node: node.name == "unclassified bacterial viruses").rank='family'
    except: 
        AttributeError 
    try:
        # set all superkingdoms to kingdom to make things simple
        for a in findall(root, lambda node: node.rank == "superkingdom"):
            a.rank = 'kingdom'
    except: 
        AttributeError 
    # set Fungi and metazoa to be direct descendants of root
    try:
        find(root, lambda node: node.name == "Metazoa").parent = root
        find(root, lambda node: node.name == "Fungi").parent = root
        find(root, lambda node: node.name == "Eukaryota").parent = None
    except: 
        AttributeError 

    # iteration to discover and fix nodes
    nl  = [n for n in LevelOrderIter(root)]
    current_node = root
    for n in nl[1:]:
        # print(n)
        if n.rank not in keep_ranks:
            # check for names we want to switch to genus
            # must be above the genus level currently
            # everything below must be species, subspecies or no rank
            # then we can set this node as a genus!
            if (re.match('environmental samples', n.name) or re.match('unclassified' , n.name)) and \
                n.parent.rank in ['superkingdom', 'kingdom', 'domain', 'phylum', 'class', 'order', 'family'] and \
                set([s.rank for s in PreOrderIter(n)]).issubset(set(['species', 'subspecies', 'no rank'])):
                        n.rank = 'genus'
                        # change name if environmental samples
                        if re.match('environmental samples', n.name):
                            n.name = n.parent.name + ': environmental samples'
                        # print(n)
            else:
                for c in n.children:
                    c.parent = n.parent
                n.parent = None
                # n.children = ()
        else:
            current_node = n

    # get data structure for each taxid
    all_paths = [[b.name for b in a.path] for a in PreOrderIter(root)]
    all_ranks = [[b.rank for b in a.path] for a in PreOrderIter(root)]
    taxids = [a.taxid for a in PreOrderIter(root)]
    names = [a.name for a in PreOrderIter(root)]

    # fill in a matrix of each tax level
    # 11 column matrix: 
    #  name, taxid, root, taxonomy(kingdom - species)
    all_filled = np.array([fill_taxonomy(a,b,c,d, working_with_uhgg) for a,b,c,d in zip(all_paths, all_ranks, taxids, names)])
    # and now its all in a beautiful and organized matrix!
    # should also have unclassified in it
    if working_with_uhgg: 
        to_add = [['unclassified', '0', '', '','','','','','','']]
    else:
        to_add = [['unclassified', '0', '', '','','','','','','','']]
    tax_array = np.append(to_add, all_filled, axis=0)
    np.savetxt(join(taxdir, 'taxonomy_array.tsv'), tax_array, delimiter='\t', fmt="%s")

    # print out a rendered version of the tree
    print(RenderTree(root, style=ContRoundStyle), file=open(join(taxdir, 'rendertree_pruned.txt'), 'w'))

    print('Done! :)')

if __name__ == '__main__':
    main()