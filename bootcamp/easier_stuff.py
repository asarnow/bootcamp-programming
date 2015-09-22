# coding: utf-8

# This file is for those of you who learned Python over the summer (you did that, right?).
# In this file, I've put all of the nitty-gritty details of what makes this website work.

# Except it doesn't work, because you need to write all the functions!

# Some of these functions will just make the w  ebsite easier to use. Some of them are
# important for the enrichment and clustering tasks that your teammates are working on.

# If you need any help, ask your team or a TA.


# (don't delete this but don't worry about it either)
import os  # a built-in module, for dealing with filenames
from collections import namedtuple
import numpy as np
from . import app  # this is part of the website guts


# These are all the files you have to work with. Go open them in a text editor so you can
# get a feel for what they look like, because you need to parse each one to turn on a
# piece of the website.

# A list of yeast genes, with standard names and short descriptions.
GENE_INFO = os.path.join(app.root_path, 'data', 'gene_info.txt')

# A file that maps from GOID to name, aspect (process/function/component), etc
GO_INFO = os.path.join(app.root_path, 'data', 'go_info.txt')

# A two-column file that maps GOID to yeast genes
GO_MEMBERSHIP = os.path.join(app.root_path, 'data', 'go_membership.txt')

# A many-columned file that contains experimental data (yeast microarrays). Each column
# (after the first) is a different experiment, and each row is a gene. The values are log2
# ratios versus untreated control.
EXPERIMENT_FILE = os.path.join(app.root_path, 'data', 'experiment_data.txt')

# Create global dictionary for data d and gene info g_infor
experiments = {}
Gene = namedtuple("Gene", "name description")
GoTerm = namedtuple("GoTerm", "term aspect definition")
genes = {}
data = {"genes": []}
go_terms = {}
go_map = {}


# return a list or dictionary that maps from the id of an experiment (an int: 0, 1, ..)
# to a list of (systematic name, fold-change value) tuples
# e.g. [[('YAL001C', -0.06), ('YAL002W', -0.3), ('YAL003W', -0.07), ... ],
#       [('YAL001C', -0.58), ('YAL002W', 0.23), ('YAL003W', -0.25), ... ],
#        ... ]
def experiment():
    # only run function if dictionary d is empty
    if len(experiments) == 0:
        # open EXPERIMENT_FILE, read only, Universal ending, as 'f'
        with open(EXPERIMENT_FILE, 'rU') as f:
            # parse f by line l
            for l in f:
                # check if l starts with 'Y'
                if l.startswith("Y"):
                    # split l into list by tabs and strip white space
                    tok = l.strip().split('\t')
                    data["genes"].append(tok[0])
                    cnt = 0
                    # parse list l
                    for val in tok[1:]:
                        # check if experiment ID Key exists in dictionary already
                        if cnt not in experiments:
                            # Create KEY if it does not exist in dictionary as empty list
                            experiments[cnt] = []
                        # add gene name and experiment value to expID key in dictionary
                        experiments[cnt].append((tok[0], float(val)))
                        # move onto next experiment ID, aka column
                        cnt += 1
        with open(EXPERIMENT_FILE, 'rU') as f:
            data["data"] = np.genfromtxt((" ".join(ln.split("\t")[1:]) for ln in f), skip_header=1)
    return experiments


def parse_genes():
    if len(genes) == 0:
        with open(GENE_INFO, 'rU') as f:
            for l in f:
                if l.startswith("Y"):
                    tok = l.strip().split('\t')
                    genes[tok[0]] = Gene(name=tok[1], description=tok[2])


# map from a gene's systematic name to its standard name
# e.g. gene_name('YGR188C') returns 'BUB1'
def gene_name(gene):
    parse_genes()
    return genes[gene]


# map from a gene's systematic name to a list of the values for that gene,
# across all of the experiments.
# e.g. gene_data('YGR188C') returns [-0.09, 0.2, -0.07, ... ]
def gene_data(gene):
    experiment()
    idx = data["genes"].index(gene)
    return list(data["data"][idx, :])


# map from a systematic name to some info about the gene (whatever you want),
# e.g  'YGR188C' -> 'Protein kinase involved in the cell cycle checkpoint into anaphase'
def gene_info(gene):
    parse_genes()
    return genes[gene].description


def parse_go():
    if len(go_map) == 0:
        with open(GO_MEMBERSHIP, 'rU') as f:
            for l in f:
                if l.startswith("Y"):
                    tok = l.strip().split('\t')
                    gene = tok[0]
                    goid = tok[1]
                    if gene not in go_map:
                        go_map[gene] = set()
                    if goid not in go_map:
                        go_map[goid] = set()
                    go_map[gene].add(goid)
                    go_map[goid].add(gene)
    if len(go_terms) == 0:
        with open(GO_INFO, 'rU') as f:
            for l in f:
                if l.startswith("GO"):
                    tok = l.strip().split('\t')
                    go_terms[tok[0]] = GoTerm(term=tok[1], aspect=tok[2], definition=tok[3])


# map from a systematic name to a list of GOIDs that the gene is associated with
# e.g. 'YGR188C' -> ['GO:0005694', 'GO:0000775', 'GO:0000778', ... ]
def gene_to_go(gene):
    parse_go()
    return list(go_map[gene])


# map from one of the GO aspects (P, F, and C, for Process, Function, Component),
# to a list of all the GOIDs in that aspect
# e.g. 'C' -> ['GO:0005737', 'GO:0005761', 'GO:0005763', ... ]
def go_aspect(aspect):
    parse_go()
    out = []
    for goid in go_terms:
        if go_terms[goid].aspect == aspect:
            out.append(goid)
    return out


# map from a GOID (e.g. GO:0005737) to a *tuple* of the term, aspect, and term definition
# e.g. 'GO:0005737' -> ('cytoplasm', 'C', 'All of the contents of a cell... (etc)'
def go_info(goid):
    parse_go()
    if goid in go_terms:
        return go_terms[goid]
    else:
        return None


# the reverse of the gene_to_go function: map from a GOID
# to a list of genes (systematic names)
# e.g. 'GO:0005737' -> ['YAL001C', 'YAL002W', 'YAL003W', ... ]
def go_to_gene(goid):
    parse_go()
    return list(go_map[goid])


def get_go_map():
    parse_go()
    return go_map

