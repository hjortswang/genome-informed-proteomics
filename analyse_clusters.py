import pandas as pd
import re
import ast
import os
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from pymsaviz import MsaViz

from Bio.Align.Applications import MuscleCommandline


def prep_cluster(clusters):
      
    clusters.loc[:, 'genomes'] = clusters['genomes'].astype(str).apply(ast.literal_eval)
    clusters.loc[:, 'seguid weight'] = clusters['seguid weight'].astype(str).apply(ast.literal_eval)
    clusters.loc[:, 'cluster sequences'] = clusters['cluster sequences'].astype(str).apply(ast.literal_eval)
    clusters.loc[:, 'cluster seguids'] = clusters['cluster seguids'].astype(str).apply(ast.literal_eval)
    return clusters

