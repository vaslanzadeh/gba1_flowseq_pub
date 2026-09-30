import os
import re
import time
import math
import warnings
import collections

import matplotlib.pyplot as plt
import numpy as np


import alignparse 
import alignparse.ccs
import alignparse.targets
import alignparse.minimap2
import alignparse.consensus
from alignparse.constants import CBPALETTE


import dms_variants
import dms_variants.utils
import dms_variants.plotnine_themes
import dms_variants.codonvarianttable

import pandas as pd
import numpy
import yaml


from IPython.display import display, HTML
from plotnine import *

from Bio.Seq import Seq


warnings.simplefilter('ignore')

theme_set(dms_variants.plotnine_themes.theme_graygrid())

print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")

with open('config.yaml') as f:
    config = yaml.safe_load(f)

os.makedirs(config['figs_dir'], exist_ok=True)
os.makedirs(config['process_ccs_dir'], exist_ok=True)

print(f"Reading amplicons from {config['amplicons']}")
print(f"Reading feature parse specs from {config['feature_parse_specs']}")

targets = alignparse.targets.Targets(
        seqsfile=config['amplicons'],
        feature_parse_specs=config['feature_parse_specs'])

fig = targets.plot(ax_width=7,
                   plots_indexing='biopython',  # numbering starts at 0
                   ax_height=2,  # height of each plot
                   hspace=1.2,  # vertical space between plots
                   )

plotfile = os.path.join(config['figs_dir'], 'amplicon.pdf')
print(f"Saving plot to {plotfile}")
fig.savefig(plotfile, bbox_inches='tight')

print(targets.feature_parse_specs('yaml'))

pacbio_runs = pd.read_csv(config['pacbio_runs'], dtype=str)
fastq=pacbio_runs['fastq']

print ("Reading amplicons from: " + (pacbio_runs['fastq'][0]))
print ("Reading ccs report from: " + pacbio_runs['report'][0])

mapper = alignparse.minimap2.Mapper(alignparse.minimap2.OPTIONS_CODON_DMS)

print(
    f"Using `minimap2` {mapper.version} with these options:\n"
    + " ".join(mapper.options)
)

readstats, aligned, filtered = targets.align_and_parse(
        df=pacbio_runs,
        mapper=mapper,
        outdir=config['process_ccs_dir'],
        name_col='run',
        group_cols=['name', 'library'],
        queryfile_col='fastq',
        overwrite=True,
        ncpus=config['max_cpus'],
        )
