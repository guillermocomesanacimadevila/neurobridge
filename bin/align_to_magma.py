#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path
import argparse

# NORMAL MAGMA -> Positionally mapped genes per locus
# take in magma ref for each trait
# take in for each locus file (input)
# take in the gene list (1XG matrix)
# 1. Compute FDR qval for each gene in MAGMA refs
# 2. Convert FDR q-vals into p-vals with significance threshold (p < 0.05)
# 3. For genes in 1XG matrix compile cross-trait FDR adjusted p-vals
# 4. If cross-trait FDR corrected p-val concordance -> add gene








