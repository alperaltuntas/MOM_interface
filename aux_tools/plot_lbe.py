#!/usr/bin/env python3

import os
import xarray as xr
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


from lbe_utils import MOM_define_layout, mpp_compute_extent

descr = """
MOM6 land block elimination preprocessor module
NOTE: this preprocessing module is no longer needed in practice (see AUTO_MASKTABLE option),
      but kept in MOM_interface for prototyping and diagnostics purposes.
"""

def read_mask_table(mask_file_path):

    with open(mask_file_path, 'r') as f:
        num_masked_blocks = int(f.readline())
        layout = [int(s) for s in f.readline().split(',')]
        masked_blocks = []
        for _ in range(num_masked_blocks):
            masked_blocks.append((int(s) for s in f.readline().split(',')))
        return layout, masked_blocks
        

def plot_mask_table(topo_file_path, mask_file_path, output_dir):

    ds_topog = xr.open_dataset(topo_file_path)
    da_mask = ds_topog.mask
    ny, nx = da_mask.shape

    layout, masked_blocks = read_mask_table(mask_file_path)
    ndivs = layout[0] * layout[1]
    idiv, jdiv = MOM_define_layout(nx, ny, ndivs)
    assert layout[0] == idiv
    assert layout[1] == jdiv

    ibegin, iend = mpp_compute_extent(1,nx,idiv)
    jbegin, jend = mpp_compute_extent(1,ny,jdiv)

    fig, ax = plt.subplots(figsize=(10,8))

    im1 = plt.imshow(da_mask, interpolation='nearest', origin='lower')

    for j in jbegin:
        plt.axhline(y = j-1, color = 'r', linestyle = '-')
    for i in ibegin:
        plt.axvline(x = i-1, color = 'r', linestyle = '-')

    for i,j in masked_blocks:
        i0, j0 = ibegin[i-1], jbegin[j-1], 
        iw, jw = iend[i-1]-i0, jend[j-1] - j0
        ax.add_patch(Rectangle((i0, jbegin[j-1]), iw, jw))    

    fig.tight_layout()
    plt.show()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('-t', metavar='topofile', type=str, required=True,
                        help='MOM6 topography file path (TOPO_FILE)')
    parser.add_argument('-m', metavar='masktable_file', type=str, required=True,
                        help='MOM6 mask_table file path to be visualized')
    parser.add_argument('-o', metavar='rundir', type=str, required=False,
                        help='Output directory')
    args = parser.parse_args()

    rundir = args.o or os.getcwd()
    plot_mask_table(args.t, args.m, rundir)
