from nd2reader import ND2Reader
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from PIL import Image

parent_dir = 'data-folder/sharepoint'
filenames = []

for file in glob.glob(os.path.join(parent_dir, '*nd2')):
    filenames.append(file)
filenames = sorted(filenames)
print(filenames)

for name in filenames:

    with ND2Reader(str(name)) as images:
        plt.figure(figsize = (10,10))
        plt.gca().set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
                hspace = 0, wspace = 0)
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.imshow(images[0],cmap = 'gray')
    ### keep metadata... e.g.
    ## savefig(fname, dpi=None, facecolor='w', edgecolor='w',
    ## orientation='portrait', papertype=None, format=None,
    ## transparent=False, bbox_inches=None, pad_inches=0.1,
    ## frameon=None, metadata=None)
        plt.savefig(str(name)+'cells.jpeg')
        plt.figure(figsize = (10,10))
        plt.gca().set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
            hspace = 0, wspace = 0)
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.imshow(images[1],cmap = 'gray')
        plt.savefig(str(name)+'foci.jpeg')
        plt.figure(figsize = (10,10))
        plt.gca().set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
                hspace = 0, wspace = 0)
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.imshow(images[2],cmap = 'gray')
        plt.savefig(str(name)+'dna.jpeg')
