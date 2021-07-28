from nd2reader import ND2Reader
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from PIL import Image

# location of the nd2 files to convert
#parent_dir = 'data-folder/from-sharepoint/OneDrive_1_01-06-2021'
parent_dir = '/Users/lmcneill/Documents/svi/imaging/data-folder/from-sharepoint/OneDrive_1_01-06-2021'
dpi_choice = 50
filenames = []
filenames_base = []

for file in glob.glob(os.path.join(parent_dir, '*nd2')):
    filenames.append(file)

filenames = sorted(filenames)

for name in filenames:
    base=os.path.basename(name)
    os.path.splitext(base)
    name_base = os.path.splitext(base)[0]
    new_name = parent_dir+'/'+str(name_base)
    print(new_name)

    try:

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
            plt.savefig(str(new_name)+'-DAPI.tif',dpi = dpi_choice)
            plt.figure(figsize = (10,10))
            plt.gca().set_axis_off()
            plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
                hspace = 0, wspace = 0)
            plt.margins(0,0)
            plt.gca().xaxis.set_major_locator(plt.NullLocator())
            plt.gca().yaxis.set_major_locator(plt.NullLocator())
            plt.imshow(images[1],cmap = 'gray')
            plt.savefig(str(new_name)+'-MLH3.tif',dpi = dpi_choice)
            plt.figure(figsize = (10,10))
            plt.gca().set_axis_off()
            plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
                    hspace = 0, wspace = 0)
            plt.margins(0,0)
            plt.gca().xaxis.set_major_locator(plt.NullLocator())
            plt.gca().yaxis.set_major_locator(plt.NullLocator())
            plt.imshow(images[2],cmap = 'gray')
            plt.savefig(str(new_name)+'-SYCP3.tif', dpi = dpi_choice)

    except:
        print("Couldn't open file"+str(name))
