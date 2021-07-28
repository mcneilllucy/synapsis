# External scripts

This directory contains extra scripts that might be required for data preparation. At the moment, these are all python scripts.

## Getting started with python

Here is a short, straightforward way to install python to launch a jupyter notebook on a mac.
If anything is unclear or doesn't work, please let us know.

For terminal comments, just type everything *after* the $ sign exactly as it appears.

### Step 1: check whether you have python

First, open a terminal and type

```r
$ which python
```

if you get something like

```r
$ which python
/Users/lmcneill/opt/anaconda3/bin/python
```

then you already have python installed. Check for anaconda.


Otherwise, if you got a blank response e.g.:

```r
$ which python
$
```
then you'll need to install python.

### Step 2: installing python with anaconda3

We recommmend installing via anaconda, where the direct download .zip can be found [here](https://www.anaconda.com/products/individual#macos).

Once you have successfully gone through the installation steps in the graphical installer, *close the terminal, open a new one*, and then repeat step 1. You should get a path to where python is now i.e.

```r
$ which python
/Users/lmcneill/opt/anaconda3/bin/python
```

## `nd2converter.ipynb`

at the moment, anaconda installed most of the packages we will need for this. Except for nd2reader. So we return to the terminal and type:

```r
$ conda config --add channels conda-forge
```
After this, we can install packages from the command line. To install the python package `nd2reader`, type:

```r
$ conda install nd2reader
```

### Step 3: Launching jupyter notebook

Once you have installed `nd2converter.ipynb`, put it somewhere that you will remember. Then, in the terminal, navigate to the folder that the file is located. For example, if I put it in a folder called jupyter-notebooks in

```r
'/Users/lmcneill/Documents/svi/jupyter-notebooks'
```
then I would open a terminal and type

```r
$ cd Documents/svi/jupyter-notebooks
```
to check that the notebook is there,  type

```r
$ ls
```
and it should return all the files you have in there, including our notebook

```r
$ ls
nd2converter.ipynb
```
OK, we are ready to launch jupyter notebook. Type:

```r
$ jupyter notebook
```
A window should open in your browser. Click on `nd2converter.ipynb`.

When you have changed the path to the images you want to analyse i.e. change

```r
parent_dir = '/Users/lmcneill/Documents/svi/imaging/data-folder/from-sharepoint/OneDrive_1_01-06-2021'
```

in the python script, click inside the cell and press Shift+Enter. The notebook should run for a few mins. You can look in parent_dir in Finder to watch the .nd2s get converted.

## Navigating the terminal

By the way, if you're new to using the mac terminal, here are some commands that will help you.

To list all files and directories in your current location, type

```r
$ ls
```

to go into a folder (or complete path) in your current directory (i.e. appears in ls)
```r
$ cd Documents
```

to go to the previous folder
```r
$ cd ..
```

to go home
```r
$ cd
```

to make a folder
```r
$ mkdir
```

to make a file
```r
$ touch test.md
```

to remove a file (a bit dangerous)
```r
$ rm test.md
```

to print working directory (helpful for getting the full path to use `parent_dir` in `nd2converter.ipynb`)
```r
$ pwd
```
