# HypoXPy
A modularized comprehensive *Python* interface for hypocenter location codes

## Disclaimer
*`hypoxpy`* is a research-oriented software package for improving earthquake locations and performing relative relocation using external location engines (e.g., HYPOINVERSE and HYPODD). It consumes phase picks and event catalogs produced by upstream tools (such as GAMMA or EQTransformer) but does not perform phase picking or event association itself. Results depend strongly on the quality of input picks, catalogs, velocity models, and parameter choices, and should be independently evaluated before scientific or operational use. This package doesn't intend to be a detailed documentation of each earthquake location method and code. The ultimate goal is to make this package a one-stop shop for earthquake location/relocation with a comprehensive collection of methods. The logic of this package is heavily inspired by, and modified and simplified from Hypo-Interface-Py (https://github.com/YijianZhou/Hypo-Interface-Py). We acknowledge the developers of Hypo-Interface-Py.

## File structure
* example: contains example scripts and jupyter notebooks running the hypoxpy workflow to relocate earthquakes.
* figs: contains figures used in this README file.
* hypoxpy: contains modules and templates of parameter files running hypoinverse and hypodd
* src: contains source codes for the earthquake lcoation/relocation methods.

## Installation
### Methods supported: 
* HypoInverse: HypoInverse earthquake location/relocation package. You can download the version from the USGS website (https://www.usgs.gov/software/hypoinverse-earthquake-location). For your convenience, a copy of the hyp1.40 codes is available under the folder `src/hypoinverse1.40`. This method is labeled as `hypoinverse` in our package.

* HypoDD: earthquake relocation using the double-difference method. developed by Felix Waldhauser. See information on the references at the end of the Readme. A copy of the version (2.1) is included in `src/hypodd2.1`, after debugging some variable types from the version on https://github.com/fwaldhauser/HypoDD. This method is labeled as `hypodd` in our package.

To use each included method, please make sure the corresponding packages are installed and working. Follow the related instructions to install each individual package.
 

### Install as a pip package
HypoXPy is available on pypl as a standalone package. It could be installed as regular pip package: `pip install hypoxpy`. The latest version (pre-release) is always available on GitHub.

### Install with local copy from GitHub
This installation method will get the latest version from GitHub.

1. Create and activate a Python virtual environment `hypox`. This is optional but recommended. This will help isolate the computational needs from other packages. This will also help avoid version incompatibility in case new packages for some codes are updated in the `base` conda environment. However, since `HypoXPy` doesn't currently need complicated packages, running under the `base` environment may just work fine. 
```
$ conda create -n hypox -c conda-forge python=3.11 jupyter numpy scipy pandas cartopy obspy mpi4py tqdm
```
Then, activate the environment with: `$ conda activate hypox`

2. Clone the github repository, in terminal under the desired directory:
```
$ git clone https://github.com/xtyangpsp/HypoXPy.git
```

3. `cd` to the repository folder in terminal and install the package with the parameters in `setup.py`.
```
$ pip install .
```

4. Create jupyter notebook kernel with the environment (after activating the environment).
```
$ pip install --user ipykernel
$ python -m ipykernel install --user --name=hypox
```

## Examples
Please run the jupyter notebooks in the `example` folder after successfully installing the package.

## References:
Klein, Fred W. 2002. User’s Guide to HYPOINVERSE-2000, a Fortran Program to Solve for Earthquake Locations and Magnitudes. https://doi.org/10.3133/ofr02171.

Waldhauser, Felix, and William L. Ellsworth, A double-difference earthquake location algorithm: Method and application to the northern Hayward fault, California, Bull. Seism. Soc. Am. 90, 1353-1368, 2000.

Waldhauser, Felix, hypoDD -- A program to compute double-difference hypocenter locations, U.S. Geological Survey Open-File Report 01-113, 2001.