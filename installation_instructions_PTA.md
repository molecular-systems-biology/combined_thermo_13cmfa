to run flux sampling with PTA : 
# Installing PTA on the Cluster

These instructions guide you through setting up a Python virtual environment and building the PTA package with C++ bindings, using cluster modules. 
pta installation via pip install creates pb with cmakeversion. We need to clone pta and install it from there, having modified the proper files (see Laurens' github). 
Then, we need to make sure that all modules we load use the same c++ compiler?
So we first load modules then create our environement, so that path like gurobi_home are exported already.  Not sure that the paths stay if we want to rebuild the package later.


## 0. Clone PTA

``` git clone --recurse-submodules <pta-repo-url> cd <pta-repo-directory>```

modify the files required 

## 1. Load necessary modules

```module purge ml Gurobi (12)  ml R/4.4.1-gfbf-2023b```


    Load modules before creating the virtual environment.

## 2. Create a new Python virtual environment

conda create -n pta_env
conda activate pta_env

or on cluster
load Python/3.11.XX/XXX (XX is GCC version, which should match other modules, see below)
python -m venv ~/venvs/pta
source ~/venvs/pta/bin/activate

## 3. Upgrade Python packaging tools

python -m pip install --upgrade pip setuptools wheel

dependencies :
- hopsy
- pyarrow (to save in feather)

for pta : 
-ninja
-gurobipy  (not on cluster)

## 5. Install PTA package
 (while in the pta_cloned folder)

5A On cluster (recap and logic):
On cluster we use the Gurobi module to avoid loading our named licence.
So we don't install gurobipy but we want the pta installation to know where our gurobi will be. creating the environement with the Gurobi module active will automatically export GUROBI_HOME=$GUROBI_HOME such that the installer will find it. Otherwise we need to load modules then activate env.

also need on cluster
```export RPY2_CFFI_MODE=API```


we also have to make sure that both modules use a same gcc library than the python we used when creating the enviornment (check with "which python" then check what is the GCC in module spider Python/3.version/GCC)
then to install package we run
python -m pip install . --no-build-isolation 


## 6. Verify the installation

python -c "import pta; print('PTA installed successfully')"


##7 troubleshooting .
sometimes we need to remove the build files to build again : 
rm -rf build/ dist/ *.egg-info

