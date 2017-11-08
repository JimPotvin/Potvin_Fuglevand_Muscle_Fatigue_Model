# Python Port

This is a direct port of the Matlab version of the model. Initially created with [smop](https://github.com/victorlei/smop) 
the code has been cleaned and modified to be readable and produce output values as close as possible to the original. (+/- 1e-13)

## Getting Started

These instructions will get you a copy of the port up and running on your local machine for development and testing purposes.

### Prerequisites

[Python 2.7](https://www.python.org/download/releases/2.7/)

We strongly recommend [pipenv](https://github.com/kennethreitz/pipenv) for Python package management.

```
pip install pipenv
```

### Installing

Use pipenv to install all the dependencies in a well isolated environment

```
cd python/
pipenv install
```

Enter that environment

```
pipenv shell
```

Run the script

NOTE: With the default values this will take over an hour and can take significantly longer than the Matlab version
```
python Python_Code_Potvin_Fuglevand_2017_Muscle_Fatigue_Model.py
```

To exit the environment

```
exit
```

## Running the tests

An output comparison tool is included to ensure the output of Matlab and Python runs are equivalent.

Artifacts from previous validation runs are included so you can immediately run the following

```
cd tests/
python compare_artifacts.py
```

If you modify input parameters in the Matlab and Python versions you can compare the output of that run as well.

To add a new run for comparison create new directories in `artifacts/`

```
mkdir -p tests/artifacts/<RUN_NAME>/matlab tests/artifacts/<RUN_NAME>/python
```

Copy the csvs output by the modified Matlab and Python runs and move them into their respective folders.

```
# Assuming you're in the python/ directory, copy Matlab output files
cp ../*.csv tests/artifacts/<RUN_NAME>/matlab/
# Copy Python output files
cp python*.csv tests/artifacts/<RUN_NAME>/python/
# Re-run comparison
python compare_artifacts.py
```

Files within `matlab/` and `python/` will be directly compared.

For each csv within matlab/ there must be an equivalent csv in python/ which
has the same name plus the prefix 'python - '.

## Authors

* **Ian Danforth** - github - [iandanforth](https://github.com/iandanforth)
