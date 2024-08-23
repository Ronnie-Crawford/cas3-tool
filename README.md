# Cas3 motif gene ranking and motif finding tool
For the Cas3 project in Parts Group, Wellcome Sanger Institute

## Environment
The exact conda environment is provided in the `environment.yml` file, which can be used to automatically
make a conda environment with `conda env create -f environment.yml`.
However, if you prefer to not use this method, the python3 packages needed are simply:
- numpy
- pandas
- pyfaidx

## Configuration
Please paste the paths to relevant files into the `config.json` file, a detailed breakdown of the option to come soon.

## Running
To run the tool, start in the main folder `cd cas3-tool`, and run `python src`.
