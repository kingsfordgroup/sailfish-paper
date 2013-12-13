This is a repository for the scripts, configuration files and data necessary to reproduce the
results reported in the Sailfish manuscript.  We are currently using [SnakeMake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) for
for the data workflows, while most of the specific file processing and plotting scripts are written in Python.  Unfortunately, the 
plotting and processing scripts are written in Python 2.7, while SnakeMake requires Python 3 --- meaning that both will be required to
run the full pipelines.  This repository will be continually updated as we automate all of the tests.
