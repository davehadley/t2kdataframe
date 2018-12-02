## T2KDataFrame

Provides an custom datasource for RDataFrame to read T2K ND280 near detector analysis files.

## Compilation

A cmake configuration is provided to compile to library standalone. Alternatively, all of the important code is in
`T2KDataSource.cxx` and `T2KDataSource.hxx` which can simply be pasted directly into an existing project.
The only dependencies are ROOT and oaAnalysis headers (which the CMake compilation will automatically generate).

To compile
```
cmake -DMAKE_PROJECT_FILE_NAME=<path_to_an_oaanalysis_file.root> /path/to/t2kdataframe
```
An oaAnalysis file is needed to generate the oaAnalysis headers using ROOT TFile::MakeProject.

A preprocessor flag `T2K_USE_GENIE` can be set to compile for GENIE input MC rather than NEUT.

## Usage

To read a set of files use `T2K::MakeT2KDataFrame`. This takes a list of input files, a boolean flag asking whether to
switch on/off truth information (when reading real data files this flag must be switch off) and optionally a list
of trees to load (commonly used trees are loaded by default).

See `example.cxx` and `example.py` for working examples.