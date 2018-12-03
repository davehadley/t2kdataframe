## T2KDataFrame

Provides an custom datasource for [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html) to read T2K ND280 near detector analysis files.

## Compilation

A cmake configuration is provided to compile to a standalone library. Alternatively, all of the important code is in
`T2KDataSource.cxx` and `T2KDataSource.hxx` which can simply be pasted directly into an existing project.
The only dependencies are ROOT and oaAnalysis headers (which the CMake compilation will automatically generate).

To compile
```
cmake -DMAKE_PROJECT_FILE_NAME=<path_to_an_oaanalysis_file.root> /path/to/t2kdataframe
```
An oaAnalysis file is needed to generate the oaAnalysis headers using ROOTs [TFile::MakeProject](https://root.cern/doc/master/classROOT_1_1RDataFrame.html).

A preprocessor flag `T2K_USE_GENIE` can be set in `T2KDataSource.hxx` compile for GENIE input MC rather than NEUT.

## Usage

To read a set of files use `T2K::MakeT2KDataFrame`. This takes a list of input files, a boolean flag asking whether to
switch on/off truth information (when reading real data files this flag must be switched off) and optionally a list
of trees to load (commonly used trees are loaded by default).

See `example.cxx` and `example.py` for working examples.

## Predefined Columns

All branches in the loaded trees have columns with the name `<Tree Name>_<Branch Name>`. For example, global reconstructed objects are called `Global_PIDs`.

MakeT2KDataFrame predefines some additional columns for convenience. These include:

* `t2kreco` a list of all global reconstructed tracks in time with the current bunch.
* `t2kbunch` an integer with the current bunch number (typically 0 to 8).
* `t2ktruthtraj`, `t2ktruthvtx` and `t2kroovtx` contain lists of all truth vertices/trajectories in the spill. For data these will be empty containers.

Contact <d.r.hadley@warwick.ac.uk> if you have any questions.
