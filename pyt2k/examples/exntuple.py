import ROOT

from argparse import ArgumentParser
from glob import glob
import os

def parsecml():
    parser = ArgumentParser()
    parser.add_argument("filelist", nargs="+", help="List of input files.", type=str)
    return parser.parse_args()

def getfilelist(flist):
    # return a sorted list of unique normalised absolute file names matching the input patterns
    return list(sorted(set(os.path.abspath(f) for p in flist for f in glob(p))))

def tostdvector(input, cls=str):
    result = ROOT.std.vector(cls)()
    for x in input:
        result.push_back(x)
    return result

def listoftrees():
    return ["HeaderDir/BasicHeader",
    "HeaderDir/BasicDataQuality",
    "HeaderDir/BeamSummaryData",
    "ReconDir/Global",
    "TruthDir/GRooTrackerVtx",
    #"TruthDir/NRooTrackerVtx",
    "TruthDir/Trajectories",
    "TruthDir/Vertices"]


def main():
    args = parsecml()
    flist = tostdvector(getfilelist(args.filelist))
    df = ROOT.T2K.MakeT2KDataFrame(flist,
                                   True,
                                   tostdvector(listoftrees()),
                                   ROOT.T2K.BunchTiming(),
                                   True,
                                   )
    return

if __name__ == "__main__":
    main()
