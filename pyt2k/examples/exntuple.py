from argparse import ArgumentParser
from glob import glob
import os

import ROOT

# Define some C++ functions that will be dynamically compiled and used later.
_CPPCODE = '''
    bool positioninfiducialvolume(const TLorentzVector &pos) {
    auto x = pos.X();
    auto y = pos.Y();
    auto z = pos.Z();
    return (abs(x) < 874.51
            && abs(y - 55.) < 874.51
            && 136.875 < z && z < 447.375);
}

    bool infiducialvolume(T2K::Reco *reco) {
        auto &pos = reco->FrontPosition;
    return positioninfiducialvolume(pos);
    }
'''

# parse command line
def parsecml():
    parser = ArgumentParser(description="Make an flat ntuple from an oaAnalysis/eventAnalysis file with a T2KDataFrame.")
    parser.add_argument("filelist", nargs="+", help="List of input files.", type=str)
    return parser.parse_args()

def getfilelist(flist):
    # return a sorted list of unique normalised absolute file names matching the input patterns
    return list(sorted(set(os.path.abspath(f) for p in flist for f in glob(p))))

def tostdvector(input, cls=str):
    # convert input iterable to type std::vector<cls> required by some ROOT methods
    result = ROOT.std.vector(cls)()
    for x in input:
        result.push_back(x)
    return result

def main():
    args = parsecml()
    # create data T2KDataFrame
    flist = tostdvector(getfilelist(args.filelist))
    df = ROOT.T2K.MakeT2KDataFrame(flist, # input file list
                                   False, # switch off truth
                                   tostdvector(["ReconDir/Global"]), # list of trees to include, only include Global
                                   )
    # Compile to C++ code defined above
    ROOT.gInterpreter.Declare(_CPPCODE)
    # Define analysis
    anal = (
        # Select all global tracks in the fiducial volume, sort then by momentum (highest first)
        df.Define("tracks_in_fv", """
            return ROOT::VecOps::Sort(Filter(t2kreco, infiducialvolume), 
            [](T2K::Reco *lhs, T2K::Reco *rhs) -> bool { return lhs->FrontMomentum > rhs->FrontMomentum; });""")
        # Only save events with at least 1 track in the fiducial volume
        .Filter("tracks_in_fv.size()>0")
        # Define variables to write out
        .Define("leading_track_momentum", "tracks_in_fv.at(0)->FrontMomentum")
        .Define("leading_track_position", "tracks_in_fv.at(0)->FrontPosition")
        )
    # Write selected events and variables to a ROOT file
    anal.Snapshot("exntuple", "exntuple.root",
                  # variables to write to disk
                  tostdvector(["leading_track_momentum", "leading_track_position"])
                  )
    return

if __name__ == "__main__":
    main()
