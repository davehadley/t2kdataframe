import ROOT

from argparse import ArgumentParser
from glob import glob
import os

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
    "TruthDir/Trajectories",
    "TruthDir/Vertices"]

def main():
    args = parsecml()
    flist = tostdvector(getfilelist(args.filelist))
    df = ROOT.T2K.MakeT2KDataFrame(flist,
                                   True,
                                   tostdvector(listoftrees()),
                                   ROOT.T2K.BunchTiming(),
                                   )
    ROOT.gInterpreter.Declare(_CPPCODE)
    anal = (df.Define("tracks_in_fv", "return ROOT::VecOps::Sort(Filter(t2kreco, infiducialvolume), [](T2K::Reco *lhs, T2K::Reco *rhs) -> bool { return lhs->FrontMomentum > rhs->FrontMomentum; });")
            .Filter("tracks_in_fv.size()>0")
            .Define("leading_track_momentum", "tracks_in_fv.at(0)->FrontMomentum")
            )
    anal.Snapshot("exntuple", "exntuple.root", tostdvector(["leading_track_momentum"]))
    return

if __name__ == "__main__":
    main()
