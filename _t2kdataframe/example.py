from argparse import ArgumentParser
import ROOT

FIDUCIALVOLUMECODE = """
bool infiducialvolume(const TLorentzVector& pos) {
    auto x = pos.X();
    auto y = pos.Y();
    auto z = pos.Z();
    return ( abs(x) < 874.51
             && abs(y-55.) < 874.51
             && 136.875 < z && z < 447.375 );
}

T2K::RecoVec tracksinfiducialvolume(const T2K::RecoVec& v) {
    return Filter(v, [](T2K::Reco* r) { return infiducialvolume(r->FrontPosition);} );
}
"""

def parsecml():
    parser = ArgumentParser()
    parser.add_argument("inputfiles", nargs="+", help="List of input ROOT files.")
    parser.add_argument("--enable-truth", action="store_true", help="Swtich on truth information")
    return parser.parse_args()

def tostdvector(l):
    result = ROOT.std.vector(ROOT.std.string)()
    for x in l:
        result.push_back(x)
    return result

def main():
    args = parsecml()
    # Load libraries
    ROOT.gSystem.Load("libT2KDataFrameMakeProject.so")
    ROOT.gSystem.Load("libT2KDataSource.so") # or .dylib on OS X
    ROOT.T2K.T2KDataSource
    ROOT.gInterpreter.Declare(FIDUCIALVOLUMECODE)
    # Run analysis
    rdf = ROOT.T2K.MakeT2KDataFrame(tostdvector(args.inputfiles), args.enable_truth)
    analysis = (rdf.Define("tracksinfv", "tracksinfiducialvolume(t2kreco)")
    # Select negative tracks
    .Define("muons", "return Filter(tracksinfv, [](T2K::Reco* r){return r->Charge<0;});")
    .Define("mup", "return Map(muons, [](T2K::Reco* r){return r->FrontMomentum;});")
    .Define("mutheta", "return Map(muons, [](T2K::Reco* r) { return r->FrontDirection.Theta(); } );"))
    # Make a histogram
    hist = analysis.Histo1D(("muonmomentum", ";p_#mu [MeV/c]", 100, 0, 2000), "mup")
    # Make an ntuple
    analysis.Snapshot("tree", "examplentuple.root", tostdvector(["mup", "mutheta"]))
    # Save the histogram we made earlier
    hist.GetValue().SaveAs("examplehistogram.root");
    return


if __name__ == "__main__":
    main()
