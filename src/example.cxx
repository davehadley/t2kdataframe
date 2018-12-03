#include "T2KDataSource.hxx"

using namespace T2K;

bool infiducialvolume(const TLorentzVector& pos) {
    auto x = pos.X();
    auto y = pos.Y();
    auto z = pos.Z();
    return ( abs(x) < 874.51
             && abs(y-55.) < 874.51
             && 136.875 < z && z < 447.375 );
}

RecoVec tracksinfiducialvolume(const RecoVec& v) {
    return Filter(v, [](Reco* r) { return infiducialvolume(r->FrontPosition);} );
}

int main(int argc, char* argv[]) {
    // Read in command line arguments
    // First argument is 1 or 0 to enable/disable truth
    bool enabletruth = std::string(argv[1])!="0";
    // Other arguments are a list of files to process
    std::vector<std::string> filelist;
    for (int i = 2; i < argc; ++i) {
        filelist.emplace_back(argv[i]);
    }

    // Switch on paralell processing
    //ROOT::EnableImplicitMT();

    auto rdf = MakeT2KDataFrame(filelist, enabletruth);

    auto analysis = rdf->Define("tracksinfv", tracksinfiducialvolume, {"t2kreco"})
                    // Select negative tracks
            .Define("muons", [](RecoVec& v) { return Filter(v, [](Reco* r){return r->Charge<0;}); }, {"tracksinfv"})
            .Define("mup", [](RecoVec& v) { return Map(v, [](Reco* r) { return r->FrontMomentum; } ); }, {"muons"} )
            .Define("mutheta", [](RecoVec& v) { return Map(v, [](Reco* r) { return r->FrontDirection.Theta(); } ); }, {"muons"} );
    // Make a histogram
    auto hist = analysis.Histo1D({"muonmomentum", ";p_#mu [MeV/c]", 100, 0, 2000}, "mup");
    // Make an ntuple
    analysis.Snapshot("tree", "examplentuple.root", {"mup", "mutheta"});
    // Save the histogram we made earlier
    hist->SaveAs("examplehistogram.root");
}