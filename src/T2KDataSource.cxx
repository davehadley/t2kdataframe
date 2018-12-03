#include <T2KDataSource.hxx>
#include <ROOT/TSeq.hxx>
#include <TClass.h>
#include <TROOT.h>         // For the gROOTMutex
#include <TVirtualMutex.h> // For the R__LOCKGUARD
#include <ROOT/RMakeUnique.hxx>

#include <algorithm>
#include <vector>

namespace T2K {

    std::vector<void *> MyRRootDS::GetColumnReadersImpl(std::string_view name, const std::type_info &id) {
        const auto colTypeName = GetTypeName(name);
        const auto &colTypeId = ROOT::Internal::RDF::TypeName2TypeID(colTypeName);
        if (id != colTypeId) {
            std::string err = "The type of column \"";
            err += name;
            err += "\" is ";
            err += colTypeName;
            err += " but a different one has been selected.";
            throw std::runtime_error(err);
        }

        const auto index =
                std::distance(fListOfBranches.begin(), std::find(fListOfBranches.begin(), fListOfBranches.end(), name));
        std::vector<void *> ret(fNSlots);
        for (auto slot : ROOT::TSeqU(fNSlots)) {
            ret[slot] = (void *) &fBranchAddresses[index][slot];
        }
        return ret;
    }

    MyRRootDS::MyRRootDS(std::string_view treeName, std::vector<std::string> fileNameGlob)
            : fTreeName(treeName), fFileNameGlob(fileNameGlob), fModelChain(std::string(treeName).c_str()) {
        for (auto &f : fFileNameGlob) { fModelChain.Add(f.c_str()); }

        const TObjArray &lob = *fModelChain.GetListOfBranches();
        //fListOfBranches.resize(lob.GetEntries());

        //TIterCategory<TObjArray> iter(&lob);
        //std::transform(iter.Begin(), iter.End(), fListOfBranches.begin(), [](TObject *o) { return o->GetName(); });
        for (int i = 0; i < lob.GetEntries(); ++i) {
            auto obj = lob.At(i);
            auto c = std::string(obj->GetName());
            if (c == "CTMTriggerPattern" || c == "SoftwareVersion" || c == "GeometryHash" || c == "AlignmentId") {
                std::cout << "Skipping" << c << std::endl;
                continue;
            }
            fListOfBranches.push_back(c);
        }
    }

    MyRRootDS::~MyRRootDS() {
        for (auto addr : fAddressesToFree) {
            delete addr;
        }
    }

    std::string MyRRootDS::GetTypeName(std::string_view colName) const {
        if (!HasColumn(colName)) {
            std::string e = "The dataset does not have column ";
            e += colName;
            throw std::runtime_error(e);
        }
        // TODO: we need to factor out the routine for the branch alone...
        // Maybe a cache for the names?
        auto typeName =
                ROOT::Internal::RDF::ColumnName2ColumnTypeName(std::string(colName), /*nsID=*/0, &fModelChain, /*ds=*/
                                                               nullptr,
                        /*isCustomCol=*/false);
        // We may not have yet loaded the library where the dictionary of this type is
        TClass::GetClass(typeName.c_str());
        return typeName;
    }

    const std::vector<std::string> &MyRRootDS::GetColumnNames() const {
        return fListOfBranches;
    }

    bool MyRRootDS::HasColumn(std::string_view colName) const {
        if (!fListOfBranches.empty())
            GetColumnNames();
        return fListOfBranches.end() != std::find(fListOfBranches.begin(), fListOfBranches.end(), colName);
    }

    void MyRRootDS::InitSlot(unsigned int slot, ULong64_t firstEntry) {
        auto chain = new TChain(fTreeName.c_str());
        chain->ResetBit(kMustCleanup);
        for (auto &f : fFileNameGlob) { chain->Add(f.c_str()); }
        TString setBranches;
        for (auto i : ROOT::TSeqU(fListOfBranches.size())) {
            auto colName = fListOfBranches[i].c_str();
            auto &addr = fBranchAddresses[i][slot];
            auto typeName = GetTypeName(colName);
            auto typeClass = TClass::GetClass(typeName.c_str());
            if (typeClass) {
                chain->SetBranchAddress(colName, &addr, nullptr, typeClass, EDataType(0), true);
            } else {
                if (!addr) {
                    addr = new double();
                    fAddressesToFree.emplace_back((double *) addr);
                }
                chain->SetBranchAddress(colName, addr);
            }
        }
        // try to increase cache size to 100MB for imporved IO
        chain->SetCacheSize(100000000);
        chain->AddBranchToCache("*", true);
        //chain->SetCacheEntryRange(firstEntry,chain->GetEntries());
        //chain->StopCacheLearningPhase();
        chain->GetEntry(firstEntry);
        fChains[slot].reset(chain);

    }

    void MyRRootDS::FinaliseSlot(unsigned int slot) {
#ifdef T2K_DEBUG_PERFORMANCE
        fChains[slot]->GetTree()->PrintCacheStats();
#endif
        fChains[slot].reset(nullptr);

    }

    std::vector<std::pair<ULong64_t, ULong64_t>> MyRRootDS::GetEntryRanges() {
        auto entryRanges(std::move(fEntryRanges)); // empty fEntryRanges
        return entryRanges;
    }

    bool MyRRootDS::SetEntry(unsigned int slot, ULong64_t entry) {
        if (entry != slotcurrententry[slot]) {
            fChains[slot]->GetEntry(entry);
            slotcurrententry[slot] = entry;
        }
        return true;
    }

    void MyRRootDS::SetNSlots(unsigned int nSlots) {
        assert(0U == fNSlots && "Setting the number of slots even if the number of slots is different from zero.");

        fNSlots = nSlots;
        slotcurrententry = std::vector<ULong64_t>(nSlots, -1);
        const auto nColumns = fListOfBranches.size();
        // Initialise the entire set of addresses
        fBranchAddresses.resize(nColumns, std::vector<void *>(fNSlots, nullptr));

        fChains.resize(fNSlots);
    }

    void MyRRootDS::Initialise() {
        const auto nentries = fModelChain.GetEntries();
        const auto chunkSize = nentries / fNSlots;
        const auto reminder = 1U == fNSlots ? 0 : nentries % fNSlots;
        auto start = 0UL;
        auto end = 0UL;
        for (auto i : ROOT::TSeqU(fNSlots)) {
            start = end;
            end += chunkSize;
            fEntryRanges.emplace_back(start, end);
            (void) i;
        }
        fEntryRanges.back().second += reminder;
    }

    ROOT::RDataFrame MakeRootDataFrame(std::string_view treeName, std::vector<std::string> fileNameGlob) {
        ROOT::RDataFrame tdf(std::make_unique<MyRRootDS>(treeName, fileNameGlob));
        return tdf;
    }

    T2KDataFrameWrapper MakeT2KDataFrame(std::vector<std::string> fileNameGlob,
                                                bool enabletruth,
                                                std::vector<std::string> treeNames,
                                                BunchTiming bunchTiming) {
        if(!enabletruth) {
            // Remove truth trees
            std::vector<std::string> t;
            std::copy_if(treeNames.begin(), treeNames.end(), std::back_inserter(t), [](std::string& s){ return s.find("TruthDir")==std::string::npos; });
            treeNames = t;
        }
        // Create RDataFrame
        auto nbunches = bunchTiming.size();
        auto datasource = std::make_unique<T2KDataSource>(treeNames, fileNameGlob, nbunches);
        auto datasourceptr = datasource.get(); // yes, I know this is super-dangerous
        ROOT::RDataFrame tdf(std::move(datasource));
        // Apply bunching condition
        auto reco = tdf.Define("t2kbunch", [nbunches](ULong64_t entry) -> int {return entry % nbunches;}, {"tdfentry_"})
                .Define("t2krunnumber", [](int runid) {return runid;}, {"Global_RunID"})
                .Define("t2kbunchtmin", [bunchTiming](int bunch, int run){ return bunchTiming.low(bunch, run); }, {"t2kbunch", "t2krunnumber"})
                .Define("t2kbunchtmax", [bunchTiming](int bunch, int run){ return bunchTiming.high(bunch, run); }, {"t2kbunch", "t2krunnumber"})
                .Define("t2kbunchcenter", "(t2kbunchtmin+t2kbunchtmax)/2.0")
                .Define("t2kallbunchcentres", [bunchTiming](int bunch, int run){ return bunchTiming.allbunchcentres(run); }, {"t2kbunch", "t2krunnumber"})
                        // Apply bunching to the reconstruction
                .Define("t2krecoraw", T2K::tcatorvec<ND::TGlobalReconModule::TGlobalPID*>, {"Global_PIDs"})
                .Define("t2krecorawtime", "return ROOT::VecOps::Map(t2krecoraw, [](ND::TGlobalReconModule::TGlobalPID* obj) -> double { return obj->FrontPosition.T(); });")
                        //.Define("t2kiscosmic", [](int a, int b){ return a||b; }, {"BasicHeader_TripTCosmicEvent", "BasicHeader_FGDCosmicEvent"})
                .Define("t2kreco", "t2krecoraw[t2kbunchtmin < t2krecorawtime && t2krecorawtime < t2kbunchtmax]")
                .Define("t2kisnewfile", [datasourceptr, nbunches](unsigned int slot, ULong64_t entry, int t2kbunch) { return datasourceptr->GetSlotChain(slot, 0).GetChainEntryNumber(0)==(entry/nbunches) && t2kbunch==0; }, {"tdfslot_", "tdfentry_", "t2kbunch"} );
        if(enabletruth) {
            //Add Truth vertices
            reco =
#ifdef T2K_USE_GENIE
                    reco.Define("t2kroovtx", T2K::tcatorvec<ND::GRooTrackerVtx*>, {"GRooTrackerVtx_Vtx"})
#else
                    reco.Define("t2kroovtx", T2K::tcatorvec<ND::GRooTrackerVtx*>, {"NRooTrackerVtx_Vtx"})
#endif
                    .Define("t2ktruthtraj", T2K::tcatorvec<TruthTraj*>, {"Trajectories_Trajectories"})
                    .Define("t2ktruthvtx", T2K::tcatorvec<TruthVtx*>, {"Vertices_Vertices"});
        }
        else {
            // No truth, fill with empty containers
            reco = reco.Define("t2kroovtx", [](){ return T2K::TruthRooVtxVec(); }, {})
                    .Define("t2ktruthtraj", []() { return T2K::TruthTrajVec(); }, {})
                    .Define("t2ktruthvtx", []() { return T2K::TruthVtxVec(); }, {});
        }
        // Wrap and return result
        return T2KDataFrameWrapper(std::move(tdf), std::move(reco));
    }

}

