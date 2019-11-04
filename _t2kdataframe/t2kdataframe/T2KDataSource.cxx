#include <t2kdataframe/T2KDataSource.hxx>
#include <ROOT/TSeq.hxx>
#include <TClass.h>
#include <TROOT.h>         // For the gROOTMutex
#include <TVirtualMutex.h> // For the R__LOCKGUARD
#include <ROOT/RMakeUnique.hxx>

#include <algorithm>
#include <vector>

namespace T2K {

    RooVtxWrapper::RooVtxWrapper(ND::RooTrackerVtxBase* vtx) : ptr(vtx), genie(nullptr), neut(nullptr) {
        if(vtx) {
            genie = dynamic_cast<ND::GRooTrackerVtx *>(vtx);
            neut = dynamic_cast<ND::NRooTrackerVtx *>(vtx);
        }
    }
    ND::RooTrackerVtxBase* RooVtxWrapper::getptr() { return ptr;}
    TObjString* RooVtxWrapper::EvtCode() { return genie?genie->EvtCode:neut->EvtCode; }
    int         RooVtxWrapper::EvtNum() const { return genie?genie->EvtNum:neut->EvtNum; }
    double      RooVtxWrapper::EvtXSec() const { return genie?genie->EvtXSec:neut->EvtXSec; }
    double      RooVtxWrapper::EvtDXSec() const { return genie?genie->EvtDXSec:neut->EvtDXSec; }
    double      RooVtxWrapper::EvtWght() const { return genie?genie->EvtWght:neut->EvtWght; }
    double      RooVtxWrapper::EvtProb() const { return genie?genie->EvtProb:neut->EvtProb; }
    double*      RooVtxWrapper::EvtVtx() { return genie?genie->EvtVtx:neut->EvtVtx; }
    int         RooVtxWrapper::StdHepN() const { return genie?genie->StdHepN:neut->StdHepN; }
    Int_t      *RooVtxWrapper::StdHepPdg() { return genie?genie->StdHepPdg:neut->StdHepPdg; }
    Int_t      *RooVtxWrapper::StdHepStatus() { return genie?genie->StdHepStatus:neut->StdHepStatus; }
    double     (*RooVtxWrapper::StdHepX4())[4] { return genie?genie->StdHepX4:neut->StdHepX4; }
    double     (*RooVtxWrapper::StdHepP4())[4] { return genie?genie->StdHepP4:neut->StdHepP4; }
    double     (*RooVtxWrapper::StdHepPolz())[3] { return genie?genie->StdHepPolz:neut->StdHepPolz; }
    Int_t      *RooVtxWrapper::StdHepFd() { return genie?genie->StdHepFd:neut->StdHepFd; }
    Int_t      *RooVtxWrapper::StdHepLd() { return genie?genie->StdHepLd:neut->StdHepLd; }
    Int_t      *RooVtxWrapper::StdHepFm() { return genie?genie->StdHepFm:neut->StdHepFm; }
    Int_t      *RooVtxWrapper::StdHepLm() { return genie?genie->StdHepLm:neut->StdHepLm; }
    int         RooVtxWrapper::G2NeutEvtCode() const { return genie?genie->G2NeutEvtCode:neut->EvtCode->GetString().Atoi(); }
    TObjString* RooVtxWrapper::GeomPath() { return genie?genie->GeomPath:neut->GeomPath; }
    TObjString* RooVtxWrapper::GeneratorName() { return genie?genie->GeneratorName:neut->GeneratorName; }
    TObjString* RooVtxWrapper::OrigFileName() { return genie?genie->OrigFileName:neut->OrigFileName; }
    TObjString* RooVtxWrapper::OrigTreeName() { return genie?genie->OrigTreeName:neut->OrigTreeName; }
    int         RooVtxWrapper::OrigEvtNum() const { return genie?genie->OrigEvtNum:neut->OrigEvtNum; }
    int         RooVtxWrapper::OrigTreeEntries() const { return genie?genie->OrigTreeEntries:neut->OrigTreeEntries; }
    double      RooVtxWrapper::OrigTreePOT() const { return genie?genie->OrigTreePOT:neut->OrigTreePOT; }
    double      RooVtxWrapper::TimeInSpill() const { return genie?genie->TimeInSpill:neut->TimeInSpill; }
    int         RooVtxWrapper::TruthVertexID() const { return genie?genie->TruthVertexID:neut->TruthVertexID; }
    int RooVtxWrapper::NeutrinoPdg() const { return genie?genie->StdHepPdg[0]:neut->StdHepPdg[0]; }
    double RooVtxWrapper::NeutrinoEnergy() const { return genie?genie->StdHepP4[0][3]:neut->StdHepP4[0][3]; }

    ROOT::VecOps::RVec<RooVtxWrapper> RooVtxWrapper::wrapclonesarray(const TClonesArray& arr) {
        ROOT::VecOps::RVec<RooVtxWrapper> v;
        for (int i = 0; i < arr.GetEntries(); ++i) {
            TObject *obj = arr.At(i);
            ND::RooTrackerVtxBase* newobj = dynamic_cast<ND::RooTrackerVtxBase*>(obj);
            v.emplace_back(RooVtxWrapper(newobj));
        }
        return std::move(v);
    }

    std::vector<void *> T2KDataSource::T2KSingleRootChainDS::GetColumnReadersImpl(std::string_view name, const std::type_info &id) {
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

    T2KDataSource::T2KSingleRootChainDS::T2KSingleRootChainDS(std::string_view treeName, std::vector<std::string> fileNameGlob)
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

    T2KDataSource::T2KSingleRootChainDS::~T2KSingleRootChainDS() {
        for (auto addr : fAddressesToFree) {
            delete addr;
        }
    }

    std::string T2KDataSource::T2KSingleRootChainDS::GetTypeName(std::string_view colName) const {
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

    const std::vector<std::string> &T2KDataSource::T2KSingleRootChainDS::GetColumnNames() const {
        return fListOfBranches;
    }

    bool T2KDataSource::T2KSingleRootChainDS::HasColumn(std::string_view colName) const {
        if (!fListOfBranches.empty())
            GetColumnNames();
        return fListOfBranches.end() != std::find(fListOfBranches.begin(), fListOfBranches.end(), colName);
    }

    void T2KDataSource::T2KSingleRootChainDS::InitSlot(unsigned int slot, ULong64_t firstEntry) {
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

    void T2KDataSource::T2KSingleRootChainDS::FinaliseSlot(unsigned int slot) {
#ifdef T2K_DEBUG_PERFORMANCE
        fChains[slot]->GetTree()->PrintCacheStats();
#endif
        fChains[slot].reset(nullptr);

    }

    std::vector<std::pair<ULong64_t, ULong64_t>> T2KDataSource::T2KSingleRootChainDS::GetEntryRanges() {
        auto entryRanges(std::move(fEntryRanges)); // empty fEntryRanges
        return entryRanges;
    }

    bool T2KDataSource::T2KSingleRootChainDS::SetEntry(unsigned int slot, ULong64_t entry) {
        if (entry != slotcurrententry[slot]) {
            fChains[slot]->GetEntry(entry);
            slotcurrententry[slot] = entry;
        }
        return true;
    }

    void T2KDataSource::T2KSingleRootChainDS::SetNSlots(unsigned int nSlots) {
        assert(0U == fNSlots && "Setting the number of slots even if the number of slots is different from zero.");

        fNSlots = nSlots;
        slotcurrententry = std::vector<ULong64_t>(nSlots, -1);
        const auto nColumns = fListOfBranches.size();
        // Initialise the entire set of addresses
        fBranchAddresses.resize(nColumns, std::vector<void *>(fNSlots, nullptr));

        fChains.resize(fNSlots);
    }

    void T2KDataSource::T2KSingleRootChainDS::Initialise() {
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

    T2KDataFrameWrapper MakeT2KDataFrame(std::vector<std::string> fileNameGlob,
                                         bool enabletruth,
                                         std::vector<std::string> treeNames,
                                         BunchTiming bunchTiming,
                                         bool truthisgenie
                                         ) {
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
            if (truthisgenie) {
                reco = reco.Define("t2kroovtx", T2K::RooVtxWrapper::wrapclonesarray, {"GRooTrackerVtx_Vtx"});
            } else {
                reco = reco.Define("t2kroovtx", T2K::RooVtxWrapper::wrapclonesarray, {"NRooTrackerVtx_Vtx"});
            }
            reco = reco.Define("t2ktruthtraj", T2K::tcatorvec<TruthTraj*>, {"Trajectories_Trajectories"})
                    .Define("t2ktruthvtx", T2K::tcatorvec<TruthVtx*>, {"Vertices_Vertices"});
        }
        else {
            // No truth, fill with empty containers
            reco = reco.Define("t2kroovtx", [](){ return T2K::TruthRooVtxVec(); }, {})
                    .Define("t2ktruthtraj", []() { return T2K::TruthTrajVec(); }, {})
                    .Define("t2ktruthvtx", []() { return T2K::TruthVtxVec(); }, {})
                    ;
        }
        // Wrap and return result
        return T2KDataFrameWrapper(std::move(tdf), std::move(reco));
    }


    ULong64_t T2KDataSource::GetEntries() { return _chainreaders.at(0)->GetEntries() * nbunches; }

    int T2KDataSource::GetNSlots() { return _nslots; }

    TChain& T2KDataSource::GetSlotChain(int slot, int treeindex) { return _chainreaders.at(treeindex)->GetSlotChain(slot); }

    std::vector<void *> T2KDataSource::GetColumnReadersImpl(std::string_view name, const std::type_info &t) {
        return _chainreaders.at(_columnnamemap.at(std::string(name)))->GetColumnReadersImpl(getBranchName(name), t);
    }


    T2KDataSource::T2KDataSource(std::vector<std::string> treenames_, std::vector<std::string> filename, int nbunches_) :
            treenames(treenames_),
            processed(0), nbunches(nbunches_) {
        for (auto &tn : treenames) {
            _chainreaders.push_back(make_unique<T2K::T2KDataSource::T2KSingleRootChainDS>(tn, filename));
        }
    }

    std::string T2KDataSource::GetTypeName(std::string_view colName) const {
        return _chainreaders.at(_columnnamemap.at(std::string(colName)))->GetTypeName(getBranchName(colName));
    }

    const std::vector<std::string>& T2KDataSource::GetColumnNames() const {
        if (columnnames.size() == 0) {
            for (int i = 0; i < treenames.size(); ++i) {
                auto tn = treenames.at(i);
                auto &chr = _chainreaders.at(i);
                for (auto &c : chr->GetColumnNames()) {
                    auto s = tn.substr(tn.find("/") + 1);
                    if (c == "CTMTriggerPattern" || c == "SoftwareVersion" || c == "GeometryHash" ||
                        c == "AlignmentId") {
                        std::cout << "Skipping" << c << std::endl;
                        continue;
                    }
                    std::string colname = s + "_" + c;
                    columnnames.push_back(colname);
                    if (_columnnamemap.insert({colname, i}).second == false) {
                        // insertion failed! Duplicate column names is not allowed
                        throw std::runtime_error("Duplicate column names " + colname + " is not allowed");
                    }

                }
            }
        }
        return columnnames;
    }

    std::string T2KDataSource::getBranchName(std::string colName) const {
        return colName.substr(colName.find("_") + 1);
    }

    std::string T2KDataSource::getBranchName(std::string_view colName) const {
        return getBranchName(std::string(colName));
    }

    bool T2KDataSource::HasColumn(std::string_view colName) const {
        for (auto &n : columnnames) {
            if (n == colName) {
                return true;
            }
        }
        return false;
    }

    void T2KDataSource::InitSlot(unsigned int slot, ULong64_t firstEntry) {
        for (auto &c : _chainreaders) { c->InitSlot(slot, firstEntry); }
    }

    void T2KDataSource::Initialise() {
        for (auto &c : _chainreaders) { c->Initialise(); }
    }

    std::vector<std::pair<ULong64_t, ULong64_t> > T2KDataSource::GetEntryRanges() {
        for (auto &c : _chainreaders) { c->GetEntryRanges(); }
        std::vector<std::pair<ULong64_t, ULong64_t> > result;
        //auto nentries = nbunches*1000;//_chainreaders.front()->GetEntries();
        auto nentries = nbunches * _chainreaders.front()->GetEntries();
        if (processed < nentries) {
            auto chunk = nentries / _nslots;
            auto start = processed;
            auto end = start;
            for (unsigned int i = 0; i < _nslots; ++i) {
                start = end;
                end += start + chunk;
                result.push_back({start, end});
            }
            result.back().second = nentries;
            processed = result.back().second;
        }
        return result;
    }

    bool T2KDataSource::SetEntry(unsigned int slot, ULong64_t entry) {
        bool result = true;
        for (auto &reader : _chainreaders) {
            result = result && reader->SetEntry(slot, GetEntryNumber(entry));
        }
        return result;
    }

    void T2KDataSource::SetNSlots(unsigned int nSlots) {
        _nslots = nSlots;
        for (auto &chr : _chainreaders) {
            chr->SetNSlots(nSlots);
        }
    }

    void T2KDataSource::Finalise() {
        for (auto &c : _chainreaders) { c->Finalise(); }
    }

    void T2KDataSource::FinaliseSlot(unsigned int i) { for (auto &c : _chainreaders) { c->FinaliseSlot(i); }}

    ULong64_t T2KDataSource::GetEntryNumber(ULong64_t entry) { return entry / this->nbunches; }

    ULong64_t T2KDataSource::GetBunchNumber(ULong64_t entry) { return entry % this->nbunches; }

    BunchTiming::BunchTiming(ROOT::VecOps::RVec<double> centres, double window) {
        _data.emplace_back(Record{numeric_limits<int>::min(), numeric_limits<int>::max(), centres, window});
    }

    BunchTiming::BunchTiming(ROOT::VecOps::RVec<double> lowedges, ROOT::VecOps::RVec<double> highedges) {
        _data.emplace_back(Record{numeric_limits<int>::min(), numeric_limits<int>::max(), lowedges, highedges});
    }

    BunchTiming::BunchTiming(double window) {
        // Taken from:
        //  https://repo.nd280.org/viewvc/ND280/psyche/psycheND280Utils/data/BunchPosition.dat?view=log
        _data = ROOT::VecOps::RVec<Record>{
                // moved below: {-1,      -1,   {2750.2,  3332.0,  3914.7,  4497.0,  5078.4,  5659.7,  6243.4,  6824.2}, window},
                {0,                          6000,                       {2839.7,  3423.5,  4005.4,  4588.6,  5172.2,  5754.6,  -1000.0, -1000.0}, window},
                {6000,                       7000,                       {2853.95, 3444.15, 4030.41, 4620.34, 5180.28, 5770.12, 6343.77, 6924.67}, window},
                {7000,                       8000,                       {3019.11, 3597.74, 4180.73, 4763.93, 5346.49, 5927.83, 6508.5,  7093.56}, window},
                {8000,                       8550,                       {3024.22, 3606.11, 4188.01, 4769.9,  5351.79, 5933.68, 6515.58, 7097.47}, window},
                {8550,                       8800,                       {3013.55, 3597.55, 4178.24, 4758.78, 5338.21, 5927.31, 6505.71, 7086.80}, window},
                {8895,                       9300,                       {3014.69, 3600.57, 4178.89, 4764.24, 5342.35, 5931.73, 6506.37, 7090.15}, window},
                {9301,                       9999,                       {3014.69, 3600.57, 4178.89, 4764.24, 5342.35, 5931.73, 6506.37, 7090.15}, window},
                {10251,                      10522,                      {3013.71, 3595.50, 4178.37, 4764.48, 5342.07, 5922.57, 6509.16, 7086.23}, window},
                {10826,                      12063,                      {3018.20, 3599.10, 4184.62, 4769.22, 5345.73, 5932.20, 6511.85, 7096.31}, window},
                {12064,                      99999,                      {2994.20, 3581.27, 4165.74, 4752.16, 5335.15, 5917.88, 6502.08, 7085.51}, window},
                {numeric_limits<int>::min(), numeric_limits<int>::max(), {2750.2,  3332.0,  3914.7,  4497.0,  5078.4,  5659.7,  6243.4,  6824.2},  window},

        };
    }

    int BunchTiming::size() const { return _data.front().timelow.size(); }

    double BunchTiming::low(int bunch, int run) const {
        for (const Record &r : this->_data) {
            if (r.runlow <= run && run < r.runhigh) {
                return r.timelow.at(bunch);
            }
        }
        throw std::runtime_error("unknown run number");
    }

    double BunchTiming::high(int bunch, int run) const {
        for (const Record &r : this->_data) {
            if (r.runlow <= run && run < r.runhigh) {
                return r.timehigh.at(bunch);
            }
        }
        throw std::runtime_error("unknown run number");
    }

    double BunchTiming::center(int bunch, int run) const { return (low(bunch, run) + high(bunch, run)) / 2.0; }

    ROOT::VecOps::RVec<double> BunchTiming::allbunchcentres(int run) const {
        for (const Record &r : this->_data) {
            if (r.runlow <= run && run < r.runhigh) {
                return r.timecenter;
            }
        }
        throw std::runtime_error("unknown run number");
    }


}

