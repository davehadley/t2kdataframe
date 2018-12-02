
#ifndef T2KDATAFRAME
#define T2KDATAFRAME

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDataSource.hxx"
//#include "ROOT/RDFUtils.hxx"
#include "ROOT/RVec.hxx"
#include <TChain.h>

using namespace std;
#include "ND__TGlobalReconModule.h"
#include "ND__TTruthTrajectoriesModule.h"
#include "ND__TTruthVerticesModule.h"
#include "ND__GRooTrackerVtx.h"

#include <memory>
#include <map>

namespace T2K {

    //Make some type aliases
    using Reco = ND::TGlobalReconModule::TGlobalPID;
    using RecoVec = ROOT::VecOps::RVec<Reco *>;
    using TruthTraj = ND::TTruthTrajectoriesModule::TTruthTrajectory;
    using TruthTrajVec = ROOT::VecOps::RVec<TruthTraj *>;
    using TruthVtx = ND::TTruthVerticesModule::TTruthVertex;
    using TruthVtxVec = ROOT::VecOps::RVec<TruthVtx *>;
    using TruthRooVtx = ND::GRooTrackerVtx;
    using TruthRooVtxVec = ROOT::VecOps::RVec<TruthRooVtx *>;

    class MyRRootDS final : public ROOT::RDF::RDataSource {
    private:
        unsigned int fNSlots = 0U;
        std::vector<ULong64_t> slotcurrententry;
        std::string fTreeName;
        std::vector<std::string> fFileNameGlob;
        mutable TChain fModelChain; // Mutable needed for getting the column type name
        std::vector<double *> fAddressesToFree;
        std::vector<std::string> fListOfBranches;
        std::vector<std::pair<ULong64_t, ULong64_t>> fEntryRanges;
        std::vector<std::vector<void *>> fBranchAddresses; // first container-> slot, second -> column;
        std::vector<std::unique_ptr<TChain>> fChains;

    public:
        std::vector<void *> GetColumnReadersImpl(std::string_view, const std::type_info &);

        ULong64_t GetEntries() { return fModelChain.GetEntries(); }

    public:
        MyRRootDS(std::string_view treeName, std::vector<std::string> fileNameGlob);

        ~MyRRootDS();

        std::string GetTypeName(std::string_view colName) const;

        const std::vector<std::string> &GetColumnNames() const;

        bool HasColumn(std::string_view colName) const;

        void InitSlot(unsigned int slot, ULong64_t firstEntry);

        void FinaliseSlot(unsigned int slot);

        std::vector<std::pair<ULong64_t, ULong64_t>> GetEntryRanges();

        bool SetEntry(unsigned int slot, ULong64_t entry);

        void SetNSlots(unsigned int nSlots);

        void Initialise();

        TChain &GetSlotChain(int slot) { return *(fChains.at(slot)); }
    };

    ROOT::RDataFrame MakeRootDataFrame(std::string_view treeName, std::vector<std::string> fileNameGlob);

    class T2KDataSource final : public ROOT::RDF::RDataSource {
    private:
        int _nslots;
        mutable std::vector<std::string> columnnames;
        std::vector<std::string> treenames;
        ULong64_t processed;
        int nbunches;
        std::vector<unique_ptr<T2K::MyRRootDS> > _chainreaders;
        mutable std::map<std::string, int> _columnnamemap;

    public:

        ULong64_t GetEntries() { return _chainreaders.at(0)->GetEntries() * nbunches; }

        int GetNSlots() { return _nslots; }

        TChain &GetSlotChain(int slot, int treeindex) { return _chainreaders.at(treeindex)->GetSlotChain(slot); }

        std::vector<void *> GetColumnReadersImpl(std::string_view name, const std::type_info &t) {
            return _chainreaders.at(_columnnamemap.at(std::string(name)))->GetColumnReadersImpl(getBranchName(name), t);
        }


        T2KDataSource(std::vector<std::string> treenames_, std::vector<std::string> filename, int nbunches_) :
                treenames(treenames_),
                processed(0), nbunches(nbunches_) {
            for (auto &tn : treenames) {
                _chainreaders.push_back(make_unique<T2K::MyRRootDS>(tn, filename));
            }
        }

        std::string GetTypeName(std::string_view colName) const {
            return _chainreaders.at(_columnnamemap.at(std::string(colName)))->GetTypeName(getBranchName(colName));
        }

        const std::vector<std::string> &GetColumnNames() const {
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

        std::string getBranchName(std::string colName) const {
            return colName.substr(colName.find("_") + 1);
        }

        std::string getBranchName(std::string_view colName) const {
            return getBranchName(std::string(colName));
        }

        bool HasColumn(std::string_view colName) const {
            for (auto &n : columnnames) {
                if (n == colName) {
                    return true;
                }
            }
            return false;
        }

        void InitSlot(unsigned int slot, ULong64_t firstEntry) {
            for (auto &c : _chainreaders) { c->InitSlot(slot, firstEntry); }
        }

        virtual void Initialise() {
            for (auto &c : _chainreaders) { c->Initialise(); }
        }

        std::vector<std::pair<ULong64_t, ULong64_t> > GetEntryRanges() {
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

        bool SetEntry(unsigned int slot, ULong64_t entry) {
            bool result = true;
            for (auto &reader : _chainreaders) {
                result = result && reader->SetEntry(slot, GetEntryNumber(entry));
            }
            return result;
        }

        void SetNSlots(unsigned int nSlots) {
            _nslots = nSlots;
            for (auto &chr : _chainreaders) {
                chr->SetNSlots(nSlots);
            }
        }

        void Finalise() {
            for (auto &c : _chainreaders) { c->Finalise(); }
        }

        void FinaliseSlot(unsigned int i) { for (auto &c : _chainreaders) { c->FinaliseSlot(i); }}

    private:
        ULong64_t GetEntryNumber(ULong64_t entry) { return entry / this->nbunches; }

        ULong64_t GetBunchNumber(ULong64_t entry) { return entry % this->nbunches; }

    };

    class BunchTiming {
    private:
        struct Record {
            int runlow;
            int runhigh;
            ROOT::VecOps::RVec<double> timelow;
            ROOT::VecOps::RVec<double> timehigh;
            ROOT::VecOps::RVec<double> timecenter;

            Record(int rl, int rh, ROOT::VecOps::RVec<double> tl, ROOT::VecOps::RVec<double> th) : runlow(rl),
                                                                                                   runhigh(rh),
                                                                                                   timelow(tl),
                                                                                                   timehigh(th),
                                                                                                   timecenter((timelow +
                                                                                                               timehigh) /
                                                                                                              2.0) {}

            Record(int rl, int rh, ROOT::VecOps::RVec<double> centres, double window) : runlow(rl), runhigh(rh) {
                timelow = centres - window / 2.0;
                timehigh = centres + window / 2.0;
                timecenter = (timelow + timehigh) / 2.0;
            }
        };

        ROOT::VecOps::RVec<Record> _data;
    public:
        BunchTiming(ROOT::VecOps::RVec<double> centres, double window) {
            _data.emplace_back(Record{numeric_limits<int>::min(), numeric_limits<int>::max(), centres, window});
        }

        BunchTiming(ROOT::VecOps::RVec<double> lowedges, ROOT::VecOps::RVec<double> highedges) {
            _data.emplace_back(Record{numeric_limits<int>::min(), numeric_limits<int>::max(), lowedges, highedges});
        }

        BunchTiming(double window = 150.0) {
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

        int size() const { return _data.front().timelow.size(); }

        double low(int bunch, int run) const {
            for (const Record &r : this->_data) {
                if (r.runlow <= run && run < r.runhigh) {
                    return r.timelow.at(bunch);
                }
            }
            throw std::runtime_error("unknown run number");
        }

        double high(int bunch, int run) const {
            for (const Record &r : this->_data) {
                if (r.runlow <= run && run < r.runhigh) {
                    return r.timehigh.at(bunch);
                }
            }
            throw std::runtime_error("unknown run number");
        }

        double center(int bunch, int run) const { return (low(bunch, run) + high(bunch, run)) / 2.0; }

        ROOT::VecOps::RVec<double> allbunchcentres(int run) const {
            for (const Record &r : this->_data) {
                if (r.runlow <= run && run < r.runhigh) {
                    return r.timecenter;
                }
            }
            throw std::runtime_error("unknown run number");
        }
    };

    template<class T>
    ROOT::VecOps::RVec<T> tcatorvec(const TClonesArray &arr) {
        ROOT::VecOps::RVec<T> v;
        for (int i = 0; i < arr.GetEntries(); ++i) {
            TObject *obj = arr.At(i);
            T newobj = dynamic_cast<T>(obj);
            if (!newobj) {
                std::string classname = std::string(obj->IsA()->GetName());
                throw std::runtime_error("tcatorvec cannot convert object " + classname + " to target class ");
            }
            v.emplace_back(newobj);
        }
        return ROOT::VecOps::RVec<T>(v);
    }

    template ROOT::VecOps::RVec<ND::TGlobalReconModule::TGlobalPID *>
    tcatorvec<ND::TGlobalReconModule::TGlobalPID *>(const TClonesArray &arr);

    class T2KDataFrameWrapper {
    private:
        ROOT::RDataFrame rdf;
        ROOT::RDF::RInterface<ROOT::RDF::RDFDetail::RLoopManager> interface;

    public:
        T2KDataFrameWrapper(ROOT::RDataFrame &&rdf,
                            ROOT::RDF::RInterface<ROOT::RDF::RDFDetail::RLoopManager> &&interface)
                : rdf(rdf), interface(interface) {
        }

        ROOT::RDF::RInterface<ROOT::RDF::RDFDetail::RLoopManager> *operator->() { return &interface; }

    };

    T2KDataFrameWrapper MakeT2KDataFrame(std::vector<std::string> fileNameGlob,
                                                bool enabletruth = true,
                                                std::vector<std::string> treeNames = {"HeaderDir/BasicHeader",
                                                                                      "HeaderDir/BasicDataQuality",
                                                                                      "HeaderDir/BeamSummaryData",
                                                                                      "ReconDir/Global",
                                                                                      "TruthDir/GRooTrackerVtx",
                                                                                      "TruthDir/Trajectories",
                                                                                      "TruthDir/Vertices",
                                                },
                                                BunchTiming bunchTiming = {});
}
#endif
