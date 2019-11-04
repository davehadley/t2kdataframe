//#define T2K_USE_GENIE

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
#include "ND__NRooTrackerVtx.h"

#include <memory>
#include <map>

namespace T2K {

    class RooVtxWrapper {
    private:
        ND::RooTrackerVtxBase* ptr;
        ND::GRooTrackerVtx* genie;
        ND::NRooTrackerVtx* neut;
    public:
    RooVtxWrapper(ND::RooTrackerVtxBase* vtx=0);
    RooVtxWrapper(const RooVtxWrapper& cpy) = default;
    ND::RooTrackerVtxBase* getptr();
    TObjString* EvtCode();
   int         EvtNum() const;
   double      EvtXSec() const;
   double      EvtDXSec() const;
   double      EvtWght() const;
   double      EvtProb() const;
   double*      EvtVtx();
   int         StdHepN() const;
   Int_t      *StdHepPdg();
   Int_t      *StdHepStatus();
   double      (*StdHepX4())[4];
   double      (*StdHepP4())[4];
   double      (*StdHepPolz())[3];
   Int_t      *StdHepFd();
   Int_t      *StdHepLd();
   Int_t      *StdHepFm();
   Int_t      *StdHepLm();
   int         G2NeutEvtCode() const;
   TObjString* GeomPath();
   TObjString* GeneratorName();
   TObjString* OrigFileName();
   TObjString* OrigTreeName();
   int         OrigEvtNum() const;
   int         OrigTreeEntries() const;
   double      OrigTreePOT() const;
   double      TimeInSpill() const;
   int         TruthVertexID() const;

   int NeutrinoPdg() const;
   double NeutrinoEnergy() const;

   static ROOT::VecOps::RVec<RooVtxWrapper> wrapclonesarray(const TClonesArray& arr);
    };

    //Make some type aliases
    using Reco = ND::TGlobalReconModule::TGlobalPID;
    using RecoVec = ROOT::VecOps::RVec<Reco *>;
    using TruthTraj = ND::TTruthTrajectoriesModule::TTruthTrajectory;
    using TruthTrajVec = ROOT::VecOps::RVec<TruthTraj *>;
    using TruthVtx = ND::TTruthVerticesModule::TTruthVertex;
    using TruthVtxVec = ROOT::VecOps::RVec<TruthVtx *>;
    using TruthRooVtx = RooVtxWrapper;

    using TruthRooVtxVec = ROOT::VecOps::RVec<TruthRooVtx *>;

    class T2KDataSource final : public ROOT::RDF::RDataSource {
    public:

        class T2KSingleRootChainDS final : public ROOT::RDF::RDataSource {
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
            T2KSingleRootChainDS(std::string_view treeName, std::vector<std::string> fileNameGlob);
            ~T2KSingleRootChainDS();
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

        int _nslots;
        mutable std::vector<std::string> columnnames;
        std::vector<std::string> treenames;
        ULong64_t processed;
        int nbunches;
        std::vector<unique_ptr<T2K::T2KDataSource::T2KSingleRootChainDS> > _chainreaders;
        mutable std::map<std::string, int> _columnnamemap;

    public:

        ULong64_t GetEntries();
        int GetNSlots();
        TChain &GetSlotChain(int slot, int treeindex);
        std::vector<void *> GetColumnReadersImpl(std::string_view name, const std::type_info &t);
        T2KDataSource(std::vector<std::string> treenames_, std::vector<std::string> filename, int nbunches_);
        std::string GetTypeName(std::string_view colName) const;
        virtual const std::vector<std::string> &GetColumnNames() const;
        std::string getBranchName(std::string colName) const;
        std::string getBranchName(std::string_view colName) const;
        bool HasColumn(std::string_view colName) const;
        void InitSlot(unsigned int slot, ULong64_t firstEntry);
        virtual void Initialise();
        std::vector<std::pair<ULong64_t, ULong64_t> > GetEntryRanges();
        bool SetEntry(unsigned int slot, ULong64_t entry);
        void SetNSlots(unsigned int nSlots);
        void Finalise();
        void FinaliseSlot(unsigned int i);

    private:
        ULong64_t GetEntryNumber(ULong64_t entry);
        ULong64_t GetBunchNumber(ULong64_t entry);

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
        BunchTiming(ROOT::VecOps::RVec<double> centres, double window);
        BunchTiming(ROOT::VecOps::RVec<double> lowedges, ROOT::VecOps::RVec<double> highedges);
        BunchTiming(double window = 150.0);
        int size() const;
        double low(int bunch, int run) const;
        double high(int bunch, int run) const;
        double center(int bunch, int run) const;
        ROOT::VecOps::RVec<double> allbunchcentres(int run) const;
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
                                                                               "TruthDir/NRooTrackerVtx",
                                                                               "TruthDir/Trajectories",
                                                                               "TruthDir/Vertices"
                                         },
                                         BunchTiming bunchTiming = {}
                                         );
}
#endif
