////////////////////////////////////////////////////////////////////////
// Class:       TrajectoryMCSNtuple
// Plugin Type: analyzer (art v2_05_00)
// File:        TrajectoryMCSNtuple_module.cc
//
// Generated at Mon Feb  6 10:06:04 2017 by Giuseppe Cerati using cetskelgen
// from cetlib version v1_21_00.
//
// Updated 2020 A. Mastbaum - MCC9, NuCC filter, truth matching
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "larcore/Geometry/Geometry.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/MCSFitResult.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardata/RecoObjects/TrackingPlaneHelper.h"
#include "lardata/RecoObjects/TrackState.h"

#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Hit.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

class TrajectoryMCSNtuple;

class TrajectoryMCSNtuple : public art::EDAnalyzer {
public:

  struct Inputs {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<std::string> inputLabel {
      Name("inputLabel"),
      Comment("Label of recob::TrackTrajectory Collection to be fit")
     };
    fhicl::Atom<std::string> MuMCSLabel {
      Name("MuMCSLabel"),
      Comment("Label of MCS fit collection with Mu hypothesis")
     };
    fhicl::Atom<std::string> PMCSLabel {
      Name("PMCSLabel"),
      Comment("Label of MCS fit collection with P hypothesis")
     };
  };

  struct Config {
    using Name = fhicl::Name;
    fhicl::Table<TrajectoryMCSNtuple::Inputs> inputs {
      Name("inputs"),
    };
    fhicl::Table<trkf::TrajectoryMCSFitter::Config> mcsfitter {
      Name("mcsfitter")
    };
    fhicl::Table<trkf::TrajectoryMCSFitter::Config> mcsfittermc {
      Name("mcsfittermc")
    };
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit TrajectoryMCSNtuple(Parameters const & p);
  ~TrajectoryMCSNtuple() {}

  TrajectoryMCSNtuple(TrajectoryMCSNtuple const &) = delete;
  TrajectoryMCSNtuple(TrajectoryMCSNtuple &&) = delete;
  TrajectoryMCSNtuple & operator = (TrajectoryMCSNtuple const &) = delete;
  TrajectoryMCSNtuple & operator = (TrajectoryMCSNtuple &&) = delete;

  void analyze(art::Event const & e) override;

  void beginJob() override;

  void resetTree();

  bool isPfpFromNeutrino(art::ValidHandle<std::vector<recob::PFParticle> > pfph, size_t ipfp);

  art::Ptr<simb::MCParticle> getAssocMCParticle(
      art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &hittruth,
      const std::vector<art::Ptr<recob::Hit>> &hits,
      float &purity);

private:
  std::string inputTracksLabel;
  std::string MuMCSLabel;
  std::string PMCSLabel;
  trkf::TrackMomentumCalculator tmc;
  trkf::TrajectoryMCSFitter mcsfitter;
  trkf::TrajectoryMCSFitter mcsfittermc;
  TTree* tree;

  int    run, subrun, eventid;
  float trkLength;
  float trkMom_MuFwd, trkMomErr_MuFwd, trkLL_MuFwd;
  float trkMom_MuBwd, trkMomErr_MuBwd, trkLL_MuBwd;
  float trkMom_Mu, trkMomErr_Mu, trkLL_Mu;
  float trkDeltaLL_Mu;
  int    trkIsBestFwd_Mu;
  float trkMom_PFwd, trkMomErr_PFwd, trkLL_PFwd;
  float trkMom_PBwd, trkMomErr_PBwd, trkLL_PBwd;
  float trkMom_P, trkMomErr_P, trkLL_P;
  float trkDeltaLL_P;
  int    trkIsBestFwd_P;
  std::vector<float> trkSegRadLengths;
  std::vector<float> trkScattAngles;

  int    trkNHits, trkIsContained;
  float trkMom_RangeMu, trkMom_RangeP;
  int    trkID_MCSRange;
  float trkStartPosX, trkStartPosY, trkStartPosZ;
  float trkStartDirX, trkStartDirY, trkStartDirZ;
  float trkEndPosX, trkEndPosY, trkEndPosZ;
  float trkEndDirX, trkEndDirY, trkEndDirZ;
  std::vector<int> trkSegNHits;
  std::vector<float> trkSegStartPosX, trkSegStartPosY, trkSegStartPosZ;
  std::vector<float> trkSegEndPosX, trkSegEndPosY, trkSegEndPosZ;
  std::vector<float> trkSegDirX, trkSegDirY, trkSegDirZ;

  float trkLength_rm10;
  float trkMom_MuFwd_rm10, trkMomErr_MuFwd_rm10, trkLL_MuFwd_rm10;
  float trkMom_MuBwd_rm10, trkMomErr_MuBwd_rm10, trkLL_MuBwd_rm10;
  float trkMom_Mu_rm10, trkMomErr_Mu_rm10, trkLL_Mu_rm10;
  float trkDeltaLL_Mu_rm10;
  int    trkIsBestFwd_Mu_rm10;

  float trkLength_rm50;
  float trkMom_MuFwd_rm50, trkMomErr_MuFwd_rm50, trkLL_MuFwd_rm50;
  float trkMom_MuBwd_rm50, trkMomErr_MuBwd_rm50, trkLL_MuBwd_rm50;
  float trkMom_Mu_rm50, trkMomErr_Mu_rm50, trkLL_Mu_rm50;
  float trkDeltaLL_Mu_rm50;
  int    trkIsBestFwd_Mu_rm50;

  float trkLength_rm80;
  float trkMom_MuFwd_rm80, trkMomErr_MuFwd_rm80, trkLL_MuFwd_rm80;
  float trkMom_MuBwd_rm80, trkMomErr_MuBwd_rm80, trkLL_MuBwd_rm80;
  float trkMom_Mu_rm80, trkMomErr_Mu_rm80, trkLL_Mu_rm80;
  float trkDeltaLL_Mu_rm80;
  int    trkIsBestFwd_Mu_rm80;

  int    trkid0;
  int    trkid1;
  int    trkid2;
  float trkpida0;
  float trkpida1;
  float trkpida2;

  float simMom, simMomEnd, simMatchPurity;
  float simLength;
  float simStartPosX, simStartPosY, simStartPosZ;
  float simStartDirX, simStartDirY, simStartDirZ;
  float simEndPosX, simEndPosY, simEndPosZ;
  float simEndDirX, simEndDirY, simEndDirZ;
  std::vector<float> simMcStepPosX, simMcStepPosY, simMcStepPosZ;
  int    simID;
  std::string simProc;
  int    simIsContained;
  int    simAndTrkSameDir;

  float simMom_MuFwd, simMomErr_MuFwd, simLL_MuFwd;
  float simMom_MuBwd, simMomErr_MuBwd, simLL_MuBwd;
  float simMom_Mu, simMomErr_Mu, simLL_Mu;
  float simDeltaLL_Mu;
  int    simIsBestFwd_Mu;

  std::vector<float> simSegRadLengths;
  std::vector<float> simScattAngles;
  std::vector<float> simL;
  std::vector<float> simP;
  std::vector<float> simBeta;
  std::vector<float> simE;
  std::vector<float> simRP;
  std::vector<float> simRBeta;
  std::vector<float> simRE;
  std::vector<float> simRMS;
  std::vector<float> simRMSHL;
  std::vector<float> simRMSR;
  std::vector<float> simA3d;
  std::vector<float> simA2dxz;
  std::vector<float> simA2dyz;
  std::vector<float> simAAvg3d;
  std::vector<float> simAAvg2dxz;
  std::vector<float> simAAvg2dyz;

  int passNuCC;
  float nuScore;
  float fmScore;
};


void TrajectoryMCSNtuple::resetTree() {
  run = -999;
  subrun = -999;
  eventid = -999;
  trkLength = -999;
  trkMom_MuFwd = -999;
  trkMomErr_MuFwd = -999;
  trkLL_MuFwd = -999;
  trkMom_MuBwd = -999;
  trkMomErr_MuBwd = -999;
  trkLL_MuBwd = -999;
  trkMom_Mu = -999;
  trkMomErr_Mu = -999;
  trkLL_Mu = -999;
  trkDeltaLL_Mu = -999;
  trkIsBestFwd_Mu = -999;
  trkMom_PFwd = -999;
  trkMomErr_PFwd = -999;
  trkLL_PFwd = -999;
  trkMom_PBwd = -999;
  trkMomErr_PBwd = -999;
  trkLL_PBwd = -999;
  trkMom_P = -999;
  trkMomErr_P = -999;
  trkLL_P = -999;
  trkDeltaLL_P = -999;
  trkIsBestFwd_P = -999;
  trkSegRadLengths.clear();
  trkScattAngles.clear();

  trkNHits = -999;
  trkIsContained = -999;
  trkMom_RangeMu = -999;
  trkMom_RangeP = -999;
  trkID_MCSRange = -999;
  trkStartPosX = -999;
  trkStartPosY = -999;
  trkStartPosZ = -999;
  trkStartDirX = -999;
  trkStartDirY = -999;
  trkStartDirZ = -999;
  trkEndPosX = -999;
  trkEndPosY = -999;
  trkEndPosZ = -999;
  trkEndDirX = -999;
  trkEndDirY = -999;
  trkEndDirZ = -999;
  trkSegNHits.clear();
  trkSegStartPosX.clear();
  trkSegStartPosY.clear();
  trkSegStartPosZ.clear();
  trkSegEndPosX.clear();
  trkSegEndPosY.clear();
  trkSegEndPosZ.clear();
  trkSegDirX.clear();
  trkSegDirY.clear();
  trkSegDirZ.clear();

  trkLength_rm10 = -999;
  trkMom_MuFwd_rm10 = -999;
  trkMomErr_MuFwd_rm10 = -999;
  trkLL_MuFwd_rm10 = -999;
  trkMom_MuBwd_rm10 = -999;
  trkMomErr_MuBwd_rm10 = -999;
  trkLL_MuBwd_rm10 = -999;
  trkMom_Mu_rm10 = -999;
  trkMomErr_Mu_rm10 = -999;
  trkLL_Mu_rm10 = -999;
  trkDeltaLL_Mu_rm10 = -999;
  trkIsBestFwd_Mu_rm10 = -999;

  trkLength_rm50 = -999;
  trkMom_MuFwd_rm50 = -999;
  trkMomErr_MuFwd_rm50 = -999;
  trkLL_MuFwd_rm50 = -999;
  trkMom_MuBwd_rm50 = -999;
  trkMomErr_MuBwd_rm50 = -999;
  trkLL_MuBwd_rm50 = -999;
  trkMom_Mu_rm50 = -999;
  trkMomErr_Mu_rm50 = -999;
  trkLL_Mu_rm50 = -999;
  trkDeltaLL_Mu_rm50 = -999;
  trkIsBestFwd_Mu_rm50 = -999;

  trkLength_rm80 = -999;
  trkMom_MuFwd_rm80 = -999;
  trkMomErr_MuFwd_rm80 = -999;
  trkLL_MuFwd_rm80 = -999;
  trkMom_MuBwd_rm80 = -999;
  trkMomErr_MuBwd_rm80 = -999;
  trkLL_MuBwd_rm80 = -999;
  trkMom_Mu_rm80 = -999;
  trkMomErr_Mu_rm80 = -999;
  trkLL_Mu_rm80 = -999;
  trkDeltaLL_Mu_rm80 = -999;
  trkIsBestFwd_Mu_rm80 = -999;

  trkid0 = -999;
  trkid1 = -999;
  trkid2 = -999;
  trkpida0 = -999;
  trkpida1 = -999;
  trkpida2 = -999;

  simMom = -999;
  simMatchPurity = -999;
  simMomEnd = -999;
  simLength = -999;
  simStartPosX = -999;
  simStartPosY = -999;
  simStartPosZ = -999;
  simStartDirX = -999;
  simStartDirY = -999;
  simStartDirZ = -999;
  simEndPosX = -999;
  simEndPosY = -999;
  simEndPosZ = -999;
  simEndDirX = -999;
  simEndDirY = -999;
  simEndDirZ = -999;
  simMcStepPosX.clear();
  simMcStepPosY.clear();
  simMcStepPosZ.clear();
  simID = -999;
  simProc = "";
  simIsContained = -999;
  simAndTrkSameDir = -999;
  simMom_MuFwd = -999;
  simMomErr_MuFwd = -999;
  simLL_MuFwd = -999;
  simMom_MuBwd = -999;
  simMomErr_MuBwd = -999;
  simLL_MuBwd = -999;
  simMom_Mu = -999;
  simMomErr_Mu = -999;
  simLL_Mu = -999;
  simDeltaLL_Mu = -999;
  simIsBestFwd_Mu = -999;
  simSegRadLengths.clear();
  simScattAngles.clear();
  simL.clear();
  simE.clear();
  simP.clear();
  simBeta.clear();
  simRE.clear();
  simRP.clear();
  simRBeta.clear();
  simRMS.clear();
  simRMSHL.clear();
  simRMSR.clear();
  simA3d.clear();
  simA2dxz.clear();
  simA2dyz.clear();
  simAAvg3d.clear();
  simAAvg2dxz.clear();
  simAAvg2dyz.clear();

  passNuCC = -999;
  nuScore = -999;
  fmScore = -999;
}


void TrajectoryMCSNtuple::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;

  tree = tfs->make<TTree>("tree", "tree");

  tree->Branch("run", &run,"run/I");
  tree->Branch("subrun", &subrun, "subrun/I");
  tree->Branch("eventid", &eventid, "eventid/I");

  tree->Branch("trkLength", &trkLength, "trkLength/F");
  tree->Branch("trkMom_MuFwd"   , &trkMom_MuFwd   , "trkMom_MuFwd/F"   );
  tree->Branch("trkMomErr_MuFwd", &trkMomErr_MuFwd, "trkMomErr_MuFwd/F");
  tree->Branch("trkLL_MuFwd"    , &trkLL_MuFwd    , "trkLL_MuFwd/F"    );
  tree->Branch("trkMom_MuBwd"   , &trkMom_MuBwd   , "trkMom_MuBwd/F"   );
  tree->Branch("trkMomErr_MuBwd", &trkMomErr_MuBwd, "trkMomErr_MuBwd/F");
  tree->Branch("trkLL_MuBwd"    , &trkLL_MuBwd    , "trkLL_MuBwd/F"    );
  tree->Branch("trkMom_Mu"      , &trkMom_Mu      , "trkMom_Mu/F"      );
  tree->Branch("trkMomErr_Mu"   , &trkMomErr_Mu   , "trkMomErr_Mu/F"   );
  tree->Branch("trkLL_Mu"       , &trkLL_Mu       , "trkLL_Mu/F"       );
  tree->Branch("trkDeltaLL_Mu"  , &trkDeltaLL_Mu  , "trkDeltaLL_Mu/F"  );
  tree->Branch("trkIsBestFwd_Mu", &trkIsBestFwd_Mu, "trkIsBestFwd_Mu/I");

  tree->Branch("trkMom_PFwd"   , &trkMom_PFwd   , "trkMom_PFwd/F"   );
  tree->Branch("trkMomErr_PFwd", &trkMomErr_PFwd, "trkMomErr_PFwd/F");
  tree->Branch("trkLL_PFwd"    , &trkLL_PFwd    , "trkLL_PFwd/F"    );
  tree->Branch("trkMom_PBwd"   , &trkMom_PBwd   , "trkMom_PBwd/F"   );
  tree->Branch("trkMomErr_PBwd", &trkMomErr_PBwd, "trkMomErr_PBwd/F");
  tree->Branch("trkLL_PBwd"    , &trkLL_PBwd    , "trkLL_PBwd/F"    );
  tree->Branch("trkMom_P"      , &trkMom_P      , "trkMom_P/F"      );
  tree->Branch("trkMomErr_P"   , &trkMomErr_P   , "trkMomErr_P/F"   );
  tree->Branch("trkLL_P"       , &trkLL_P       , "trkLL_P/F"       );
  tree->Branch("trkDeltaLL_P"  , &trkDeltaLL_P  , "trkDeltaLL_P/F"  );
  tree->Branch("trkIsBestFwd_P", &trkIsBestFwd_P, "trkIsBestFwd_P/I");

  tree->Branch("trkSegRadLengths", &trkSegRadLengths);
  tree->Branch("trkScattAngles"  , &trkScattAngles  );

  tree->Branch("trkNHits", &trkNHits, "trkNHits/I");
  tree->Branch("trkIsContained", &trkIsContained, "trkIsContained/I");
  tree->Branch("trkMom_RangeMu", &trkMom_RangeMu, "trkMom_RangeMu/F");
  tree->Branch("trkMom_RangeP" , &trkMom_RangeP , "trkMom_RangeP/F" );
  tree->Branch("trkID_MCSRange", &trkID_MCSRange, "trkID_MCSRange/I");
  tree->Branch("trkStartPosX", &trkStartPosX, "trkStartPosX/F");
  tree->Branch("trkStartPosY", &trkStartPosY, "trkStartPosY/F");
  tree->Branch("trkStartPosZ", &trkStartPosZ, "trkStartPosZ/F");
  tree->Branch("trkStartDirX", &trkStartDirX, "trkStartDirX/F");
  tree->Branch("trkStartDirY", &trkStartDirY, "trkStartDirY/F");
  tree->Branch("trkStartDirZ", &trkStartDirZ, "trkStartDirZ/F");
  tree->Branch("trkEndPosX"  , &trkEndPosX  , "trkEndPosX/F"  );
  tree->Branch("trkEndPosY"  , &trkEndPosY  , "trkEndPosY/F"  );
  tree->Branch("trkEndPosZ"  , &trkEndPosZ  , "trkEndPosZ/F"  );
  tree->Branch("trkEndDirX"  , &trkEndDirX  , "trkEndDirX/F"  );
  tree->Branch("trkEndDirY"  , &trkEndDirY  , "trkEndDirY/F"  );
  tree->Branch("trkEndDirZ"  , &trkEndDirZ  , "trkEndDirZ/F"  );

  tree->Branch("trkSegNHits"    , &trkSegNHits    );
  tree->Branch("trkSegStartPosX", &trkSegStartPosX);
  tree->Branch("trkSegStartPosY", &trkSegStartPosY);
  tree->Branch("trkSegStartPosZ", &trkSegStartPosZ);
  tree->Branch("trkSegEndPosX"  , &trkSegEndPosX  );
  tree->Branch("trkSegEndPosY"  , &trkSegEndPosY  );
  tree->Branch("trkSegEndPosZ"  , &trkSegEndPosZ  );
  tree->Branch("trkSegDirX"     , &trkSegDirX          );
  tree->Branch("trkSegDirY"     , &trkSegDirY          );
  tree->Branch("trkSegDirZ"     , &trkSegDirZ     );

  tree->Branch("trkLength_rm10", &trkLength_rm10, "trkLength_rm10/F");
  tree->Branch("trkMom_MuFwd_rm10"   , &trkMom_MuFwd_rm10   , "trkMom_MuFwd_rm10/F"   );
  tree->Branch("trkMomErr_MuFwd_rm10", &trkMomErr_MuFwd_rm10, "trkMomErr_MuFwd_rm10/F");
  tree->Branch("trkLL_MuFwd_rm10"    , &trkLL_MuFwd_rm10    , "trkLL_MuFwd_rm10/F"    );
  tree->Branch("trkMom_MuBwd_rm10"   , &trkMom_MuBwd_rm10   , "trkMom_MuBwd_rm10/F"   );
  tree->Branch("trkMomErr_MuBwd_rm10", &trkMomErr_MuBwd_rm10, "trkMomErr_MuBwd_rm10/F");
  tree->Branch("trkLL_MuBwd_rm10"    , &trkLL_MuBwd_rm10    , "trkLL_MuBwd_rm10/F"    );
  tree->Branch("trkMom_Mu_rm10"      , &trkMom_Mu_rm10      , "trkMom_Mu_rm10/F"      );
  tree->Branch("trkMomErr_Mu_rm10"   , &trkMomErr_Mu_rm10   , "trkMomErr_Mu_rm10/F"   );
  tree->Branch("trkLL_Mu_rm10"       , &trkLL_Mu_rm10       , "trkLL_Mu_rm10/F"       );
  tree->Branch("trkDeltaLL_Mu_rm10"  , &trkDeltaLL_Mu_rm10  , "trkDeltaLL_Mu_rm10/F"  );
  tree->Branch("trkIsBestFwd_Mu_rm10", &trkIsBestFwd_Mu_rm10, "trkIsBestFwd_Mu_rm10/I");

  tree->Branch("trkLength_rm50", &trkLength_rm50, "trkLength_rm50/F");
  tree->Branch("trkMom_MuFwd_rm50"   , &trkMom_MuFwd_rm50   , "trkMom_MuFwd_rm50/F"   );
  tree->Branch("trkMomErr_MuFwd_rm50", &trkMomErr_MuFwd_rm50, "trkMomErr_MuFwd_rm50/F");
  tree->Branch("trkLL_MuFwd_rm50"    , &trkLL_MuFwd_rm50    , "trkLL_MuFwd_rm50/F"    );
  tree->Branch("trkMom_MuBwd_rm50"   , &trkMom_MuBwd_rm50   , "trkMom_MuBwd_rm50/F"   );
  tree->Branch("trkMomErr_MuBwd_rm50", &trkMomErr_MuBwd_rm50, "trkMomErr_MuBwd_rm50/F");
  tree->Branch("trkLL_MuBwd_rm50"    , &trkLL_MuBwd_rm50    , "trkLL_MuBwd_rm50/F"    );
  tree->Branch("trkMom_Mu_rm50"      , &trkMom_Mu_rm50      , "trkMom_Mu_rm50/F"      );
  tree->Branch("trkMomErr_Mu_rm50"   , &trkMomErr_Mu_rm50   , "trkMomErr_Mu_rm50/F"   );
  tree->Branch("trkLL_Mu_rm50"       , &trkLL_Mu_rm50       , "trkLL_Mu_rm50/F"       );
  tree->Branch("trkDeltaLL_Mu_rm50"  , &trkDeltaLL_Mu_rm50  , "trkDeltaLL_Mu_rm50/F"  );
  tree->Branch("trkIsBestFwd_Mu_rm50", &trkIsBestFwd_Mu_rm50, "trkIsBestFwd_Mu_rm50/I");

  tree->Branch("trkLength_rm80", &trkLength_rm80, "trkLength_rm80/F");
  tree->Branch("trkMom_MuFwd_rm80"   , &trkMom_MuFwd_rm80   , "trkMom_MuFwd_rm80/F"   );
  tree->Branch("trkMomErr_MuFwd_rm80", &trkMomErr_MuFwd_rm80, "trkMomErr_MuFwd_rm80/F");
  tree->Branch("trkLL_MuFwd_rm80"    , &trkLL_MuFwd_rm80    , "trkLL_MuFwd_rm80/F"    );
  tree->Branch("trkMom_MuBwd_rm80"   , &trkMom_MuBwd_rm80   , "trkMom_MuBwd_rm80/F"   );
  tree->Branch("trkMomErr_MuBwd_rm80", &trkMomErr_MuBwd_rm80, "trkMomErr_MuBwd_rm80/F");
  tree->Branch("trkLL_MuBwd_rm80"    , &trkLL_MuBwd_rm80    , "trkLL_MuBwd_rm80/F"    );
  tree->Branch("trkMom_Mu_rm80"      , &trkMom_Mu_rm80      , "trkMom_Mu_rm80/F"      );
  tree->Branch("trkMomErr_Mu_rm80"   , &trkMomErr_Mu_rm80   , "trkMomErr_Mu_rm80/F"   );
  tree->Branch("trkLL_Mu_rm80"       , &trkLL_Mu_rm80       , "trkLL_Mu_rm80/F"       );
  tree->Branch("trkDeltaLL_Mu_rm80"  , &trkDeltaLL_Mu_rm80  , "trkDeltaLL_Mu_rm80/F"  );
  tree->Branch("trkIsBestFwd_Mu_rm80", &trkIsBestFwd_Mu_rm80, "trkIsBestFwd_Mu_rm80/I");

  tree->Branch("trkid0"  , &trkid0  , "trkid0/I"  );
  tree->Branch("trkid1"  , &trkid1  , "trkid1/I"  );
  tree->Branch("trkid2"  , &trkid2  , "trkid2/I"  );
  tree->Branch("trkpida0", &trkpida0, "trkpida0/F");
  tree->Branch("trkpida1", &trkpida1, "trkpida1/F");
  tree->Branch("trkpida2", &trkpida2, "trkpida2/F");

  tree->Branch("simMom"   , &simMom   , "simMom/F"   );
  tree->Branch("simMatchPurity", &simMatchPurity, "simMatchPurity/F");
  tree->Branch("simMomEnd", &simMomEnd, "simMomEnd/F"   );
  tree->Branch("simLength", &simLength, "simLength/F");
  tree->Branch("simStartPosX", &simStartPosX, "simStartPosX/F");
  tree->Branch("simStartPosY", &simStartPosY, "simStartPosY/F");
  tree->Branch("simStartPosZ", &simStartPosZ, "simStartPosZ/F");
  tree->Branch("simStartDirX", &simStartDirX, "simStartDirX/F");
  tree->Branch("simStartDirY", &simStartDirY, "simStartDirY/F");
  tree->Branch("simStartDirZ", &simStartDirZ, "simStartDirZ/F");
  tree->Branch("simEndPosX"  , &simEndPosX  , "simEndPosX/F"  );
  tree->Branch("simEndPosY"  , &simEndPosY  , "simEndPosY/F"  );
  tree->Branch("simEndPosZ"  , &simEndPosZ  , "simEndPosZ/F"  );
  tree->Branch("simEndDirX"  , &simEndDirX  , "simEndDirX/F"  );
  tree->Branch("simEndDirY"  , &simEndDirY  , "simEndDirY/F"  );
  tree->Branch("simEndDirZ"  , &simEndDirZ  , "simEndDirZ/F"  );
  tree->Branch("simMcStepPosX"  , &simMcStepPosX  );
  tree->Branch("simMcStepPosY"  , &simMcStepPosY  );
  tree->Branch("simMcStepPosZ"  , &simMcStepPosZ  );
  tree->Branch("simID", &simID, "simID/I");
  tree->Branch("simProc", &simProc);
  tree->Branch("simIsContained", &simIsContained, "simIsContained/I");
  tree->Branch("simAndTrkSameDir", &simAndTrkSameDir, "simAndTrkSameDir/I");

  tree->Branch("simMom_MuFwd"   , &simMom_MuFwd   , "simMom_MuFwd/F"   );
  tree->Branch("simMomErr_MuFwd", &simMomErr_MuFwd, "simMomErr_MuFwd/F");
  tree->Branch("simLL_MuFwd"    , &simLL_MuFwd    , "simLL_MuFwd/F"    );
  tree->Branch("simMom_MuBwd"   , &simMom_MuBwd   , "simMom_MuBwd/F"   );
  tree->Branch("simMomErr_MuBwd", &simMomErr_MuBwd, "simMomErr_MuBwd/F");
  tree->Branch("simLL_MuBwd"    , &simLL_MuBwd    , "simLL_MuBwd/F"    );
  tree->Branch("simMom_Mu"      , &simMom_Mu      , "simMom_Mu/F"      );
  tree->Branch("simMomErr_Mu"   , &simMomErr_Mu   , "simMomErr_Mu/F"   );
  tree->Branch("simLL_Mu"       , &simLL_Mu       , "simLL_Mu/F"       );
  tree->Branch("simDeltaLL_Mu"  , &simDeltaLL_Mu  , "simDeltaLL_Mu/F"  );
  tree->Branch("simIsBestFwd_Mu", &simIsBestFwd_Mu, "simIsBestFwd_Mu/I");

  tree->Branch("simSegRadLengths", &simSegRadLengths);
  tree->Branch("simScattAngles"  , &simScattAngles  );
  tree->Branch("simL"  , &simL  );
  tree->Branch("simE"  , &simE  );
  tree->Branch("simP"  , &simP  );
  tree->Branch("simBeta"  , &simBeta  );
  tree->Branch("simRE"  , &simRE  );
  tree->Branch("simRP"  , &simRP  );
  tree->Branch("simRBeta"  , &simRBeta  );
  tree->Branch("simRMS"  , &simRMS  );
  tree->Branch("simRMSHL"  , &simRMSHL  );
  tree->Branch("simRMSR"  , &simRMSR  );
  tree->Branch("simA3d"  , &simA3d  );
  tree->Branch("simA2dxz"  , &simA2dxz  );
  tree->Branch("simA2dyz"  , &simA2dyz  );
  tree->Branch("simAAvg3d"  , &simAAvg3d  );
  tree->Branch("simAAvg2dxz"  , &simAAvg2dxz  );
  tree->Branch("simAAvg2dyz"  , &simAAvg2dyz  );

  tree->Branch("passNuCC", &passNuCC, "passNuCC/I");
  tree->Branch("nuScore", &nuScore, "nuScore/F");
  tree->Branch("fmScore", &fmScore, "fmScore/F");
}


TrajectoryMCSNtuple::TrajectoryMCSNtuple(Parameters const& p)
  : EDAnalyzer(p)
  , inputTracksLabel(p().inputs().inputLabel())
  , MuMCSLabel(p().inputs().MuMCSLabel())
  , PMCSLabel(p().inputs().PMCSLabel())
  , mcsfitter(p().mcsfitter)
  , mcsfittermc(p().mcsfittermc) {}


void TrajectoryMCSNtuple::analyze(art::Event const & e) {
  using namespace std;
  using namespace trkf;
  using namespace recob;
  using namespace recob::tracking;

  // Truth matching info
  art::ValidHandle<std::vector<simb::MCParticle> >* simParticles = nullptr;
  art::Handle<std::vector<simb::MCParticle> > test;
  art::InputTag SimParticleInputTag("largeant");

  art::ValidHandle<std::vector<recob::Hit> > hitHandle = e.getValidHandle<std::vector<recob::Hit> >("gaushit");
  art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>* simHits = nullptr;
  art::InputTag SimHitInputTag("gaushitTruthMatch");

  if (e.getByLabel(SimParticleInputTag, test)) {
    try {
      simParticles = \
        new art::ValidHandle<std::vector<simb::MCParticle> >(
          e.getValidHandle<std::vector<simb::MCParticle> >(SimParticleInputTag));
      simHits = \
        new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(
          hitHandle, e, SimHitInputTag);
    } catch (...) {
      std::cout << "MC event with invalid simb::MCParticle collection... skipping" << std::endl;
      return;
    }
  }

  // Pandora info
  art::InputTag PFInputTag("pandora");
  art::InputTag TrackInputTag(inputTracksLabel);

  art::ValidHandle<std::vector<recob::PFParticle> > Pfps = e.getValidHandle<std::vector<recob::PFParticle> >(PFInputTag);
  auto const& assocTracks = *e.getValidHandle<art::Assns<recob::PFParticle, recob::Track> >(TrackInputTag);
  art::ValidHandle<std::vector<recob::Track> > Tracks = e.getValidHandle<std::vector<recob::Track> >(TrackInputTag);

  art::FindManyP<larpandoraobj::PFParticleMetadata> pfPartToMetadataAssoc(Pfps, e, PFInputTag);
  art::FindManyP<recob::Hit> assocHits(Tracks, e, PFInputTag);

  // NuCC filter info
  art::InputTag NuCCInputTag("NuCCproducer");
  art::Handle<art::Assns<recob::PFParticle,anab::T0,void> > nucc_t0a;
  e.getByLabel(NuCCInputTag, nucc_t0a);

  art::InputTag fmtag("flashmatch");
  art::FindManyP<anab::T0> pfp_t0_assn_v(Pfps, e, fmtag);

  // MCS info
  art::InputTag MuMCSInputTag(MuMCSLabel);
  art::ValidHandle<std::vector<recob::MCSFitResult> > MCSMu = e.getValidHandle<std::vector<recob::MCSFitResult> >(MuMCSInputTag);

  art::InputTag PMCSInputTag(PMCSLabel);
  art::ValidHandle<std::vector<recob::MCSFitResult> > MCSP = e.getValidHandle<std::vector<recob::MCSFitResult> >(PMCSInputTag);

  assert(Tracks->size() == MCSMu->size());

  // Process NuCC filtered events only
  if (!nucc_t0a.isValid()) {
    return;
  }

  // Loop over reconstructed tracks
  for (size_t iTrack=0; iTrack<Tracks->size(); iTrack++) {
    const recob::Track& track = Tracks->at(iTrack);
    const recob::MCSFitResult& mcsMu = MCSMu->at(iTrack);
    const recob::MCSFitResult& mcsP = MCSP->at(iTrack);
    std::vector<art::Ptr<recob::Hit> > hits = assocHits.at(iTrack);

    // Get PFP for this track
    size_t assocpfpkey = 99999;
    for (auto at : assocTracks) {
      if (at.second.key()==iTrack) {
        assocpfpkey = at.first.key();
        break;
      }
    }

    // Track must be associated with a neutrino PFP
    if (!isPfpFromNeutrino(Pfps, assocpfpkey)) continue;

    // Track must be from the NuCC neutrino candidate PFP
    bool nucc_nu = false;
    for (auto ap : *nucc_t0a) {
      if (ap.first.key() == assocpfpkey) {
        nucc_nu = true;
        break;
      }
    }
    
    if (!nucc_nu) continue;

    // Extract neutrino scores
    float _fmscore = -999;
    float _topo_score = -999;
    bool nuScoreSet = false;

    for (size_t ip=0; ip<Pfps->size(); ip++) {
      auto const& pfp = Pfps->at(ip);

      if (!pfp.IsPrimary() ||
          !(std::abs(pfp.PdgCode()) == 12 ||
            std::abs(pfp.PdgCode()) == 14 ||
            std::abs(pfp.PdgCode()) == 16)) {
        continue;
      }

      auto const& pfParticleMetadataList(pfPartToMetadataAssoc.at(ip));

      if (!pfParticleMetadataList.empty()) {
        for (unsigned int j=0; j<pfParticleMetadataList.size(); j++) {
          const art::Ptr<larpandoraobj::PFParticleMetadata>& pfParticleMetadata(pfParticleMetadataList.at(j));
          const larpandoraobj::PFParticleMetadata::PropertiesMap& pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
          if (pfParticlePropertiesMap.find("NuScore") != pfParticlePropertiesMap.end()) {
            assert(!nuScoreSet);
            _topo_score = pfParticlePropertiesMap.at("NuScore");
            nuScoreSet = true;
          }
        }
      }

      _fmscore = pfp_t0_assn_v.at(ip).at(0)->TriggerConfidence();
    }

    // Fill output tree
    resetTree();
    nuScore = _topo_score;
    fmScore = _fmscore;
    passNuCC = (int) nucc_t0a.isValid();

    run = e.run();
    subrun = e.subRun();
    eventid = e.event();

    trkLength = track.Length();

    trkMom_MuFwd    = mcsMu.fwdMomentum();
    trkMomErr_MuFwd = mcsMu.fwdMomUncertainty();
    trkLL_MuFwd     = mcsMu.fwdLogLikelihood();
    trkMom_MuBwd    = mcsMu.bwdMomentum();
    trkMomErr_MuBwd = mcsMu.bwdMomUncertainty();
    trkLL_MuBwd     = mcsMu.bwdLogLikelihood();
    trkMom_Mu       = mcsMu.bestMomentum();
    trkMomErr_Mu    = mcsMu.bestMomUncertainty();
    trkLL_Mu        = mcsMu.bestLogLikelihood();
    trkDeltaLL_Mu   = mcsMu.deltaLogLikelihood();
    trkIsBestFwd_Mu = mcsMu.isBestFwd();

    trkMom_PFwd    = mcsP.fwdMomentum();
    trkMomErr_PFwd = mcsP.fwdMomUncertainty();
    trkLL_PFwd     = mcsP.fwdLogLikelihood();
    trkMom_PBwd    = mcsP.bwdMomentum();
    trkMomErr_PBwd = mcsP.bwdMomUncertainty();
    trkLL_PBwd     = mcsP.bwdLogLikelihood();
    trkMom_P       = mcsP.bestMomentum();
    trkMomErr_P    = mcsP.bestMomUncertainty();
    trkLL_P        = mcsP.bestLogLikelihood();
    trkDeltaLL_P   = mcsP.deltaLogLikelihood();
    trkIsBestFwd_P = mcsP.isBestFwd();

    trkSegRadLengths = mcsMu.segmentRadLengths();
    trkScattAngles   = mcsMu.scatterAngles();

    trkNHits = track.NPoints();

    Point_t start(track.Start().X(),track.Start().Y(),track.Start().Z());
    Point_t end(track.End().X(),track.End().Y(),track.End().Z());
    trkIsContained = (start.X()>30.  && start.X()<230.  && end.X()>30.  && end.X()<230. &&
                      start.Y()>-85. && start.Y()<85.   && end.Y()>-85. && end.Y()<85.  &&
                      start.Z()>30.  && start.Z()<1010. && end.Z()>30.  && end.Z()<1010.);
    trkMom_RangeMu = tmc.GetTrackMomentum(trkLength,13  );
    trkMom_RangeP  = tmc.GetTrackMomentum(trkLength,2212);
    trkID_MCSRange = ( ( fabs(trkMom_RangeMu-trkMom_Mu)<fabs(trkMom_RangeP-trkMom_P) ) ? 13 : 2212);
    trkStartPosX = start.X();
    trkStartPosY = start.Y();
    trkStartPosZ = start.Z();
    trkStartDirX = track.StartDirection().X();
    trkStartDirY = track.StartDirection().Y();
    trkStartDirZ = track.StartDirection().Z();
    trkEndPosX = end.X();
    trkEndPosY = end.Y();
    trkEndPosZ = end.Z();
    trkEndDirX = track.EndDirection().X();
    trkEndDirY = track.EndDirection().Y();
    trkEndDirZ = track.EndDirection().Z();
    vector<size_t> breakpoints;
    vector<float> segradlengths;
    vector<float> cumseglens;
    mcsfitter.breakTrajInSegments(track.Trajectory(), breakpoints, segradlengths, cumseglens);

    std::vector<int> trkSegNHitsTmp;
    std::vector<float> trkSegStartPosXtmp;
    std::vector<float> trkSegStartPosYtmp;
    std::vector<float> trkSegStartPosZtmp;
    std::vector<float> trkSegEndPosXtmp;
    std::vector<float> trkSegEndPosYtmp;
    std::vector<float> trkSegEndPosZtmp;
    std::vector<float> trkSegDirXtmp;
    std::vector<float> trkSegDirYtmp;
    std::vector<float> trkSegDirZtmp;

    Vector_t pcdir;
    for (unsigned int p = 0; p<segradlengths.size(); p++) {
      mcsfitter.linearRegression(track.Trajectory(), breakpoints[p], breakpoints[p+1], pcdir);

      trkSegNHitsTmp.push_back(breakpoints[p+1]-breakpoints[p]);
      trkSegStartPosXtmp.push_back( track.Trajectory().LocationAtPoint(breakpoints[p]).X() );
      trkSegStartPosYtmp.push_back( track.Trajectory().LocationAtPoint(breakpoints[p]).Y() );
      trkSegStartPosZtmp.push_back( track.Trajectory().LocationAtPoint(breakpoints[p]).Z() );
      trkSegEndPosXtmp.push_back( track.Trajectory().LocationAtPoint(breakpoints[p+1]-1).X());
      trkSegEndPosYtmp.push_back( track.Trajectory().LocationAtPoint(breakpoints[p+1]-1).Y());
      trkSegEndPosZtmp.push_back( track.Trajectory().LocationAtPoint(breakpoints[p+1]-1).Z());

      trkSegDirXtmp.push_back(pcdir.X());
      trkSegDirYtmp.push_back(pcdir.Y());
      trkSegDirZtmp.push_back(pcdir.Z());
    }

    trkSegNHits     = trkSegNHitsTmp;
    trkSegStartPosX = trkSegStartPosXtmp;
    trkSegStartPosY = trkSegStartPosYtmp;
    trkSegStartPosZ = trkSegStartPosZtmp;
    trkSegEndPosX   = trkSegEndPosXtmp;
    trkSegEndPosY   = trkSegEndPosYtmp;
    trkSegEndPosZ   = trkSegEndPosZtmp;
    trkSegDirX      = trkSegDirXtmp;
    trkSegDirY      = trkSegDirYtmp;
    trkSegDirZ      = trkSegDirZtmp;

    // remove last 10 cm and repeat the mcs fit
    vector<Point_t> posv_rm10 = track.Trajectory().Trajectory().Positions();
    vector<Vector_t> momv_rm10 = track.Trajectory().Trajectory().Momenta();
    vector<TrajectoryPointFlags> flgv_rm10 = track.Trajectory().Flags();

    size_t init_size_rm10 = posv_rm10.size();
    for (size_t i=0;i<init_size_rm10;++i) {
      if (posv_rm10[init_size_rm10-1-i].X()>-1.) break;
      posv_rm10.pop_back();
      momv_rm10.pop_back();
      flgv_rm10.pop_back();
    }

    init_size_rm10 = posv_rm10.size();
    for (size_t i=0;i<init_size_rm10;++i) {
      if ((track.End()-posv_rm10[init_size_rm10-1-i]).R()>10.) break;
      posv_rm10.pop_back();
      momv_rm10.pop_back();
      flgv_rm10.pop_back();
    }

    if (posv_rm10.size()>4) {
      recob::TrackTrajectory tt_rm10(std::move(posv_rm10),std::move(momv_rm10),std::move(flgv_rm10),track.HasMomentum());
      auto mcsMu_rm10 = mcsfitter.fitMcs(tt_rm10,13);
      trkLength_rm10 = tt_rm10.Length();
      trkMom_MuFwd_rm10    = mcsMu_rm10.fwdMomentum();
      trkMomErr_MuFwd_rm10 = mcsMu_rm10.fwdMomUncertainty();
      trkLL_MuFwd_rm10     = mcsMu_rm10.fwdLogLikelihood();
      trkMom_MuBwd_rm10    = mcsMu_rm10.bwdMomentum();
      trkMomErr_MuBwd_rm10 = mcsMu_rm10.bwdMomUncertainty();
      trkLL_MuBwd_rm10     = mcsMu_rm10.bwdLogLikelihood();
      trkMom_Mu_rm10       = mcsMu_rm10.bestMomentum();
      trkMomErr_Mu_rm10    = mcsMu_rm10.bestMomUncertainty();
      trkLL_Mu_rm10        = mcsMu_rm10.bestLogLikelihood();
      trkDeltaLL_Mu_rm10   = mcsMu_rm10.deltaLogLikelihood();
      trkIsBestFwd_Mu_rm10 = mcsMu_rm10.isBestFwd();
    }

    // remove last 50 cm and repeat the mcs fit
    vector<Point_t> posv_rm50 = track.Trajectory().Trajectory().Positions();
    vector<Vector_t> momv_rm50 = track.Trajectory().Trajectory().Momenta();
    vector<TrajectoryPointFlags> flgv_rm50 = track.Trajectory().Flags();

    size_t init_size_rm50 = posv_rm50.size();
    for (size_t i=0;i<init_size_rm50;++i) {
      if (posv_rm50[init_size_rm50-1-i].X()>-1.) break;
      posv_rm50.pop_back();
      momv_rm50.pop_back();
      flgv_rm50.pop_back();
    }

    init_size_rm50 = posv_rm50.size();
    for (size_t i=0;i<init_size_rm50;++i) {
      if ((track.End()-posv_rm50[init_size_rm50-1-i]).R()>50.) break;
      posv_rm50.pop_back();
      momv_rm50.pop_back();
      flgv_rm50.pop_back();
    }

    if (posv_rm50.size()>4) {
      recob::TrackTrajectory tt_rm50(std::move(posv_rm50),std::move(momv_rm50),std::move(flgv_rm50),track.HasMomentum());
      auto mcsMu_rm50 = mcsfitter.fitMcs(tt_rm50,13);
      trkLength_rm50 = tt_rm50.Length();
      trkMom_MuFwd_rm50    = mcsMu_rm50.fwdMomentum();
      trkMomErr_MuFwd_rm50 = mcsMu_rm50.fwdMomUncertainty();
      trkLL_MuFwd_rm50     = mcsMu_rm50.fwdLogLikelihood();
      trkMom_MuBwd_rm50    = mcsMu_rm50.bwdMomentum();
      trkMomErr_MuBwd_rm50 = mcsMu_rm50.bwdMomUncertainty();
      trkLL_MuBwd_rm50     = mcsMu_rm50.bwdLogLikelihood();
      trkMom_Mu_rm50       = mcsMu_rm50.bestMomentum();
      trkMomErr_Mu_rm50    = mcsMu_rm50.bestMomUncertainty();
      trkLL_Mu_rm50        = mcsMu_rm50.bestLogLikelihood();
      trkDeltaLL_Mu_rm50   = mcsMu_rm50.deltaLogLikelihood();
      trkIsBestFwd_Mu_rm50 = mcsMu_rm50.isBestFwd();
    }

    // remove last 80 cm and repeat the mcs fit
    vector<Point_t> posv_rm80 = track.Trajectory().Trajectory().Positions();
    vector<Vector_t> momv_rm80 = track.Trajectory().Trajectory().Momenta();
    vector<TrajectoryPointFlags> flgv_rm80 = track.Trajectory().Flags();

    size_t init_size_rm80 = posv_rm80.size();
    for (size_t i=0;i<init_size_rm80;++i) {
      if (posv_rm80[init_size_rm80-1-i].X()>-1.) break;
      posv_rm80.pop_back();
      momv_rm80.pop_back();
      flgv_rm80.pop_back();
    }

    init_size_rm80 = posv_rm80.size();
    for (size_t i=0;i<init_size_rm80;++i) {
      if ((track.End()-posv_rm80[init_size_rm80-1-i]).R()>80.) break;
      posv_rm80.pop_back();
      momv_rm80.pop_back();
      flgv_rm80.pop_back();
    }

    if (posv_rm80.size()>4) {
      recob::TrackTrajectory tt_rm80(std::move(posv_rm80),std::move(momv_rm80),std::move(flgv_rm80),track.HasMomentum());
      auto mcsMu_rm80 = mcsfitter.fitMcs(tt_rm80,13);
      trkLength_rm80 = tt_rm80.Length();
      trkMom_MuFwd_rm80    = mcsMu_rm80.fwdMomentum();
      trkMomErr_MuFwd_rm80 = mcsMu_rm80.fwdMomUncertainty();
      trkLL_MuFwd_rm80     = mcsMu_rm80.fwdLogLikelihood();
      trkMom_MuBwd_rm80    = mcsMu_rm80.bwdMomentum();
      trkMomErr_MuBwd_rm80 = mcsMu_rm80.bwdMomUncertainty();
      trkLL_MuBwd_rm80     = mcsMu_rm80.bwdLogLikelihood();
      trkMom_Mu_rm80       = mcsMu_rm80.bestMomentum();
      trkMomErr_Mu_rm80    = mcsMu_rm80.bestMomUncertainty();
      trkLL_Mu_rm80        = mcsMu_rm80.bestLogLikelihood();
      trkDeltaLL_Mu_rm80   = mcsMu_rm80.deltaLogLikelihood();
      trkIsBestFwd_Mu_rm80 = mcsMu_rm80.isBestFwd();
    }

    // MCParticle truth matching
    art::Ptr<simb::MCParticle> pmcp;
    if (simParticles && simHits) {
      pmcp = getAssocMCParticle(*simHits, hits, simMatchPurity);
    }

    // Fill simulated fields
    if (pmcp.isNonnull()) {
      const simb::MCParticle& mcp = *pmcp;

      Vector_t mcstartmom(mcp.Momentum(0).X() * 0.001,
                          mcp.Momentum(0).Y() * 0.001,
                          mcp.Momentum(0).Z() * 0.001);

      Vector_t mcstartdir = mcstartmom.Unit();

      float dotvtx = (track.VertexDirection().X() * mcstartdir.X() +
                      track.VertexDirection().Y() * mcstartdir.Y() +
                      track.VertexDirection().Z() * mcstartdir.Z());

      Point_t mcstart(mcp.Position(0).X(),
                      mcp.Position(0).Y(),
                      mcp.Position(0).Z());

      Point_t mcend(mcp.EndPosition().X(),
                    mcp.EndPosition().Y(),
                    mcp.EndPosition().Z());

      bool mccontained = (mcstart.X() >   30 &&
                          mcstart.X() <  230 &&
                          mcend.X()   >   30 && 
                          mcend.X()   <  230 &&
                          mcstart.Y() >  -85 &&
                          mcstart.Y() <   85 &&
                          mcend.Y()   >  -85 && 
                          mcend.Y()   <   85 &&
                          mcstart.Z() >   30 && 
                          mcstart.Z() < 1010 && 
                          mcend.Z()   >   30 &&
                          mcend.Z()   < 1010);

      std::vector<Point_t> positions;
      std::vector<Vector_t> momenta;

      float thislen = 0;
      std::vector<float> lengths;
      std::vector<float> moms;
      std::vector<float> betas;
      std::vector<float> rmss;
      std::vector<float> rmsshl;
      std::vector<float> es;
      std::vector<float> rmoms;
      std::vector<float> rbetas;
      std::vector<float> res;
      std::vector<float> rmsRs;
      std::vector<Vector_t> startsegdirs;

      bool newseg = true;
      const float thisSegLen = 14.;
      const float m = mcsfittermc.mass(13);

      float mclen = 0;
      const simb::MCTrajectory& mctrack = mcp.Trajectory();

      std::vector<float> simMcStepPosXtmp;
      std::vector<float> simMcStepPosYtmp;
      std::vector<float> simMcStepPosZtmp;

      for (size_t imc=0; imc<mctrack.size(); imc++) {
        if (mctrack.Position(imc).X() <    0 || mctrack.Position(imc).X() >  260) continue;
        if (mctrack.Position(imc).Y() < -115 || mctrack.Position(imc).Y() >  115) continue;
        if (mctrack.Position(imc).Z() <    0 || mctrack.Position(imc).Z() > 1040) continue;

        positions.push_back({
          mctrack.Position(imc).X(), mctrack.Position(imc).Y(), mctrack.Position(imc).Z()
        });

        simMcStepPosXtmp.push_back(mctrack.Position(imc).X());
        simMcStepPosYtmp.push_back(mctrack.Position(imc).Y());
        simMcStepPosZtmp.push_back(mctrack.Position(imc).Z());

        Vector_t mcmom(mctrack.Momentum(imc).X() * 0.001,
                       mctrack.Momentum(imc).Y() * 0.001,
                       mctrack.Momentum(imc).Z() * 0.001);
        momenta.push_back(mcmom);
        Vector_t mcdir = mcmom.Unit();

        if (newseg) {
          startsegdirs.push_back(mcdir);
          newseg = false;
          moms.push_back(mcmom.R());
          float beta = sqrt( 1. - ((m*m)/(moms.back()*moms.back() + m*m)) );
          betas.push_back(beta);
          es.push_back(mctrack.Momentum(imc).E()*0.001);

          float e = mcsfittermc.GetE(std::sqrt(mcstartmom.R()*mcstartmom.R() - m*m), mclen, m);
          float p = std::sqrt(e*e - m*m);
          rmoms.push_back(p);
          float betaR = sqrt( 1. - ((m*m)/(rmoms.back()*rmoms.back() + m*m)) );
          rbetas.push_back(betaR);
          res.push_back(e);
        }

        if (imc < (mctrack.size() - 0.5)) {
          if (!(mctrack.Position(imc+1).X() <    0 || mctrack.Position(imc+1).X() >  260 ||
                mctrack.Position(imc+1).Y() < -115 || mctrack.Position(imc+1).Y() >  115 ||
                mctrack.Position(imc+1).Z() <    0 || mctrack.Position(imc+1).Z() > 1040)) {
            float dx = mctrack.Position(imc+1).X() - mctrack.Position(imc).X();
            float dy = mctrack.Position(imc+1).Y() - mctrack.Position(imc).Y();
            float dz = mctrack.Position(imc+1).Z() - mctrack.Position(imc).Z();
            mclen += sqrt(dx*dx + dy*dy + dz*dz);

            // check segment length along the initial direction
            thislen += startsegdirs.back().Dot(Vector_t(dx, dy, dz));
          }
        }

        if (thislen > thisSegLen-1.0) {
          lengths.push_back(thislen);

          // momentum is before the segment, but rms needs the actual length so goes here
          float beta = betas.back();
          float rms = ( mcsfittermc.HighlandFirstTerm(moms.back()) / (moms.back()*beta) ) * ( 1.0 + 0.038 * std::log( thislen/14. ) ) * sqrt( thislen/14. );
          rmss.push_back(rms);
          float rmshl = ( 13.6 / (moms.back()*beta) ) * ( 1.0 + 0.038 * std::log( thislen/14. ) ) * sqrt( thislen/14. );
          rmsshl.push_back(rmshl);

          float betaR = rbetas.back();
          float rmsR = ( mcsfittermc.HighlandFirstTerm(rmoms.back()) / (rmoms.back()*betaR) ) * ( 1.0 + 0.038 * std::log( thislen/14. ) ) * sqrt( thislen/14. );
          rmsRs.push_back(rmsR);
          thislen = 0;
          newseg = true;
        }

        // update values until the last point
        simEndPosX = mctrack.Position(imc).X();
        simEndPosY = mctrack.Position(imc).Y();
        simEndPosZ = mctrack.Position(imc).Z();
        simEndDirX = mcdir.X();
        simEndDirY = mcdir.Y();
        simEndDirZ = mcdir.Z();
        simMomEnd = mcmom.R();
      }

      if (positions.size() < 2) continue;
      simMom = mcstartmom.R();
      simLength = mclen;
      simStartPosX = mcstart.X();
      simStartPosY = mcstart.Y();
      simStartPosZ = mcstart.Z();
      simStartDirX = mcstartdir.X();
      simStartDirY = mcstartdir.Y();
      simStartDirZ = mcstartdir.Z();
      simMcStepPosX = simMcStepPosXtmp;
      simMcStepPosY = simMcStepPosYtmp;
      simMcStepPosZ = simMcStepPosZtmp;
      simID = std::abs(mcp.PdgCode());
      simProc = mcp.Process();
      simIsContained = mccontained;
      simAndTrkSameDir = dotvtx > 0;

      recob::TrackTrajectory::Flags_t flags(positions.size());
      const recob::TrackTrajectory mctj(std::move(positions),std::move(momenta),std::move(flags),true);
      const auto mcsmcMu = mcsfittermc.fitMcs(mctj, 13);
      simMom_MuFwd    = mcsmcMu.fwdMomentum();
      simMomErr_MuFwd = mcsmcMu.fwdMomUncertainty();
      simLL_MuFwd     = mcsmcMu.fwdLogLikelihood();
      simMom_MuBwd    = mcsmcMu.bwdMomentum();
      simMomErr_MuBwd = mcsmcMu.bwdMomUncertainty();
      simLL_MuBwd     = mcsmcMu.bwdLogLikelihood();
      simMom_Mu       = mcsmcMu.bestMomentum();
      simMomErr_Mu    = mcsmcMu.bestMomUncertainty();
      simLL_Mu        = mcsmcMu.bestLogLikelihood();
      simDeltaLL_Mu   = mcsmcMu.deltaLogLikelihood();
      simIsBestFwd_Mu = mcsmcMu.isBestFwd();

      simSegRadLengths = mcsmcMu.segmentRadLengths();
      simScattAngles   = mcsmcMu.scatterAngles();
      simL = lengths;
      simE = es;
      simP = moms;
      simBeta = betas;
      simRE = res;
      simRP = rmoms;
      simRBeta = rbetas;
      simRMS = rmss;
      simRMSHL = rmsshl;
      simRMSR = rmsRs;

      std::vector<float> dtheta3d;
      std::vector<float> dtheta2dxz;
      std::vector<float> dtheta2dyz;
      for (size_t p=1; p<startsegdirs.size(); p++) {
        auto pcdir0 = startsegdirs[p-1];
        auto pcdir1 = startsegdirs[p];

        const float cosval3d = pcdir0.X()*pcdir1.X()+pcdir0.Y()*pcdir1.Y()+pcdir0.Z()*pcdir1.Z();

        // units are mrad
        float dt3d = 1000.*acos(cosval3d);  // should we try to use expansion for small angles?
        dtheta3d.push_back(dt3d);

        // this is the new basis with pcdir along z
        auto pcdir0z = pcdir0;
        auto pcdir0y = pcdir0.Cross(Vector_t(0,0,1)).Unit();
        auto pcdir0x = pcdir0y.Cross(pcdir0z).Unit();

        // now find projections of pcdir1 along the axes of the new basis
        auto pcdir1x = pcdir0x.Dot(pcdir1);
        auto pcdir1y = pcdir0y.Dot(pcdir1);
        auto pcdir1z = pcdir0z.Dot(pcdir1);

        // compute angles with respect to z axis
        dtheta2dxz.push_back(1000.*atan2(pcdir1x,pcdir1z));
        dtheta2dyz.push_back(1000.*atan2(pcdir1y,pcdir1z));
      }

      simA3d = dtheta3d;
      simA2dxz = dtheta2dxz;
      simA2dyz = dtheta2dyz;

      vector<size_t> breakpoints_mc;
      vector<float> segradlengths_mc;
      vector<float> cumseglens_mc;
      mcsfittermc.breakTrajInSegments(mctj, breakpoints_mc, segradlengths_mc, cumseglens_mc);
      Vector_t pcdir_mc;

      std::vector<Vector_t> segdirs_mc;
      for (size_t p=0; p<segradlengths_mc.size(); p++) {
        mcsfittermc.linearRegression(mctj, breakpoints_mc[p], breakpoints_mc[p+1], pcdir_mc);
        segdirs_mc.push_back(pcdir_mc);
      }

      std::vector<float> dthetaAvg3d;
      std::vector<float> dthetaAvg2dxz;
      std::vector<float> dthetaAvg2dyz;
      for (size_t p=1; p<segdirs_mc.size(); p++) {
        auto pcdir_mc0 = segdirs_mc[p-1];
        auto pcdir_mc1 = segdirs_mc[p];

        const float cosval3d = pcdir_mc0.X()*pcdir_mc1.X()+pcdir_mc0.Y()*pcdir_mc1.Y()+pcdir_mc0.Z()*pcdir_mc1.Z();

        //units are mrad
        float dt3d = 1000.*acos(cosval3d);  // should we try to use expansion for small angles?
        dthetaAvg3d.push_back(dt3d);

        // this is the new basis with pcdir_mc along z
        auto pcdir_mc0z = pcdir_mc0;
        auto pcdir_mc0y = pcdir_mc0.Cross(Vector_t(0,0,1)).Unit();
        auto pcdir_mc0x = pcdir_mc0y.Cross(pcdir_mc0z).Unit();
        //
        // now find projections of pcdir_mc1 along the axes of the new basis
        auto pcdir_mc1x = pcdir_mc0x.Dot(pcdir_mc1);
        auto pcdir_mc1y = pcdir_mc0y.Dot(pcdir_mc1);
        auto pcdir_mc1z = pcdir_mc0z.Dot(pcdir_mc1);
        //
        // compute angles with respect to z axis
        dthetaAvg2dxz.push_back(1000.*atan2(pcdir_mc1x,pcdir_mc1z));
        dthetaAvg2dyz.push_back(1000.*atan2(pcdir_mc1y,pcdir_mc1z));
      }

      simAAvg3d = dthetaAvg3d;
      simAAvg2dxz = dthetaAvg2dxz;
      simAAvg2dyz = dthetaAvg2dyz;
    }

    tree->Fill();
  }

  delete simParticles;
  delete simHits;
}

bool TrajectoryMCSNtuple::isPfpFromNeutrino(
    art::ValidHandle<std::vector<recob::PFParticle> > pfph, size_t ipfp) {
  const recob::PFParticle& pf = pfph->at(ipfp);
  if (pf.IsPrimary() && std::abs(pf.PdgCode()) != 14) return false;
  if (pf.IsPrimary()) return true;
  size_t parentid = pf.Parent();

  for (size_t iPF=0; iPF<pfph->size(); iPF++) {
    const recob::PFParticle& parent = pfph->at(iPF);
    if (parent.Self() != parentid) continue;
    return isPfpFromNeutrino(pfph, iPF);
  }

  std::cout << "argh you should never end up here... "
            << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << std::endl;
  return false;
}


// BackTrack a single hit collection (i.e. from a PFParticle)
// https://github.com/ubneutrinos/searchingfornues/blob/develop/Selection/CommonDefs/BacktrackingFuncs.h#L38-L81
art::Ptr<simb::MCParticle> TrajectoryMCSNtuple::getAssocMCParticle(
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &hittruth,
    const std::vector<art::Ptr<recob::Hit> >& hits, float& purity) {
  // store total charge from hits
  float pfpcharge = 0;  // total hit charge from clusters
  float maxcharge = 0;  // charge backtracked to best match

  // credit: Wes Ketchum
  std::unordered_map<int, double> trkide;
  std::unordered_map<int, float> trkq;
  double maxe = -1, tote = 0;
  art::Ptr<simb::MCParticle> maxp_me;  // pointer for the particle match we will calculate

  for (auto h : hits) {
    pfpcharge += h->Integral();
    std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth.at(h.key());
    std::vector<anab::BackTrackerHitMatchingData const *> match_vec = hittruth.data(h.key());

    // loop over particles
    for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
      trkide[particle_vec[i_p]->TrackId()] += match_vec[i_p]->energy;
      trkq[particle_vec[i_p]->TrackId()] += h->Integral() * match_vec[i_p]->ideFraction;
      tote += match_vec[i_p]->energy;
      if (trkide[particle_vec[i_p]->TrackId()] > maxe) {
        // keep track of maximum
        maxe = trkide[particle_vec[i_p]->TrackId()];
        maxp_me = particle_vec[i_p];
        maxcharge = trkq[particle_vec[i_p]->TrackId()];
      }
    }
  }

  purity = maxcharge / pfpcharge;

  return maxp_me;
}

DEFINE_ART_MODULE(TrajectoryMCSNtuple)

