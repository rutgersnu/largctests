////////////////////////////////////////////////////////////////////////
// Class:       SplitTrackNtuplizer
// Plugin Type: analyzer (art v2_05_00)
// File:        SplitTrackNtuplizer_module.cc
//
// Generated at Mon Feb  6 10:06:04 2017 by Giuseppe Cerati using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcore/Geometry/Geometry.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingPlane.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"

#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "lardata/RecoObjects/TrackingPlaneHelper.h"
#include "lardata/RecoObjects/TrackState.h"
#include "lardata/RecoObjects/TrackStatePropagator.h"

#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

class SplitTrackNtuplizer;

class SplitTrackNtuplizer : public art::EDAnalyzer {
public:
  explicit SplitTrackNtuplizer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SplitTrackNtuplizer(SplitTrackNtuplizer const &) = delete;
  SplitTrackNtuplizer(SplitTrackNtuplizer &&) = delete;
  SplitTrackNtuplizer & operator = (SplitTrackNtuplizer const &) = delete;
  SplitTrackNtuplizer & operator = (SplitTrackNtuplizer &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  void beginJob() override;
  void resetTree();
private:
  //
  std::string inputPFLabel;
  std::string inputTracksLabel;
  std::string inputTracksLabel1st;
  std::string inputTracksLabel2nd;
  bool doSplitVertices;
  //
  trkf::TrajectoryMCSFitter mcsFitMu;
  trkf::TrajectoryMCSFitter mcsFitP;
  //
  TTree* tree;
  //
  // Event
  //
  int    run, subrun, eventid;
  int passSelII;
  //
  // Original track
  //
  int   tk_id     ;
  int   tk_contain;
  int   tk_nhits  ;
  int   tk_nvhits ;
  int   tk_nhits_u;
  int   tk_nhits_v;
  int   tk_nhits_y;
  int   tk_ndof   ;
  float tk_chi2   ;
  float tk_length ;
  float tk_drvtx   ;
  float tk_mommumcs;
  float tk_mompmcs ;
  float tk_mommumcserr;
  float tk_mompmcserr ;
  float tk_mommumcsdll;
  float tk_mompmcsdll ;
  float tk_mommurng;
  float tk_momprng ;
  std::vector<float> tk_vtx_gpar  ;
  std::vector<float> tk_vtx_lpar  ;
  std::vector<float> tk_vtx_gcov  ;
  std::vector<float> tk_vtx_lcov  ;
  std::vector<float> tk_end_gpar  ;
  std::vector<float> tk_end_lpar  ;
  std::vector<float> tk_end_gcov  ;
  std::vector<float> tk_end_lcov  ;
  std::vector<float> tk_mid_gpar  ;
  // std::vector<float> tk_mid_lpar  ;
  // std::vector<float> tk_mid_gcov  ;
  // std::vector<float> tk_mid_lcov  ;
  std::vector<float> tk_vtxmc_gpar;
  std::vector<float> tk_vtxmc_lpar;
  std::vector<float> tk_vtxmc_gcov;
  std::vector<float> tk_vtxmc_lcov;
  std::vector<float> tk_endmc_gpar;
  std::vector<float> tk_endmc_lpar;
  std::vector<float> tk_endmc_gcov;
  std::vector<float> tk_endmc_lcov;
  std::vector<float> tk_midmc_gpar;
  // std::vector<float> tk_midmc_lpar;
  // std::vector<float> tk_midmc_gcov;
  // std::vector<float> tk_midmc_lcov;
  std::vector<float> tk_hits_x    ;
  std::vector<float> tk_hits_rms  ;
  std::vector<int>   tk_hits_plane;
  std::vector<int>   tk_hits_wire ;
  std::vector<float> vx_pos;
  int vx_ntks;
  int vx_tkid;
  //
  // MC track
  //
  int   mc_pid    ;
  int   mc_mpid   ;
  int   mc_contain;
  float mc_length ;
  float mc_mom    ;
  std::vector<float> mc_vtx_gpar  ;
  std::vector<float> mc_end_gpar  ;
  std::vector<float> mc_vtxtk_gpar;
  std::vector<float> mc_endtk_gpar;
  std::vector<float> mc_midtk_gpar;
  std::vector<float> mc_vtxtk1_gpar;
  std::vector<float> mc_endtk1_gpar;
  std::vector<float> mc_vtxtk2_gpar;
  std::vector<float> mc_endtk2_gpar;
  //
  // Split track 1
  //
  int   tk1_id     ;
  int   tk1_nhits  ;
  int   tk1_nvhits ;
  int   tk1_nhits_u;
  int   tk1_nhits_v;
  int   tk1_nhits_y;
  int   tk1_ndof   ;
  float tk1_chi2   ;
  float tk1_length ;
  std::vector<float> tk1_vtx_gpar  ;
  std::vector<float> tk1_vtx_lpar  ;
  std::vector<float> tk1_vtx_gcov  ;
  std::vector<float> tk1_vtx_lcov  ;
  std::vector<float> tk1_end_gpar  ;
  std::vector<float> tk1_end_lpar  ;
  std::vector<float> tk1_end_gcov  ;
  std::vector<float> tk1_end_lcov  ;
  std::vector<float> tk1_vtxmc_gpar;
  std::vector<float> tk1_vtxmc_lpar;
  std::vector<float> tk1_vtxmc_gcov;
  std::vector<float> tk1_vtxmc_lcov;
  std::vector<float> tk1_midmc_gpar;
  std::vector<float> tk1_endmc_gpar;
  std::vector<float> tk1_endmc_lpar;
  std::vector<float> tk1_endmc_gcov;
  std::vector<float> tk1_endmc_lcov;
  std::vector<float> tk1_vtxtk_gpar;
  std::vector<float> tk1_vtxtk_lpar;
  std::vector<float> tk1_vtxtk_gcov;
  std::vector<float> tk1_vtxtk_lcov;
  std::vector<float> tk1_midtk_gpar;
  std::vector<float> tk1_midtk_lpar;
  std::vector<float> tk1_midtk_gcov;
  std::vector<float> tk1_midtk_lcov;
  std::vector<float> tk1_endtk_gpar;
  std::vector<float> tk1_endtk_lpar;
  std::vector<float> tk1_endtk_gcov;
  std::vector<float> tk1_endtk_lcov;
  std::vector<float> tk1_hits_x    ;
  std::vector<float> tk1_hits_rms  ;
  std::vector<int>   tk1_hits_plane;
  std::vector<int>   tk1_hits_wire ;
  std::vector<float> vx1_pos;
  int vx1_ntks;
  int vx1_tkid;
  float vx1_dist;
  //
  // Split track 2
  //
  int   tk2_id     ;
  int   tk2_nhits  ;
  int   tk2_nvhits ;
  int   tk2_nhits_u;
  int   tk2_nhits_v;
  int   tk2_nhits_y;
  int   tk2_ndof   ;
  float tk2_chi2   ;
  float tk2_length ;
  std::vector<float> tk2_vtx_gpar  ;
  std::vector<float> tk2_vtx_lpar  ;
  std::vector<float> tk2_vtx_gcov  ;
  std::vector<float> tk2_vtx_lcov  ;
  std::vector<float> tk2_end_gpar  ;
  std::vector<float> tk2_end_lpar  ;
  std::vector<float> tk2_end_gcov  ;
  std::vector<float> tk2_end_lcov  ;
  std::vector<float> tk2_vtxmc_gpar;
  std::vector<float> tk2_vtxmc_lpar;
  std::vector<float> tk2_vtxmc_gcov;
  std::vector<float> tk2_vtxmc_lcov;
  std::vector<float> tk2_midmc_gpar;
  std::vector<float> tk2_endmc_gpar;
  std::vector<float> tk2_endmc_lpar;
  std::vector<float> tk2_endmc_gcov;
  std::vector<float> tk2_endmc_lcov;
  std::vector<float> tk2_vtxtk_gpar;
  std::vector<float> tk2_vtxtk_lpar;
  std::vector<float> tk2_vtxtk_gcov;
  std::vector<float> tk2_vtxtk_lcov;
  std::vector<float> tk2_midtk_gpar;
  std::vector<float> tk2_midtk_lpar;
  std::vector<float> tk2_midtk_gcov;
  std::vector<float> tk2_midtk_lcov;
  std::vector<float> tk2_endtk_gpar;
  std::vector<float> tk2_endtk_lpar;
  std::vector<float> tk2_endtk_gcov;
  std::vector<float> tk2_endtk_lcov;
  std::vector<float> tk2_vtxtk1_gpar;
  std::vector<float> tk2_vtxtk1_lpar;
  std::vector<float> tk2_vtxtk1_gcov;
  std::vector<float> tk2_vtxtk1_lcov;
  std::vector<float> tk2_endtk1_gpar;
  std::vector<float> tk2_endtk1_lpar;
  std::vector<float> tk2_endtk1_gcov;
  std::vector<float> tk2_endtk1_lcov;
  std::vector<float> tk2_hits_x    ;
  std::vector<float> tk2_hits_rms  ;
  std::vector<int>   tk2_hits_plane;
  std::vector<int>   tk2_hits_wire ;
  std::vector<float> vx2_pos;
  int vx2_ntks;
  int vx2_tkid;
  float vx2_dist;
};

void SplitTrackNtuplizer::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  //
  tree = tfs->make<TTree>("tree", "tree");
  //
  // Event
  //
  tree->Branch("run", &run,"run/I");
  tree->Branch("subrun", &subrun, "subrun/I");
  tree->Branch("eventid", &eventid, "eventid/I");
  tree->Branch("passSelII", &passSelII, "passSelII/I");
  //
  // Original track
  //
  tree->Branch("tk_id"     , &tk_id     , "tk_id/I"     );
  tree->Branch("tk_contain", &tk_contain, "tk_contain/I");
  tree->Branch("tk_nhits"  , &tk_nhits  , "tk_nhits/I"  );
  tree->Branch("tk_nvhits" , &tk_nvhits , "tk_nvhits/I" );
  tree->Branch("tk_nhits_u", &tk_nhits_u, "tk_nhits_u/I");
  tree->Branch("tk_nhits_v", &tk_nhits_v, "tk_nhits_v/I");
  tree->Branch("tk_nhits_y", &tk_nhits_y, "tk_nhits_y/I");
  tree->Branch("tk_ndof"   , &tk_ndof   , "tk_ndof/I"   );
  tree->Branch("tk_chi2"   , &tk_chi2   , "tk_chi2/F"   );
  tree->Branch("tk_length" , &tk_length , "tk_length/F" );
  tree->Branch("tk_drvtx"    , &tk_drvtx    , "tk_drvtx/F"    );
  tree->Branch("tk_mommumcs" , &tk_mommumcs , "tk_mommumcs/F" );
  tree->Branch("tk_mompmcs"  , &tk_mompmcs  , "tk_mompmcs/F"  );
  tree->Branch("tk_mommumcserr" , &tk_mommumcserr , "tk_mommumcserr/F" );
  tree->Branch("tk_mompmcserr"  , &tk_mompmcserr  , "tk_mompmcserr/F"  );
  tree->Branch("tk_mommumcsdll" , &tk_mommumcsdll , "tk_mommumcsdll/F" );
  tree->Branch("tk_mompmcsdll"  , &tk_mompmcsdll  , "tk_mompmcsdll/F"  );
  tree->Branch("tk_mommurng" , &tk_mommurng , "tk_mommurng/F" );
  tree->Branch("tk_momprng"  , &tk_momprng  , "tk_momprng/F"  );
  tree->Branch("tk_vtx_gpar"  , &tk_vtx_gpar  );
  tree->Branch("tk_vtx_lpar"  , &tk_vtx_lpar  );
  tree->Branch("tk_vtx_gcov"  , &tk_vtx_gcov  );
  tree->Branch("tk_vtx_lcov"  , &tk_vtx_lcov  );
  tree->Branch("tk_end_gpar"  , &tk_end_gpar  );
  tree->Branch("tk_end_lpar"  , &tk_end_lpar  );
  tree->Branch("tk_end_gcov"  , &tk_end_gcov  );
  tree->Branch("tk_end_lcov"  , &tk_end_lcov  );
  tree->Branch("tk_mid_gpar"  , &tk_mid_gpar  );
  // tree->Branch("tk_mid_lpar", &tk_mid_lpar  );
  // tree->Branch("tk_mid_gcov", &tk_mid_gcov  );
  // tree->Branch("tk_mid_lcov", &tk_mid_lcov  );
  tree->Branch("tk_vtxmc_gpar", &tk_vtxmc_gpar);
  tree->Branch("tk_vtxmc_lpar", &tk_vtxmc_lpar);
  tree->Branch("tk_vtxmc_gcov", &tk_vtxmc_gcov);
  tree->Branch("tk_vtxmc_lcov", &tk_vtxmc_lcov);
  tree->Branch("tk_endmc_gpar", &tk_endmc_gpar);
  tree->Branch("tk_endmc_lpar", &tk_endmc_lpar);
  tree->Branch("tk_endmc_gcov", &tk_endmc_gcov);
  tree->Branch("tk_endmc_lcov", &tk_endmc_lcov);
  tree->Branch("tk_midmc_gpar", &tk_midmc_gpar);
  // tree->Branch("tk_midmc_lpar", &tk_midmc_lpar);
  // tree->Branch("tk_midmc_gcov", &tk_midmc_gcov);
  // tree->Branch("tk_midmc_lcov", &tk_midmc_lcov);
  tree->Branch("tk_hits_x"    , &tk_hits_x    );
  tree->Branch("tk_hits_rms"  , &tk_hits_rms  );
  tree->Branch("tk_hits_plane", &tk_hits_plane);
  tree->Branch("tk_hits_wire" , &tk_hits_wire );
  tree->Branch("vx_pos"  , &vx_pos );
  tree->Branch("vx_ntks", &vx_ntks, "vx_ntks/I");
  tree->Branch("vx_tkid", &vx_tkid, "vx_tkid/I");
  //
  // MC track
  //
  tree->Branch("mc_pid"    , &mc_pid    ,"mc_pid/I"    );
  tree->Branch("mc_mpid"   , &mc_mpid   ,"mc_mpid/I"   );
  tree->Branch("mc_contain", &mc_contain,"mc_contain/I");
  tree->Branch("mc_length" , &mc_length ,"mc_length/F" );
  tree->Branch("mc_mom"    , &mc_mom    ,"mc_mom/F"    );
  tree->Branch("mc_vtx_gpar"  , &mc_vtx_gpar  );
  tree->Branch("mc_end_gpar"  , &mc_end_gpar  );
  tree->Branch("mc_vtxtk_gpar", &mc_vtxtk_gpar);
  tree->Branch("mc_endtk_gpar", &mc_endtk_gpar);
  tree->Branch("mc_midtk_gpar", &mc_midtk_gpar);
  tree->Branch("mc_vtxtk1_gpar", &mc_vtxtk1_gpar);
  tree->Branch("mc_endtk1_gpar", &mc_endtk1_gpar);
  tree->Branch("mc_vtxtk2_gpar", &mc_vtxtk2_gpar);
  tree->Branch("mc_endtk2_gpar", &mc_endtk2_gpar);
  //
  // Split track 1
  //
  tree->Branch("tk1_id"     , &tk1_id     , "tk1_id/I"     );
  tree->Branch("tk1_nhits"  , &tk1_nhits  , "tk1_nhits/I"  );
  tree->Branch("tk1_nvhits" , &tk1_nvhits , "tk1_nvhits/I" );
  tree->Branch("tk1_nhits_u", &tk1_nhits_u, "tk1_nhits_u/I");
  tree->Branch("tk1_nhits_v", &tk1_nhits_v, "tk1_nhits_v/I");
  tree->Branch("tk1_nhits_y", &tk1_nhits_y, "tk1_nhits_y/I");
  tree->Branch("tk1_ndof"   , &tk1_ndof   , "tk1_ndof/I"   );
  tree->Branch("tk1_chi2"   , &tk1_chi2   , "tk1_chi2/F"   );
  tree->Branch("tk1_length" , &tk1_length , "tk1_length/F" );
  tree->Branch("tk1_vtx_gpar"  , &tk1_vtx_gpar  );
  tree->Branch("tk1_vtx_lpar"  , &tk1_vtx_lpar  );
  tree->Branch("tk1_vtx_gcov"  , &tk1_vtx_gcov  );
  tree->Branch("tk1_vtx_lcov"  , &tk1_vtx_lcov  );
  tree->Branch("tk1_end_gpar"  , &tk1_end_gpar  );
  tree->Branch("tk1_end_lpar"  , &tk1_end_lpar  );
  tree->Branch("tk1_end_gcov"  , &tk1_end_gcov  );
  tree->Branch("tk1_end_lcov"  , &tk1_end_lcov  );
  tree->Branch("tk1_vtxmc_gpar", &tk1_vtxmc_gpar);
  tree->Branch("tk1_vtxmc_lpar", &tk1_vtxmc_lpar);
  tree->Branch("tk1_vtxmc_gcov", &tk1_vtxmc_gcov);
  tree->Branch("tk1_vtxmc_lcov", &tk1_vtxmc_lcov);
  tree->Branch("tk1_midmc_gpar", &tk1_midmc_gpar);
  tree->Branch("tk1_endmc_gpar", &tk1_endmc_gpar);
  tree->Branch("tk1_endmc_lpar", &tk1_endmc_lpar);
  tree->Branch("tk1_endmc_gcov", &tk1_endmc_gcov);
  tree->Branch("tk1_endmc_lcov", &tk1_endmc_lcov);
  tree->Branch("tk1_vtxtk_gpar", &tk1_vtxtk_gpar);
  tree->Branch("tk1_vtxtk_lpar", &tk1_vtxtk_lpar);
  tree->Branch("tk1_vtxtk_gcov", &tk1_vtxtk_gcov);
  tree->Branch("tk1_vtxtk_lcov", &tk1_vtxtk_lcov);
  tree->Branch("tk1_midtk_gpar", &tk1_midtk_gpar);
  tree->Branch("tk1_midtk_lpar", &tk1_midtk_lpar);
  tree->Branch("tk1_midtk_gcov", &tk1_midtk_gcov);
  tree->Branch("tk1_midtk_lcov", &tk1_midtk_lcov);
  tree->Branch("tk1_endtk_gpar", &tk1_endtk_gpar);
  tree->Branch("tk1_endtk_lpar", &tk1_endtk_lpar);
  tree->Branch("tk1_endtk_gcov", &tk1_endtk_gcov);
  tree->Branch("tk1_endtk_lcov", &tk1_endtk_lcov);
  tree->Branch("tk1_hits_x"    , &tk1_hits_x    );
  tree->Branch("tk1_hits_rms"  , &tk1_hits_rms  );
  tree->Branch("tk1_hits_plane", &tk1_hits_plane);
  tree->Branch("tk1_hits_wire" , &tk1_hits_wire );
  tree->Branch("vx1_pos" , &vx1_pos );
  tree->Branch("vx1_ntks", &vx1_ntks, "vx1_ntks/I");
  tree->Branch("vx1_tkid", &vx1_tkid, "vx1_tkid/I");
  tree->Branch("vx1_dist", &vx1_dist, "vx1_dist/F");
  //
  // Split track 2
  //
  tree->Branch("tk2_id"     , &tk2_id     , "tk2_id/I"     );
  tree->Branch("tk2_nhits"  , &tk2_nhits  , "tk2_nhits/I"  );
  tree->Branch("tk2_nvhits" , &tk2_nvhits , "tk2_nvhits/I" );
  tree->Branch("tk2_nhits_u", &tk2_nhits_u, "tk2_nhits_u/I");
  tree->Branch("tk2_nhits_v", &tk2_nhits_v, "tk2_nhits_v/I");
  tree->Branch("tk2_nhits_y", &tk2_nhits_y, "tk2_nhits_y/I");
  tree->Branch("tk2_ndof"   , &tk2_ndof   , "tk2_ndof/I"   );
  tree->Branch("tk2_chi2"   , &tk2_chi2   , "tk2_chi2/F"   );
  tree->Branch("tk2_length" , &tk2_length , "tk2_length/F" );
  tree->Branch("tk2_vtx_gpar"  , &tk2_vtx_gpar  );
  tree->Branch("tk2_vtx_lpar"  , &tk2_vtx_lpar  );
  tree->Branch("tk2_vtx_gcov"  , &tk2_vtx_gcov  );
  tree->Branch("tk2_vtx_lcov"  , &tk2_vtx_lcov  );
  tree->Branch("tk2_end_gpar"  , &tk2_end_gpar  );
  tree->Branch("tk2_end_lpar"  , &tk2_end_lpar  );
  tree->Branch("tk2_end_gcov"  , &tk2_end_gcov  );
  tree->Branch("tk2_end_lcov"  , &tk2_end_lcov  );
  tree->Branch("tk2_vtxmc_gpar", &tk2_vtxmc_gpar);
  tree->Branch("tk2_vtxmc_lpar", &tk2_vtxmc_lpar);
  tree->Branch("tk2_vtxmc_gcov", &tk2_vtxmc_gcov);
  tree->Branch("tk2_vtxmc_lcov", &tk2_vtxmc_lcov);
  tree->Branch("tk2_midmc_gpar", &tk2_midmc_gpar);
  tree->Branch("tk2_endmc_gpar", &tk2_endmc_gpar);
  tree->Branch("tk2_endmc_lpar", &tk2_endmc_lpar);
  tree->Branch("tk2_endmc_gcov", &tk2_endmc_gcov);
  tree->Branch("tk2_endmc_lcov", &tk2_endmc_lcov);
  tree->Branch("tk2_vtxtk_gpar", &tk2_vtxtk_gpar);
  tree->Branch("tk2_vtxtk_lpar", &tk2_vtxtk_lpar);
  tree->Branch("tk2_vtxtk_gcov", &tk2_vtxtk_gcov);
  tree->Branch("tk2_vtxtk_lcov", &tk2_vtxtk_lcov);
  tree->Branch("tk2_midtk_gpar", &tk2_midtk_gpar);
  tree->Branch("tk2_midtk_lpar", &tk2_midtk_lpar);
  tree->Branch("tk2_midtk_gcov", &tk2_midtk_gcov);
  tree->Branch("tk2_midtk_lcov", &tk2_midtk_lcov);
  tree->Branch("tk2_endtk_gpar", &tk2_endtk_gpar);
  tree->Branch("tk2_endtk_lpar", &tk2_endtk_lpar);
  tree->Branch("tk2_endtk_gcov", &tk2_endtk_gcov);
  tree->Branch("tk2_endtk_lcov", &tk2_endtk_lcov);
  tree->Branch("tk2_vtxtk1_gpar", &tk2_vtxtk1_gpar);
  tree->Branch("tk2_vtxtk1_lpar", &tk2_vtxtk1_lpar);
  tree->Branch("tk2_vtxtk1_gcov", &tk2_vtxtk1_gcov);
  tree->Branch("tk2_vtxtk1_lcov", &tk2_vtxtk1_lcov);
  tree->Branch("tk2_endtk1_gpar", &tk2_endtk1_gpar);
  tree->Branch("tk2_endtk1_lpar", &tk2_endtk1_lpar);
  tree->Branch("tk2_endtk1_gcov", &tk2_endtk1_gcov);
  tree->Branch("tk2_endtk1_lcov", &tk2_endtk1_lcov);
  tree->Branch("tk2_hits_x"    , &tk2_hits_x    );
  tree->Branch("tk2_hits_rms"  , &tk2_hits_rms  );
  tree->Branch("tk2_hits_plane", &tk2_hits_plane);
  tree->Branch("tk2_hits_wire" , &tk2_hits_wire );
  tree->Branch("vx2_pos" , &vx2_pos );
  tree->Branch("vx2_ntks", &vx2_ntks, "vx2_ntks/I");
  tree->Branch("vx2_tkid", &vx2_tkid, "vx2_tkid/I");
  tree->Branch("vx2_dist", &vx2_dist, "vx2_dist/F");
}

SplitTrackNtuplizer::SplitTrackNtuplizer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  inputPFLabel(p.get<std::string>("inputPFLabel","pandora")),
  inputTracksLabel(p.get<std::string>("inputTracksLabel")),
  inputTracksLabel1st(p.get<std::string>("inputTracksLabel1st")),
  inputTracksLabel2nd(p.get<std::string>("inputTracksLabel2nd")),
  doSplitVertices(p.get<bool>("doSplitVertices")),
  mcsFitMu(fhicl::Table<trkf::TrajectoryMCSFitter::Config>(p.get<fhicl::ParameterSet>("mcsfitmu"))),
  mcsFitP(fhicl::Table<trkf::TrajectoryMCSFitter::Config>(p.get<fhicl::ParameterSet>("mcsfitp")))
{}

void SplitTrackNtuplizer::resetTree() {
  //
  // Event
  //
  run = -999;
  subrun = -999;
  eventid = -999;
  passSelII = -999;
  //
  // Original track
  //  
  tk_id      = -999;
  tk_contain = -999;
  tk_nhits   = -999;
  tk_nvhits  = -999;
  tk_nhits_u = -999;
  tk_nhits_v = -999;
  tk_nhits_y = -999;
  tk_ndof    = -999;
  tk_chi2    = -999;
  tk_length  = -999;
  tk_drvtx    = -999;
  tk_mommumcs = -999;
  tk_mompmcs  = -999;
  tk_mommumcserr = -999;
  tk_mompmcserr  = -999;
  tk_mommumcsdll = -999;
  tk_mompmcsdll  = -999;
  tk_mommurng = -999;
  tk_momprng  = -999;
  tk_vtx_gpar.clear();  
  tk_vtx_lpar.clear();
  tk_vtx_gcov.clear();
  tk_vtx_lcov.clear();
  tk_end_gpar.clear();
  tk_end_lpar.clear();
  tk_end_gcov.clear();
  tk_end_lcov.clear();
  tk_mid_gpar.clear();
  // tk_mid_lpar.clear();
  // tk_mid_gcov.clear();
  // tk_mid_lcov.clear();
  tk_vtxmc_gpar.clear();
  tk_vtxmc_lpar.clear();
  tk_vtxmc_gcov.clear();
  tk_vtxmc_lcov.clear();
  tk_endmc_gpar.clear();
  tk_endmc_lpar.clear();
  tk_endmc_gcov.clear();
  tk_endmc_lcov.clear();
  tk_midmc_gpar.clear();
  // tk_midmc_lpar.clear();
  // tk_midmc_gcov.clear();
  // tk_midmc_lcov.clear();
  tk_hits_x.clear();
  tk_hits_rms.clear();
  tk_hits_plane.clear();
  tk_hits_wire.clear();
  vx_pos.clear();
  vx_ntks = -999;
  vx_tkid = -999;
  //
  // MC track
  //
  mc_pid     = -999;
  mc_mpid    = -999;
  mc_contain = -999;
  mc_length  = -999;
  mc_mom     = -999;
  mc_vtx_gpar  .clear();
  mc_end_gpar  .clear();
  mc_vtxtk_gpar.clear();
  mc_endtk_gpar.clear();
  mc_midtk_gpar.clear();
  mc_vtxtk1_gpar.clear();
  mc_endtk1_gpar.clear();
  mc_vtxtk2_gpar.clear();
  mc_endtk2_gpar.clear();
  //
  // Split track 1
  //  
  tk1_id      = -999;
  tk1_nhits   = -999;
  tk1_nvhits  = -999;
  tk1_nhits_u = -999;
  tk1_nhits_v = -999;
  tk1_nhits_y = -999;
  tk1_ndof    = -999;
  tk1_chi2    = -999;
  tk1_length  = -999;
  tk1_vtx_gpar.clear();  
  tk1_vtx_lpar.clear();
  tk1_vtx_gcov.clear();
  tk1_vtx_lcov.clear();
  tk1_end_gpar.clear();
  tk1_end_lpar.clear();
  tk1_end_gcov.clear();
  tk1_end_lcov.clear();
  tk1_vtxmc_gpar.clear();
  tk1_vtxmc_lpar.clear();
  tk1_vtxmc_gcov.clear();
  tk1_vtxmc_lcov.clear();
  tk1_midmc_gpar.clear();
  tk1_endmc_gpar.clear();
  tk1_endmc_lpar.clear();
  tk1_endmc_gcov.clear();
  tk1_endmc_lcov.clear();
  tk1_vtxtk_gpar.clear();
  tk1_vtxtk_lpar.clear();
  tk1_vtxtk_gcov.clear();
  tk1_vtxtk_lcov.clear();
  tk1_midtk_gpar.clear();
  tk1_midtk_lpar.clear();
  tk1_midtk_gcov.clear();
  tk1_midtk_lcov.clear();
  tk1_endtk_gpar.clear();
  tk1_endtk_lpar.clear();
  tk1_endtk_gcov.clear();
  tk1_endtk_lcov.clear();
  tk1_hits_x.clear();
  tk1_hits_rms.clear();
  tk1_hits_plane.clear();
  tk1_hits_wire.clear();
  vx1_pos.clear();
  vx1_ntks = -999;
  vx1_tkid = -999;
  vx1_dist = -999;
  //
  // Split track 2
  //  
  tk2_id      = -999;
  tk2_nhits   = -999;
  tk2_nvhits  = -999;
  tk2_nhits_u = -999;
  tk2_nhits_v = -999;
  tk2_nhits_y = -999;
  tk2_ndof    = -999;
  tk2_chi2    = -999;
  tk2_length  = -999;
  tk2_vtx_gpar.clear();  
  tk2_vtx_lpar.clear();
  tk2_vtx_gcov.clear();
  tk2_vtx_lcov.clear();
  tk2_end_gpar.clear();
  tk2_end_lpar.clear();
  tk2_end_gcov.clear();
  tk2_end_lcov.clear();
  tk2_vtxmc_gpar.clear();
  tk2_vtxmc_lpar.clear();
  tk2_vtxmc_gcov.clear();
  tk2_vtxmc_lcov.clear();
  tk2_midmc_gpar.clear();
  tk2_endmc_gpar.clear();
  tk2_endmc_lpar.clear();
  tk2_endmc_gcov.clear();
  tk2_endmc_lcov.clear();
  tk2_vtxtk_gpar.clear();
  tk2_vtxtk_lpar.clear();
  tk2_vtxtk_gcov.clear();
  tk2_vtxtk_lcov.clear();
  tk2_midtk_gpar.clear();
  tk2_midtk_lpar.clear();
  tk2_midtk_gcov.clear();
  tk2_midtk_lcov.clear();
  tk2_endtk_gpar.clear();
  tk2_endtk_lpar.clear();
  tk2_endtk_gcov.clear();
  tk2_endtk_lcov.clear();
  tk2_vtxtk1_gpar.clear();
  tk2_vtxtk1_lpar.clear();
  tk2_vtxtk1_gcov.clear();
  tk2_vtxtk1_lcov.clear();
  tk2_endtk1_gpar.clear();
  tk2_endtk1_lpar.clear();
  tk2_endtk1_gcov.clear();
  tk2_endtk1_lcov.clear();
  tk2_hits_x.clear();
  tk2_hits_rms.clear();
  tk2_hits_plane.clear();
  tk2_hits_wire.clear();
  vx2_pos.clear();
  vx2_ntks = -999;
  vx2_tkid = -999;
  vx2_dist = -999;
}

void SplitTrackNtuplizer::analyze(art::Event const & e)
{
  //
  using namespace std;
  using namespace trkf;
  using namespace recob::tracking;
  //
  // art::ValidHandle<art::TriggerResults> filter = e.getValidHandle<art::TriggerResults>("TriggerResults");
  // size_t ntp =  art::ServiceHandle<art::TriggerNamesService>()->size();
  // size_t ftp = ntp;
  // for (size_t itp=0;itp<ntp;itp++) {
  //   //std::cout << art::ServiceHandle<art::TriggerNamesService>()->getTrigPath(itp) << " " << filter->at(itp).accept()  << std::endl;
  //   if (art::ServiceHandle<art::TriggerNamesService>()->getTrigPath(itp)=="filtpro") ftp = itp; 
  // }
  // assert(ftp<ntp);
  //
  detinfo::DetectorProperties const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
  detinfo::DetectorClocks const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  TrackStatePropagator prop(1.0,0.1,10,10.,0.01,false);
  //
  art::InputTag PFInputTag(inputPFLabel);
  art::InputTag TrackInputTag(inputTracksLabel);
  art::ValidHandle<std::vector<recob::PFParticle> > inputPFParticle = e.getValidHandle<std::vector<recob::PFParticle> >(PFInputTag);
  auto assocTracks = std::unique_ptr<art::FindManyP<recob::Track> >(new art::FindManyP<recob::Track>(inputPFParticle, e, TrackInputTag));
  auto assocVertices = std::unique_ptr<art::FindManyP<recob::Vertex> >(new art::FindManyP<recob::Vertex>(inputPFParticle, e, PFInputTag));
  //
  art::InputTag TrackInputTag1st(inputTracksLabel1st);
  art::ValidHandle<std::vector<recob::Track> > Tracks1st = e.getValidHandle<std::vector<recob::Track> >(TrackInputTag1st);
  auto const& tkHitsAssn1st = *e.getValidHandle<art::Assns<recob::Track, recob::Hit> >(TrackInputTag1st);
  art::Handle<std::vector<recob::Vertex> > Vertices1st;
  std::unique_ptr<art::FindManyP<recob::Track> > assocVtxTracks1st;
  if (doSplitVertices) {
    e.getByLabel<std::vector<recob::Vertex> >(TrackInputTag1st, Vertices1st);
    assocVtxTracks1st = std::unique_ptr<art::FindManyP<recob::Track> >(new art::FindManyP<recob::Track>(Vertices1st, e, TrackInputTag1st));
  }
  //
  art::InputTag TrackInputTag2nd(inputTracksLabel2nd);
  art::ValidHandle<std::vector<recob::Track> > Tracks2nd = e.getValidHandle<std::vector<recob::Track> >(TrackInputTag2nd);
  auto const& tkHitsAssn2nd = *e.getValidHandle<art::Assns<recob::Track, recob::Hit> >(TrackInputTag2nd);
  art::Handle<std::vector<recob::Vertex> > Vertices2nd;
  std::unique_ptr<art::FindManyP<recob::Track> > assocVtxTracks2nd;
  if (doSplitVertices) {
    e.getByLabel<std::vector<recob::Vertex> >(TrackInputTag2nd, Vertices2nd);
    assocVtxTracks2nd = std::unique_ptr<art::FindManyP<recob::Track> >(new art::FindManyP<recob::Track>(Vertices2nd, e, TrackInputTag2nd));
  }
  //
  art::InputTag SimTrackInputTag("mcreco");
  const std::vector<sim::MCTrack>* simTracks = 0;
  if (e.isRealData()==0) {
    simTracks = e.getValidHandle<std::vector<sim::MCTrack> >(SimTrackInputTag).product();
  }
  //
  TrackMomentumCalculator tmc;
  //
  cout << inputTracksLabel << " " << inputTracksLabel1st << " " << inputTracksLabel2nd << endl;
  for (unsigned int iPF = 0; iPF < inputPFParticle->size(); ++iPF) {
    //
    const std::vector<art::Ptr<recob::Track> >& Tracks = assocTracks->at(iPF);
    auto const& tkHitsAssn = *e.getValidHandle<art::Assns<recob::Track, recob::Hit> >(TrackInputTag);
    //    
    cout << "pf #" << iPF << " Ntracks=" << Tracks.size() << " primary=" << inputPFParticle->at(iPF).IsPrimary() << " parent=" << inputPFParticle->at(iPF).Parent() << endl;
    //
    for (unsigned int iTrack = 0; iTrack < Tracks.size(); ++iTrack) {
      //
      resetTree();
      //
      run = e.run();
      subrun = e.subRun();
      eventid = e.event();
      // passSelII = filter->at(ftp).accept();
      //
      Point_t vertex;
      if (inputPFParticle->at(iPF).IsPrimary()==0) {
	auto parentPF = inputPFParticle->at(iPF).Parent();
	int nvtxtk = 0;
	for (auto idaug : inputPFParticle->at(parentPF).Daughters()) {
	  if (inputPFParticle->at(idaug).PdgCode()!=13) continue;
	  nvtxtk++;
	}
	if (nvtxtk>1) {
	  const std::vector<art::Ptr<recob::Vertex> >& parentVertices = assocVertices->at(parentPF);
	  for (unsigned int iVertex = 0; iVertex < parentVertices.size(); ++iVertex) {
	    art::Ptr<recob::Vertex> pvertex = parentVertices[iVertex];
	    double xyz[3];
	    pvertex->XYZ(xyz);
	    vertex = Point_t(xyz[0], xyz[1], xyz[2]);
	    break;
	  }
	  vx_pos.push_back(vertex.X());
	  vx_pos.push_back(vertex.Y());
	  vx_pos.push_back(vertex.Z());
	  vx_ntks = nvtxtk;
	  int dcount = 0;
	  for (auto idaug : inputPFParticle->at(parentPF).Daughters()) {
	    if (inputPFParticle->at(idaug).PdgCode()!=13) continue;
	    if (idaug == inputPFParticle->at(iPF).Self()) break;
	    dcount++;
	  }
	  vx_tkid = dcount;
	  //cout << "pf vertex " << vertex << " vx_ntks=" << vx_ntks << " vx_tkid=" << vx_tkid << endl;
	}
      }
      //
      art::Ptr<recob::Track> ptrack = Tracks[iTrack];
      std::vector<art::Ptr<recob::Hit> > inHits_tk;
      for (auto it = tkHitsAssn.begin(); it!=tkHitsAssn.end(); ++it) {
	if (it->first == ptrack) inHits_tk.push_back(it->second);
	else if (inHits_tk.size()>0) break;
      }
      //cout << inHits_tk.size() << endl;
      //
      auto id = ptrack->ID();
      //
      const recob::Track* ptrack1 = 0;
      std::vector<art::Ptr<recob::Hit> > inHits_tk1;
      for (unsigned int iTrack1 = 0; iTrack1 < Tracks1st->size(); ++iTrack1) {
	art::Ptr<recob::Track> ptrack1tmp(Tracks1st, iTrack1);
	if (ptrack1tmp->ID()!=id) continue;
	if (ptrack1tmp->CountValidPoints()<3) continue;
	ptrack1 = ptrack1tmp.get();
	for (auto it = tkHitsAssn1st.begin(); it!=tkHitsAssn1st.end(); ++it) {
	  if (it->first == ptrack1tmp) inHits_tk1.push_back(it->second);
	  else if (inHits_tk1.size()>0) break;
	}
	break;
      }
      if (ptrack1 && doSplitVertices) {
	Point_t vertex1;
	double mindist1 = 999999.;
	int ntrk = -1;
	int tkid = -1;
	for (unsigned int iVertex1 = 0; iVertex1 < Vertices1st->size(); ++iVertex1) {
	  art::Ptr<recob::Vertex> pvertex1(Vertices1st, iVertex1);
	  const std::vector<art::Ptr<recob::Track> >& vtxTracks1st = assocVtxTracks1st->at(iVertex1);
	  if (vtxTracks1st.size()<2) continue;
	  bool found = false;
	  int count = 0;
	  for (auto it : vtxTracks1st) {
	    if (it.get() == ptrack1) {
	      found = true;
	      break;
	    }
	    count++;
	  }
	  if (!found) continue;
	  ntrk = vtxTracks1st.size();
	  tkid = count;
	  double xyz[3];
	  pvertex1->XYZ(xyz);
	  Point_t tmpvertex1 = Point_t(xyz[0], xyz[1], xyz[2]);
	  double dist = (vertex-tmpvertex1).R();
	  if ( dist<mindist1 ) {
	    mindist1 = dist;
	    vertex1 = tmpvertex1;
	  }
	}
	if (ntrk>=0) {
	  vx1_pos.push_back(vertex1.X());
	  vx1_pos.push_back(vertex1.Y());
	  vx1_pos.push_back(vertex1.Z());
	  vx1_ntks = ntrk;
	  vx1_tkid = tkid;
	  vx1_dist = mindist1;
	  //cout << "vertex1=" << vertex1 << " track start=" << ptrack1->Trajectory().Start() << " vx1_ntks=" << vx1_ntks << " vx1_tkid=" << vx1_tkid << endl;
	}
      }
      //
      const recob::Track* ptrack2 = 0;
      std::vector<art::Ptr<recob::Hit> > inHits_tk2;
      for (unsigned int iTrack2 = 0; iTrack2 < Tracks2nd->size(); ++iTrack2) {
	art::Ptr<recob::Track> ptrack2tmp(Tracks2nd, iTrack2);
	if (ptrack2tmp->ID()!=id) continue;
	if (ptrack2tmp->CountValidPoints()<3) continue;
	ptrack2 = ptrack2tmp.get();
	for (auto it = tkHitsAssn2nd.begin(); it!=tkHitsAssn2nd.end(); ++it) {
	  if (it->first == ptrack2tmp) inHits_tk2.push_back(it->second);
	  else if (inHits_tk2.size()>0) break;
	}
	break;
      }
      if (ptrack2 && doSplitVertices) {
	Point_t vertex2;
	double mindist2 = 999999.;
	int ntrk = -1;
	int tkid = -1;
	for (unsigned int iVertex2 = 0; iVertex2 < Vertices2nd->size(); ++iVertex2) {
	  art::Ptr<recob::Vertex> pvertex2(Vertices2nd, iVertex2);
	  const std::vector<art::Ptr<recob::Track> >& vtxTracks2nd = assocVtxTracks2nd->at(iVertex2);
	  if (vtxTracks2nd.size()<2) continue;
	  bool found = false;
	  int count = 0;
	  for (auto it : vtxTracks2nd) {
	    if (it.get() == ptrack2) {
	      found = true;
	      break;
	    }
	    count++;
	  }
	  if (!found) continue;
	  ntrk = vtxTracks2nd.size();
	  tkid = count;
	  double xyz[3];
	  pvertex2->XYZ(xyz);
	  Point_t tmpvertex2 = Point_t(xyz[0], xyz[1], xyz[2]);
	  double dist = (vertex-tmpvertex2).R();
	  if ( dist<mindist2 ) {
	    mindist2 = dist;
	    vertex2 = tmpvertex2;
	  }
	}
	if (ntrk>=0) {
	  vx2_pos.push_back(vertex2.X());
	  vx2_pos.push_back(vertex2.Y());
	  vx2_pos.push_back(vertex2.Z());
	  vx2_ntks = ntrk;
	  vx2_tkid = tkid;
	  vx2_dist = mindist2;
	  //cout << "vertex2=" << vertex2 << " track start=" << ptrack2->Trajectory().Start() << " vx2_ntks=" << vx2_ntks << " vx2_tkid=" << vx2_tkid << endl;
	}
      }
      //
      if (ptrack1==0 || ptrack2==0) continue;
      //
      // Original track
      //
      int nhitstk_0 = 0;
      int nhitstk_1 = 0;
      int nhitstk_2 = 0;
      for (auto h : inHits_tk) {
	if (h->WireID().Plane==0) nhitstk_0++;
	if (h->WireID().Plane==1) nhitstk_1++;
	if (h->WireID().Plane==2) nhitstk_2++;
	double t = h->PeakTime();
	double trms = h->RMS();
	double x = theDetector->ConvertTicksToX(t, h->WireID().Plane, h->WireID().TPC, h->WireID().Cryostat);
	double xrms = trms * theDetector->GetXTicksCoefficient();
	tk_hits_x.push_back(x);
	tk_hits_rms.push_back(xrms);
	tk_hits_plane.push_back(h->WireID().Plane);
	tk_hits_wire.push_back(h->WireID().Wire);
      }
      //
      tk_id      = ptrack->ID();
      tk_nhits   = ptrack->NumberTrajectoryPoints();
      tk_nvhits  = ptrack->CountValidPoints();
      // std::cout << "ntuple full track start=" << ptrack->Trajectory().Start() << " dir" <<  ptrack->Trajectory().StartDirection() << " nh=" << ptrack->CountValidPoints() << " id=" << ptrack->ID() << std::endl;
      tk_nhits_u = nhitstk_0;
      tk_nhits_v = nhitstk_1;
      tk_nhits_y = nhitstk_2;
      tk_ndof    = ptrack->Ndof();
      tk_chi2    = ptrack->Chi2();
      tk_length  = ptrack->Length();
      //
      Point_t starttk(ptrack->Start().X(),ptrack->Start().Y(),ptrack->Start().Z());
      Point_t endtk(ptrack->End().X(),ptrack->End().Y(),ptrack->End().Z());
      Vector_t startdirtk(ptrack->StartDirection().X(),ptrack->StartDirection().Y(),ptrack->StartDirection().Z());    
      Vector_t enddirtk(ptrack->EndDirection().X(),ptrack->EndDirection().Y(),ptrack->EndDirection().Z());
      //
      tk_drvtx    = (vertex-starttk).R();
      auto mcsm = mcsFitMu.fitMcs(*ptrack);
      auto mcsp = mcsFitP.fitMcs(*ptrack);
      tk_mommumcs = mcsm.bestMomentum();
      tk_mompmcs  = mcsp.bestMomentum();
      tk_mommumcserr = mcsm.bestMomUncertainty();
      tk_mompmcserr  = mcsp.bestMomUncertainty();
      tk_mommumcsdll = mcsm.deltaLogLikelihood();
      tk_mompmcsdll  = mcsp.deltaLogLikelihood();
      tk_mommurng = tmc.GetTrackMomentum(tk_length,13);
      tk_momprng  = tmc.GetTrackMomentum(tk_length,2212);
      //
      std::cout << "starttk=" << starttk << " startdirtk=" << startdirtk << " endtk=" << endtk << " enddirtk=" << enddirtk << endl;
      //
      tk_contain = (starttk.X()>30.  && starttk.X()<230.  && endtk.X()>30.  && endtk.X()<230. &&
		    starttk.Y()>-85. && starttk.Y()<85.   && endtk.Y()>-85. && endtk.Y()<85.  &&
		    starttk.Z()>30.  && starttk.Z()<1010. && endtk.Z()>30.  && endtk.Z()<1010.);
      //
      auto vtx_gpar_tk = ptrack->VertexParametersGlobal6D();
      tk_vtx_gpar.assign( vtx_gpar_tk.begin(), vtx_gpar_tk.end() );
      auto vtx_gcov_tk = ptrack->VertexCovarianceGlobal6D();
      tk_vtx_gcov.assign( vtx_gcov_tk.begin(), vtx_gcov_tk.end() );
      auto vtx_lpar_tk = ptrack->VertexParametersLocal5D();
      tk_vtx_lpar.assign( vtx_lpar_tk.begin(), vtx_lpar_tk.end() );
      auto vtx_lcov_tk = ptrack->VertexCovarianceLocal5D();
      tk_vtx_lcov.assign( vtx_lcov_tk.begin(), vtx_lcov_tk.end() );
      Plane tk_startplane(starttk,startdirtk);
      TrackState startstate_tk(vtx_lpar_tk,vtx_lcov_tk,tk_startplane,true,13);
      //
      auto end_gpar_tk = ptrack->EndParametersGlobal6D();
      tk_end_gpar.assign( end_gpar_tk.begin(), end_gpar_tk.end() );
      auto end_gcov_tk = ptrack->EndCovarianceGlobal6D();
      tk_end_gcov.assign( end_gcov_tk.begin(), end_gcov_tk.end() );
      auto end_lpar_tk = ptrack->EndParametersLocal5D();
      tk_end_lpar.assign( end_lpar_tk.begin(), end_lpar_tk.end() );
      auto end_lcov_tk = ptrack->EndCovarianceLocal5D();
      tk_end_lcov.assign( end_lcov_tk.begin(), end_lcov_tk.end() );
      Plane tk_endplane(endtk,enddirtk);
      TrackState endstate_tk(end_lpar_tk,end_lcov_tk,tk_endplane,true,13);
      //
      Point_t midtk = ptrack->Trajectory().LocationAtPoint(ptrack->CountValidPoints()*0.5);
      Vector_t middirtk = ptrack->Trajectory().DirectionAtPoint(ptrack->CountValidPoints()*0.5);
      Plane tk_midplane(midtk,middirtk);
      tk_mid_gpar.push_back( midtk.X() );
      tk_mid_gpar.push_back( midtk.Y() );
      tk_mid_gpar.push_back( midtk.Z() );
      tk_mid_gpar.push_back( middirtk.X() );
      tk_mid_gpar.push_back( middirtk.Y() );
      tk_mid_gpar.push_back( middirtk.Z() );
      //
      // Split track 1
      //
      Point_t starttk1(ptrack1->Start().X(),ptrack1->Start().Y(),ptrack1->Start().Z());
      Point_t endtk1(ptrack1->End().X(),ptrack1->End().Y(),ptrack1->End().Z());
      Vector_t startdirtk1(ptrack1->StartDirection().X(),ptrack1->StartDirection().Y(),ptrack1->StartDirection().Z());
      Vector_t enddirtk1(ptrack1->EndDirection().X(),ptrack1->EndDirection().Y(),ptrack1->EndDirection().Z());
      std::cout << "starttk1=" << starttk1 << " startdirtk1=" << startdirtk1 << " endtk1=" << endtk1 << " enddirtk1=" << enddirtk1 << endl;
      int nhitstk1_0 = 0;
      int nhitstk1_1 = 0;
      int nhitstk1_2 = 0;
      for (auto h : inHits_tk1) {
	if (h->WireID().Plane==0) nhitstk1_0++;
	if (h->WireID().Plane==1) nhitstk1_1++;
	if (h->WireID().Plane==2) nhitstk1_2++;
	double t = h->PeakTime();
	double trms = h->RMS();
	double x = theDetector->ConvertTicksToX(t, h->WireID().Plane, h->WireID().TPC, h->WireID().Cryostat);
	double xrms = trms * theDetector->GetXTicksCoefficient();
	tk1_hits_x.push_back(x);
	tk1_hits_rms.push_back(xrms);
	tk1_hits_plane.push_back(h->WireID().Plane);
	tk1_hits_wire.push_back(h->WireID().Wire);
      }
      //
      tk1_id      = ptrack1->ID();
      tk1_nhits   = ptrack1->NumberTrajectoryPoints();
      tk1_nvhits  = ptrack1->CountValidPoints();
      // std::cout << "ntuple split track start=" << ptrack1->Trajectory().Start() << " dir" <<  ptrack1->Trajectory().StartDirection() << " nh=" << ptrack1->CountValidPoints() << " id=" << ptrack1->ID() << std::endl;
      tk1_nhits_u = nhitstk1_0;
      tk1_nhits_v = nhitstk1_1;
      tk1_nhits_y = nhitstk1_2;
      tk1_ndof    = ptrack1->Ndof();
      tk1_chi2    = ptrack1->Chi2();
      tk1_length  = ptrack1->Length();
      //
      auto vtx_gpar_tk1 = ptrack1->VertexParametersGlobal6D();
      tk1_vtx_gpar.assign( vtx_gpar_tk1.begin(), vtx_gpar_tk1.end() );
      auto vtx_gcov_tk1 = ptrack1->VertexCovarianceGlobal6D();
      tk1_vtx_gcov.assign( vtx_gcov_tk1.begin(), vtx_gcov_tk1.end() );
      auto vtx_lpar_tk1 = ptrack1->VertexParametersLocal5D();
      tk1_vtx_lpar.assign( vtx_lpar_tk1.begin(), vtx_lpar_tk1.end() );
      auto vtx_lcov_tk1 = ptrack1->VertexCovarianceLocal5D();
      tk1_vtx_lcov.assign( vtx_lcov_tk1.begin(), vtx_lcov_tk1.end() );
      TrackState startstate_tk1(vtx_lpar_tk1,vtx_lcov_tk1,Plane(starttk1,startdirtk1),true,13);
      //
      auto end_gpar_tk1 = ptrack1->EndParametersGlobal6D();
      tk1_end_gpar.assign( end_gpar_tk1.begin(), end_gpar_tk1.end() );
      auto end_gcov_tk1 = ptrack1->EndCovarianceGlobal6D();
      tk1_end_gcov.assign( end_gcov_tk1.begin(), end_gcov_tk1.end() );
      auto end_lpar_tk1 = ptrack1->EndParametersLocal5D();
      tk1_end_lpar.assign( end_lpar_tk1.begin(), end_lpar_tk1.end() );
      auto end_lcov_tk1 = ptrack1->EndCovarianceLocal5D();
      tk1_end_lcov.assign( end_lcov_tk1.begin(), end_lcov_tk1.end() );
      TrackState endstate_tk1(end_lpar_tk1,end_lcov_tk1,Plane(endtk1,enddirtk1),true,13);
      //
      //cout << "vtx_gpar_tk1=" << vtx_gpar_tk1 << " end_gpar_tk1=" << end_gpar_tk1 << endl;
      //
      bool propok_start_tk1 = true;
      auto dist_start_start_tk1 = prop.distanceToPlane(propok_start_tk1, startstate_tk1.position(), startstate_tk1.momentum().Unit(), tk_startplane);
      auto dist_start_end_tk1   = prop.distanceToPlane(propok_start_tk1, endstate_tk1.position(), endstate_tk1.momentum().Unit(), tk_startplane);
      bool flip_start_tk1 = fabs(dist_start_end_tk1)<fabs(dist_start_start_tk1);
      auto propState_start_tk1 = prop.propagateToPlane(propok_start_tk1, (flip_start_tk1 ? endstate_tk1 : startstate_tk1), tk_startplane, true, true, TrackStatePropagator::UNKNOWN);
      if (propok_start_tk1) {
	auto p5 = propState_start_tk1.parameters();
	auto c5 = propState_start_tk1.covariance();
	auto p6 = propState_start_tk1.parameters6D();
	auto c6 = propState_start_tk1.covariance6D();
	tk1_vtxtk_gpar.assign( p6.begin(), p6.end() );
	tk1_vtxtk_lpar.assign( p5.begin(), p5.end() );
	tk1_vtxtk_gcov.assign( c6.begin(), c6.end() );
	tk1_vtxtk_lcov.assign( c5.begin(), c5.end() );
      }    
      //
      bool propok_mid_tk1 = true;
      auto dist_mid_start_tk1 = prop.distanceToPlane(propok_mid_tk1, startstate_tk1.position(), startstate_tk1.momentum().Unit(), tk_midplane);
      auto dist_mid_end_tk1   = prop.distanceToPlane(propok_mid_tk1, endstate_tk1.position(), endstate_tk1.momentum().Unit(), tk_midplane);
      bool flip_mid_tk1 = fabs(dist_mid_end_tk1)<fabs(dist_mid_start_tk1);
      auto propState_mid_tk1 = prop.propagateToPlane(propok_mid_tk1, (flip_mid_tk1 ? endstate_tk1 : startstate_tk1), tk_midplane, true, true, TrackStatePropagator::UNKNOWN);
      if (propok_mid_tk1) {
	auto p5 = propState_mid_tk1.parameters();
	auto c5 = propState_mid_tk1.covariance();
	auto p6 = propState_mid_tk1.parameters6D();
	auto c6 = propState_mid_tk1.covariance6D();
	tk1_midtk_gpar.assign( p6.begin(), p6.end() );
	tk1_midtk_lpar.assign( p5.begin(), p5.end() );
	tk1_midtk_gcov.assign( c6.begin(), c6.end() );
	tk1_midtk_lcov.assign( c5.begin(), c5.end() );
      }    
      //
      bool propok_end_tk1 = true;
      auto dist_end_start_tk1 = prop.distanceToPlane(propok_end_tk1, startstate_tk1.position(), startstate_tk1.momentum().Unit(), tk_endplane);
      auto dist_end_end_tk1   = prop.distanceToPlane(propok_end_tk1, endstate_tk1.position(), endstate_tk1.momentum().Unit(), tk_endplane);
      bool flip_end_tk1 = fabs(dist_end_end_tk1)<fabs(dist_end_start_tk1);
      auto propState_end_tk1 = prop.propagateToPlane(propok_end_tk1, (flip_end_tk1 ? endstate_tk1 : startstate_tk1), tk_endplane, true, true, TrackStatePropagator::UNKNOWN);
      if (propok_end_tk1) {
	auto p5 = propState_end_tk1.parameters();
	auto c5 = propState_end_tk1.covariance();
	auto p6 = propState_end_tk1.parameters6D();
	auto c6 = propState_end_tk1.covariance6D();
	tk1_endtk_gpar.assign( p6.begin(), p6.end() );
	tk1_endtk_lpar.assign( p5.begin(), p5.end() );
	tk1_endtk_gcov.assign( c6.begin(), c6.end() );
	tk1_endtk_lcov.assign( c5.begin(), c5.end() );
      }    
      //
      // Split track 2
      //
      Point_t starttk2(ptrack2->Start().X(),ptrack2->Start().Y(),ptrack2->Start().Z());
      Point_t endtk2(ptrack2->End().X(),ptrack2->End().Y(),ptrack2->End().Z());
      Vector_t startdirtk2(ptrack2->StartDirection().X(),ptrack2->StartDirection().Y(),ptrack2->StartDirection().Z());
      Vector_t enddirtk2(ptrack2->EndDirection().X(),ptrack2->EndDirection().Y(),ptrack2->EndDirection().Z());
      std::cout << "starttk2=" << starttk2 << " startdirtk2=" << startdirtk2 << " endtk2=" << endtk2 << " enddirtk2=" << enddirtk2 << endl;
      int nhitstk2_0 = 0;
      int nhitstk2_1 = 0;
      int nhitstk2_2 = 0;
      for (auto h : inHits_tk2) {
	if (h->WireID().Plane==0) nhitstk2_0++;
	if (h->WireID().Plane==1) nhitstk2_1++;
	if (h->WireID().Plane==2) nhitstk2_2++;
	double t = h->PeakTime();
	double trms = h->RMS();
	double x = theDetector->ConvertTicksToX(t, h->WireID().Plane, h->WireID().TPC, h->WireID().Cryostat);
	double xrms = trms * theDetector->GetXTicksCoefficient();
	tk2_hits_x.push_back(x);
	tk2_hits_rms.push_back(xrms);
	tk2_hits_plane.push_back(h->WireID().Plane);
	tk2_hits_wire.push_back(h->WireID().Wire);
      }
      //
      tk2_id      = ptrack2->ID();
      tk2_nhits   = ptrack2->NumberTrajectoryPoints();
      tk2_nvhits  = ptrack2->CountValidPoints();
      tk2_nhits_u = nhitstk2_0;
      tk2_nhits_v = nhitstk2_1;
      tk2_nhits_y = nhitstk2_2;
      tk2_ndof    = ptrack2->Ndof();
      tk2_chi2    = ptrack2->Chi2();
      tk2_length  = ptrack2->Length();
      //
      auto vtx_gpar_tk2 = ptrack2->VertexParametersGlobal6D();
      tk2_vtx_gpar.assign( vtx_gpar_tk2.begin(), vtx_gpar_tk2.end() );
      auto vtx_gcov_tk2 = ptrack2->VertexCovarianceGlobal6D();
      tk2_vtx_gcov.assign( vtx_gcov_tk2.begin(), vtx_gcov_tk2.end() );
      auto vtx_lpar_tk2 = ptrack2->VertexParametersLocal5D();
      tk2_vtx_lpar.assign( vtx_lpar_tk2.begin(), vtx_lpar_tk2.end() );
      auto vtx_lcov_tk2 = ptrack2->VertexCovarianceLocal5D();
      tk2_vtx_lcov.assign( vtx_lcov_tk2.begin(), vtx_lcov_tk2.end() );
      TrackState startstate_tk2(vtx_lpar_tk2,vtx_lcov_tk2,Plane(starttk2,startdirtk2),true,13);
      //
      auto end_gpar_tk2 = ptrack2->EndParametersGlobal6D();
      tk2_end_gpar.assign( end_gpar_tk2.begin(), end_gpar_tk2.end() );
      auto end_gcov_tk2 = ptrack2->EndCovarianceGlobal6D();
      tk2_end_gcov.assign( end_gcov_tk2.begin(), end_gcov_tk2.end() );
      auto end_lpar_tk2 = ptrack2->EndParametersLocal5D();
      tk2_end_lpar.assign( end_lpar_tk2.begin(), end_lpar_tk2.end() );
      auto end_lcov_tk2 = ptrack2->EndCovarianceLocal5D();
      tk2_end_lcov.assign( end_lcov_tk2.begin(), end_lcov_tk2.end() );
      TrackState endstate_tk2(end_lpar_tk2,end_lcov_tk2,Plane(endtk2,enddirtk2),true,13);
      //
      bool swap = startdirtk2.Dot(startdirtk1)<0.;
      //
      //cout << "vtx_gpar_tk2=" << vtx_gpar_tk2 << " end_gpar_tk2=" << end_gpar_tk2 << endl;
      bool propok_start_tk1_tk2 = true;
      auto propState_start_tk1_tk2 = prop.propagateToPlane(propok_start_tk1_tk2, (swap ? endstate_tk2 : startstate_tk2), startstate_tk1.plane(), true, true, TrackStatePropagator::UNKNOWN);
      //cout << "propState_start_tk1_tk2 par6d=" << propState_start_tk1_tk2.parameters6D() << endl;
      if (propok_start_tk1_tk2) {
	auto p5 = propState_start_tk1_tk2.parameters();
	auto c5 = propState_start_tk1_tk2.covariance();
	auto p6 = propState_start_tk1_tk2.parameters6D();
	auto c6 = propState_start_tk1_tk2.covariance6D();
	tk2_vtxtk1_gpar.assign( p6.begin(), p6.end() );
	tk2_vtxtk1_lpar.assign( p5.begin(), p5.end() );
	tk2_vtxtk1_gcov.assign( c6.begin(), c6.end() );
	tk2_vtxtk1_lcov.assign( c5.begin(), c5.end() );
      }    
      //
      //cout << "end_gpar_tk2=" << end_gpar_tk2 << " end_gpar_tk2=" << end_gpar_tk2 << endl;
      bool propok_end_tk1_tk2 = true;
      auto propState_end_tk1_tk2 = prop.propagateToPlane(propok_end_tk1_tk2, (swap ? startstate_tk2 : endstate_tk2), endstate_tk1.plane(), true, true, TrackStatePropagator::UNKNOWN);
      //cout << "propState_end_tk1_tk2 par6d=" << propState_end_tk1_tk2.parameters6D() << endl;
      if (propok_end_tk1_tk2) {
	auto p5 = propState_end_tk1_tk2.parameters();
	auto c5 = propState_end_tk1_tk2.covariance();
	auto p6 = propState_end_tk1_tk2.parameters6D();
	auto c6 = propState_end_tk1_tk2.covariance6D();
	tk2_endtk1_gpar.assign( p6.begin(), p6.end() );
	tk2_endtk1_lpar.assign( p5.begin(), p5.end() );
	tk2_endtk1_gcov.assign( c6.begin(), c6.end() );
	tk2_endtk1_lcov.assign( c5.begin(), c5.end() );
      }    
      //
      bool propok_start_tk2 = true;
      auto dist_start_start_tk2 = prop.distanceToPlane(propok_start_tk2, startstate_tk2.position(), startstate_tk2.momentum().Unit(), tk_startplane);
      auto dist_start_end_tk2   = prop.distanceToPlane(propok_start_tk2, endstate_tk2.position(), endstate_tk2.momentum().Unit(), tk_startplane);
      bool flip_start_tk2 = fabs(dist_start_end_tk2)<fabs(dist_start_start_tk2);
      auto propState_start_tk2 = prop.propagateToPlane(propok_start_tk2, (flip_start_tk2 ? endstate_tk2 : startstate_tk2), tk_startplane, true, true, TrackStatePropagator::UNKNOWN);
      if (propok_start_tk2) {
	auto p5 = propState_start_tk2.parameters();
	auto c5 = propState_start_tk2.covariance();
	auto p6 = propState_start_tk2.parameters6D();
	auto c6 = propState_start_tk2.covariance6D();
	tk2_vtxtk_gpar.assign( p6.begin(), p6.end() );
	tk2_vtxtk_lpar.assign( p5.begin(), p5.end() );
	tk2_vtxtk_gcov.assign( c6.begin(), c6.end() );
	tk2_vtxtk_lcov.assign( c5.begin(), c5.end() );
      }    
      //
      bool propok_mid_tk2 = true;
      auto dist_mid_start_tk2 = prop.distanceToPlane(propok_mid_tk2, startstate_tk2.position(), startstate_tk2.momentum().Unit(), tk_midplane);
      auto dist_mid_end_tk2   = prop.distanceToPlane(propok_mid_tk2, endstate_tk2.position(), endstate_tk2.momentum().Unit(), tk_midplane);
      bool flip_mid_tk2 = fabs(dist_mid_end_tk2)<fabs(dist_mid_start_tk2);
      auto propState_mid_tk2 = prop.propagateToPlane(propok_mid_tk2, (flip_mid_tk2 ? endstate_tk2 : startstate_tk2), tk_midplane, true, true, TrackStatePropagator::UNKNOWN);
      if (propok_mid_tk2) {
	auto p5 = propState_mid_tk2.parameters();
	auto c5 = propState_mid_tk2.covariance();
	auto p6 = propState_mid_tk2.parameters6D();
	auto c6 = propState_mid_tk2.covariance6D();
	tk2_midtk_gpar.assign( p6.begin(), p6.end() );
	tk2_midtk_lpar.assign( p5.begin(), p5.end() );
	tk2_midtk_gcov.assign( c6.begin(), c6.end() );
	tk2_midtk_lcov.assign( c5.begin(), c5.end() );
      }    
      //
      bool propok_end_tk2 = true;
      auto dist_end_start_tk2 = prop.distanceToPlane(propok_end_tk2, startstate_tk2.position(), startstate_tk2.momentum().Unit(), tk_endplane);
      auto dist_end_end_tk2   = prop.distanceToPlane(propok_end_tk2, endstate_tk2.position(), endstate_tk2.momentum().Unit(), tk_endplane);
      bool flip_end_tk2 = fabs(dist_end_end_tk2)<fabs(dist_end_start_tk2);
      auto propState_end_tk2 = prop.propagateToPlane(propok_end_tk2, (flip_end_tk2 ? endstate_tk2 : startstate_tk2), tk_endplane, true, true, TrackStatePropagator::UNKNOWN);
      if (propok_end_tk2) {
	auto p5 = propState_end_tk2.parameters();
	auto c5 = propState_end_tk2.covariance();
	auto p6 = propState_end_tk2.parameters6D();
	auto c6 = propState_end_tk2.covariance6D();
	tk2_endtk_gpar.assign( p6.begin(), p6.end() );
	tk2_endtk_lpar.assign( p5.begin(), p5.end() );
	tk2_endtk_gcov.assign( c6.begin(), c6.end() );
	tk2_endtk_lcov.assign( c5.begin(), c5.end() );
      }    
      //
      // MC track
      //
      const sim::MCTrack* mctk = 0;
      Point_t mcstart, mcend;
      Vector_t mcstartmom, mcendmom;
      Vector_t mcstartdir, mcenddir;
      if (e.isRealData()==0) {
	for (unsigned int iMC = 0; iMC < (simTracks)->size(); ++iMC) {
	  const sim::MCTrack& mctrack = (simTracks)->at(iMC);
	  //
	  mcstartmom = Vector_t(mctrack.Start().Momentum().X()*0.001,mctrack.Start().Momentum().Y()*0.001,mctrack.Start().Momentum().Z()*0.001);
	  mcstartdir = mcstartmom.Unit();
	  mcendmom = Vector_t(mctrack.End().Momentum().X()*0.001,mctrack.End().Momentum().Y()*0.001,mctrack.End().Momentum().Z()*0.001);
	  mcenddir = mcendmom.Unit();
	  double dotvtx = startdirtk.X()*mcstartdir.X()+startdirtk.Y()*mcstartdir.Y()+startdirtk.Z()*mcstartdir.Z();
	  double dotend = enddirtk.X()*mcstartdir.X()+enddirtk.Y()*mcstartdir.Y()+enddirtk.Z()*mcstartdir.Z();
	  //
	  double g4Ticks = detClocks->TPCG4Time2Tick(mctrack.Start().T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
	  auto scecorr = SCE->GetPosOffsets( {mctrack.Start().Position().X(),mctrack.Start().Position().Y(),mctrack.Start().Position().Z()} );
	  double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
	  double yOffset = scecorr.Y();
	  double zOffset = scecorr.Z();
	  mcstart = Point_t(mctrack.Start().Position().X()+xOffset,mctrack.Start().Position().Y()+yOffset,mctrack.Start().Position().Z()+zOffset);
	  scecorr = SCE->GetPosOffsets( {mctrack.End().Position().X(),mctrack.End().Position().Y(),mctrack.End().Position().Z()} );
	  xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
	  yOffset = scecorr.Y();
	  zOffset = scecorr.Z();
	  mcend = Point_t(mctrack.End().Position().X()+xOffset,mctrack.End().Position().Y()+yOffset,mctrack.End().Position().Z()+zOffset);
	  //
	  bool match = false;
	  if ( (mcstart-starttk).R()<5 && dotvtx>0.9 )  match = true;
	  if ( (mcstart-endtk).R()<5   && dotend<-0.9 ) match = true;
	  if (!match) continue;
	  mctk = &mctrack;
	  break;
	}
      }
      //
      // cout << endl;
      // cout << "starttk=" << starttk << " startdirtk=" << startdirtk << " midtk=" << midtk << " middirtk=" << middirtk << " endtk=" << endtk << " enddirtk=" << enddirtk <<endl;
      // cout << "starttk1=" << starttk1 << " startdirtk1=" << startdirtk1 << " endtk1=" << endtk1 << " enddirtk1=" << enddirtk1 <<endl;
      // cout << "starttk2=" << starttk2 << " startdirtk2=" << startdirtk2 << " endtk2=" << endtk2 << " enddirtk2=" << enddirtk2 <<endl;
      // if (mctk) cout << "mcstart=" << mcstart << " mcstartdir=" << mcstartdir << " mcend=" << mcend << " mcenddir=" << mcenddir <<endl;
      //
      if (mctk) {
	//cout << "mcstart=" << mcstart << " mcstartdir=" << mcstartdir << " mcend=" << mcend << " mcenddir=" << mcenddir <<endl;
	//
	double g4Ticks = detClocks->TPCG4Time2Tick(mctk->Start().T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
	//double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0);
	//
	double mclen = 0.;
	//
	int istep_start = -1; double istep_start_dist = 9999999.;
	int istep_end   = -1; double istep_end_dist   = 9999999.;
	int istep_mid   = -1; double istep_mid_dist   = 9999999.;
	int istep_start1 = -1; double istep_start1_dist = 9999999.;
	int istep_start2 = -1; double istep_start2_dist = 9999999.;
	int istep_end1   = -1; double istep_end1_dist   = 9999999.;
	int istep_end2   = -1; double istep_end2_dist   = 9999999.;
	//
	for (unsigned int imc=0; imc<mctk->size(); ++imc) {
	  auto& mctrack = *mctk;
	  auto scecorr = SCE->GetPosOffsets( {mctrack[imc].X(),mctrack[imc].Y(),mctrack[imc].Z()} );
	  double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
	  double yOffset = scecorr.Y();
	  double zOffset = scecorr.Z();
	  Point_t mcsteppos(mctrack[imc].X()+xOffset,mctrack[imc].Y()+yOffset,mctrack[imc].Z()+zOffset);
	  //
	  double start_dist = (starttk-mcsteppos).R(); if (start_dist<istep_start_dist) { istep_start = imc; istep_start_dist = start_dist; }
	  double end_dist = (endtk-mcsteppos).R(); if (end_dist<istep_end_dist) { istep_end = imc; istep_end_dist = end_dist; }
	  double mid_dist = (midtk-mcsteppos).R(); if (mid_dist<istep_mid_dist) { istep_mid = imc; istep_mid_dist = mid_dist; }
	  double start1_dist = (starttk1-mcsteppos).R(); if (start1_dist<istep_start1_dist) { istep_start1 = imc; istep_start1_dist = start1_dist; }
	  double start2_dist = (starttk2-mcsteppos).R(); if (start2_dist<istep_start2_dist) { istep_start2 = imc; istep_start2_dist = start2_dist; }
	  double end1_dist = (endtk1-mcsteppos).R(); if (end1_dist<istep_end1_dist) { istep_end1 = imc; istep_end1_dist = end1_dist; }
	  double end2_dist = (endtk2-mcsteppos).R(); if (end2_dist<istep_end2_dist) { istep_end2 = imc; istep_end2_dist = end2_dist; }
	  //
	  if (imc<(mctk->size()-1)) {
	    mclen+=sqrt( (mctrack[imc+1].X()-mctrack[imc].X())*(mctrack[imc+1].X()-mctrack[imc].X()) +
			 (mctrack[imc+1].Y()-mctrack[imc].Y())*(mctrack[imc+1].Y()-mctrack[imc].Y()) +
			 (mctrack[imc+1].Z()-mctrack[imc].Z())*(mctrack[imc+1].Z()-mctrack[imc].Z()) );
	  }
	}
	//
	mc_pid     = mctk->PdgCode();
	mc_mpid    = mctk->MotherPdgCode();
	mc_contain = (mcstart.X()>30.  && mcstart.X()<230.  && mcend.X()>30.  && mcend.X()<230. &&
		      mcstart.Y()>-85. && mcstart.Y()<85.   && mcend.Y()>-85. && mcend.Y()<85.  &&
		      mcstart.Z()>30.  && mcstart.Z()<1010. && mcend.Z()>30.  && mcend.Z()<1010.);
	mc_length  = mclen;
	mc_mom     = mcstartmom.R();
	mc_vtx_gpar  .push_back(mcstart.X());
	mc_vtx_gpar  .push_back(mcstart.Y());
	mc_vtx_gpar  .push_back(mcstart.Z());
	mc_vtx_gpar  .push_back(mcstartmom.X());
	mc_vtx_gpar  .push_back(mcstartmom.Y());
	mc_vtx_gpar  .push_back(mcstartmom.Z());
	mc_end_gpar  .push_back(mcend.X());
	mc_end_gpar  .push_back(mcend.Y());
	mc_end_gpar  .push_back(mcend.Z());
	mc_end_gpar  .push_back(mcendmom.X());
	mc_end_gpar  .push_back(mcendmom.Y());
	mc_end_gpar  .push_back(mcendmom.Z());
	//
	if (istep_start>=0) {
	  auto scecorr = SCE->GetPosOffsets( {mctk->at(istep_start).X(),mctk->at(istep_start).Y(),mctk->at(istep_start).Z()} );
	  double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
	  double yOffset = scecorr.Y();
	  double zOffset = scecorr.Z();
	  Point_t mcstarttkpos(mctk->at(istep_start).X()+xOffset,mctk->at(istep_start).Y()+yOffset,mctk->at(istep_start).Z()+zOffset);
	  Vector_t mcstarttkmom(mctk->at(istep_start).Px(),mctk->at(istep_start).Py(),mctk->at(istep_start).Pz());
	  Vector_t mcstarttkdir = mcstarttkmom.Unit();
	  //cout << "istep_start=" << istep_start << " istep_start_dist=" << istep_start_dist << " mcstarttkpos=" << mcstarttkpos << " mcstarttkdir=" << mcstarttkdir << endl;
	  mc_vtxtk_gpar.push_back(mcstarttkpos.X());
	  mc_vtxtk_gpar.push_back(mcstarttkpos.Y());
	  mc_vtxtk_gpar.push_back(mcstarttkpos.Z());
	  mc_vtxtk_gpar.push_back(mcstarttkmom.X());
	  mc_vtxtk_gpar.push_back(mcstarttkmom.Y());
	  mc_vtxtk_gpar.push_back(mcstarttkmom.Z());
	  Plane mcstarttkplane(mcstarttkpos,mcstarttkdir);
	  bool propok = true;
	  auto propState = prop.propagateToPlane(propok, startstate_tk, mcstarttkplane, true, true, TrackStatePropagator::UNKNOWN);
	  if (propok) {
	    auto p5 = propState.parameters();
	    auto c5 = propState.covariance();
	    auto p6 = propState.parameters6D();
	    auto c6 = propState.covariance6D();
	    tk_vtxmc_gpar.assign( p6.begin(), p6.end() );
	    tk_vtxmc_lpar.assign( p5.begin(), p5.end() );
	    tk_vtxmc_gcov.assign( c6.begin(), c6.end() );
	    tk_vtxmc_lcov.assign( c5.begin(), c5.end() );
	  }
	}
	if (istep_mid>=0) {
	  auto scecorr = SCE->GetPosOffsets( {mctk->at(istep_mid).X(),mctk->at(istep_mid).Y(),mctk->at(istep_mid).Z()} );
	  double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
	  double yOffset = scecorr.Y();
	  double zOffset = scecorr.Z();
	  Point_t mcmidtkpos(mctk->at(istep_mid).X()+xOffset,mctk->at(istep_mid).Y()+yOffset,mctk->at(istep_mid).Z()+zOffset);
	  Vector_t mcmidtkmom(mctk->at(istep_mid).Px(),mctk->at(istep_mid).Py(),mctk->at(istep_mid).Pz());
	  Vector_t mcmidtkdir = mcmidtkmom.Unit();
	  //cout << "istep_mid=" << istep_mid << " istep_mid_dist=" << istep_mid_dist << " mcmidtkpos=" << mcmidtkpos << " mcmidtkdir=" << mcmidtkdir << endl;
	  mc_midtk_gpar.push_back(mcmidtkpos.X());
	  mc_midtk_gpar.push_back(mcmidtkpos.Y());
	  mc_midtk_gpar.push_back(mcmidtkpos.Z());
	  mc_midtk_gpar.push_back(mcmidtkmom.X());
	  mc_midtk_gpar.push_back(mcmidtkmom.Y());
	  mc_midtk_gpar.push_back(mcmidtkmom.Z());
	  Plane mcmidtkplane(mcmidtkpos,mcmidtkdir);
	  bool propok = true;
	  double d = prop.distanceToPlane (propok, midtk, middirtk, mcmidtkplane);
	  auto pos = prop.propagatedPosByDistance (midtk, middirtk, d);
	  if (propok) {
	    tk_midmc_gpar.push_back( pos.X() );
	    tk_midmc_gpar.push_back( pos.Y() );
	    tk_midmc_gpar.push_back( pos.Z() );
	    tk_midmc_gpar.push_back( middirtk.X() );
	    tk_midmc_gpar.push_back( middirtk.Y() );
	    tk_midmc_gpar.push_back( middirtk.Z() );
	  }
	  //
	  bool propok_midmc_tk1 = true;
	  auto dist_midmc_start_tk1 = prop.distanceToPlane(propok_midmc_tk1, startstate_tk1.position(), startstate_tk1.momentum().Unit(), mcmidtkplane);
	  auto dist_midmc_end_tk1   = prop.distanceToPlane(propok_midmc_tk1, endstate_tk1.position(), endstate_tk1.momentum().Unit(), mcmidtkplane);
	  bool flip_midmc_tk1 = fabs(dist_midmc_end_tk1)<fabs(dist_midmc_start_tk1);
	  auto pos_midmc_tk1 = (flip_midmc_tk1 ? prop.propagatedPosByDistance (endstate_tk1.position(), endstate_tk1.momentum().Unit(), dist_midmc_end_tk1)
				: prop.propagatedPosByDistance (startstate_tk1.position(), startstate_tk1.momentum().Unit(), dist_midmc_start_tk1) );
	  if (propok_midmc_tk1) {
	    tk1_midmc_gpar.push_back( pos_midmc_tk1.X() );
	    tk1_midmc_gpar.push_back( pos_midmc_tk1.Y() );
	    tk1_midmc_gpar.push_back( pos_midmc_tk1.Z() );
	    tk1_midmc_gpar.push_back( (flip_midmc_tk1 ? endstate_tk1.momentum().Unit() : startstate_tk1.momentum().Unit()).X() );
	    tk1_midmc_gpar.push_back( (flip_midmc_tk1 ? endstate_tk1.momentum().Unit() : startstate_tk1.momentum().Unit()).Y() );
	    tk1_midmc_gpar.push_back( (flip_midmc_tk1 ? endstate_tk1.momentum().Unit() : startstate_tk1.momentum().Unit()).Z() );
	  }
	  //
	  bool propok_midmc_tk2 = true;
	  auto dist_midmc_start_tk2 = prop.distanceToPlane(propok_midmc_tk2, startstate_tk2.position(), startstate_tk2.momentum().Unit(), mcmidtkplane);
	  auto dist_midmc_end_tk2   = prop.distanceToPlane(propok_midmc_tk2, endstate_tk2.position(), endstate_tk2.momentum().Unit(), mcmidtkplane);
	  bool flip_midmc_tk2 = fabs(dist_midmc_end_tk2)<fabs(dist_midmc_start_tk2);
	  auto pos_midmc_tk2 = (flip_midmc_tk2 ? prop.propagatedPosByDistance (endstate_tk2.position(), endstate_tk2.momentum().Unit(), dist_midmc_end_tk2)
				: prop.propagatedPosByDistance (startstate_tk2.position(), startstate_tk2.momentum().Unit(), dist_midmc_start_tk2) );
	  if (propok_midmc_tk2) {
	    tk2_midmc_gpar.push_back( pos_midmc_tk2.X() );
	    tk2_midmc_gpar.push_back( pos_midmc_tk2.Y() );
	    tk2_midmc_gpar.push_back( pos_midmc_tk2.Z() );
	    tk2_midmc_gpar.push_back( (flip_midmc_tk2 ? endstate_tk2.momentum().Unit() : startstate_tk2.momentum().Unit()).X() );
	    tk2_midmc_gpar.push_back( (flip_midmc_tk2 ? endstate_tk2.momentum().Unit() : startstate_tk2.momentum().Unit()).Y() );
	    tk2_midmc_gpar.push_back( (flip_midmc_tk2 ? endstate_tk2.momentum().Unit() : startstate_tk2.momentum().Unit()).Z() );
	  }
	}
	if (istep_end>=0) {
	  auto scecorr = SCE->GetPosOffsets( {mctk->at(istep_end).X(),mctk->at(istep_end).Y(),mctk->at(istep_end).Z()} );
	  double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
	  double yOffset = scecorr.Y();
	  double zOffset = scecorr.Z();
	  Point_t mcendtkpos(mctk->at(istep_end).X()+xOffset,mctk->at(istep_end).Y()+yOffset,mctk->at(istep_end).Z()+zOffset);
	  Vector_t mcendtkmom(mctk->at(istep_end).Px(),mctk->at(istep_end).Py(),mctk->at(istep_end).Pz());
	  Vector_t mcendtkdir = mcendtkmom.Unit();
	  //cout << "istep_end=" << istep_end << " istep_end_dist=" << istep_end_dist << " mcendtkpos=" << mcendtkpos << " mcendtkdir=" << mcendtkdir << endl;
	  mc_endtk_gpar.push_back(mcendtkpos.X());
	  mc_endtk_gpar.push_back(mcendtkpos.Y());
	  mc_endtk_gpar.push_back(mcendtkpos.Z());
	  mc_endtk_gpar.push_back(mcendtkmom.X());
	  mc_endtk_gpar.push_back(mcendtkmom.Y());
	  mc_endtk_gpar.push_back(mcendtkmom.Z());
	  Plane mcendtkplane(mcendtkpos,mcendtkdir);
	  bool propok = true;
	  auto propState = prop.propagateToPlane(propok, endstate_tk, mcendtkplane, true, true, TrackStatePropagator::UNKNOWN);
	  if (propok) {
	    auto p5 = propState.parameters();
	    auto c5 = propState.covariance();
	    auto p6 = propState.parameters6D();
	    auto c6 = propState.covariance6D();
	    tk_endmc_gpar.assign( p6.begin(), p6.end() );
	    tk_endmc_lpar.assign( p5.begin(), p5.end() );
	    tk_endmc_gcov.assign( c6.begin(), c6.end() );
	    tk_endmc_lcov.assign( c5.begin(), c5.end() );
	  }
	}
	if (istep_start1>=0) {
	  auto scecorr = SCE->GetPosOffsets( {mctk->at(istep_start1).X(),mctk->at(istep_start1).Y(),mctk->at(istep_start1).Z()} );
	  double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
	  double yOffset = scecorr.Y();
	  double zOffset = scecorr.Z();
	  Point_t mcstarttk1pos(mctk->at(istep_start1).X()+xOffset,mctk->at(istep_start1).Y()+yOffset,mctk->at(istep_start1).Z()+zOffset);
	  Vector_t mcstarttk1mom(mctk->at(istep_start1).Px(),mctk->at(istep_start1).Py(),mctk->at(istep_start1).Pz());
	  Vector_t mcstarttk1dir = mcstarttk1mom.Unit();
	  //cout << "istep_start1=" << istep_start1 << " istep_start1_dist=" << istep_start1_dist << " mcstarttk1pos=" << mcstarttk1pos << " mcstarttk1dir=" << mcstarttk1dir << endl;
	  mc_vtxtk1_gpar.push_back(mcstarttk1pos.X());
	  mc_vtxtk1_gpar.push_back(mcstarttk1pos.Y());
	  mc_vtxtk1_gpar.push_back(mcstarttk1pos.Z());
	  mc_vtxtk1_gpar.push_back(mcstarttk1mom.X());
	  mc_vtxtk1_gpar.push_back(mcstarttk1mom.Y());
	  mc_vtxtk1_gpar.push_back(mcstarttk1mom.Z());
	  Plane mcstarttk1plane(mcstarttk1pos,mcstarttk1dir);
	  bool propok = true;
	  //cout << "startstate_tk1 pos=" << startstate_tk1.position() << " dir=" << startstate_tk1.momentum().Unit() << endl;
	  auto propState = prop.propagateToPlane(propok, startstate_tk1, mcstarttk1plane, true, true, TrackStatePropagator::UNKNOWN);
	  if (propok) {
	    auto p5 = propState.parameters();
	    auto c5 = propState.covariance();
	    auto p6 = propState.parameters6D();
	    auto c6 = propState.covariance6D();
	    //cout << "prop pos=" << propState.position() << " mom=" << propState.momentum()<< " dir=" << propState.momentum().Unit() << endl;
	    tk1_vtxmc_gpar.assign( p6.begin(), p6.end() );
	    tk1_vtxmc_lpar.assign( p5.begin(), p5.end() );
	    tk1_vtxmc_gcov.assign( c6.begin(), c6.end() );
	    tk1_vtxmc_lcov.assign( c5.begin(), c5.end() );
	  }
	}
	if (istep_end1>=0) {
	  auto scecorr = SCE->GetPosOffsets( {mctk->at(istep_end1).X(),mctk->at(istep_end1).Y(),mctk->at(istep_end1).Z()} );
	  double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
	  double yOffset = scecorr.Y();
	  double zOffset = scecorr.Z();
	  Point_t mcendtk1pos(mctk->at(istep_end1).X()+xOffset,mctk->at(istep_end1).Y()+yOffset,mctk->at(istep_end1).Z()+zOffset);
	  Vector_t mcendtk1mom(mctk->at(istep_end1).Px(),mctk->at(istep_end1).Py(),mctk->at(istep_end1).Pz());
	  Vector_t mcendtk1dir = mcendtk1mom.Unit();
	  //cout << "istep_end1=" << istep_end1 << " istep_end1_dist=" << istep_end1_dist << " mcendtk1pos=" << mcendtk1pos << " mcendtk1dir=" << mcendtk1dir << endl;
	  mc_endtk1_gpar.push_back(mcendtk1pos.X());
	  mc_endtk1_gpar.push_back(mcendtk1pos.Y());
	  mc_endtk1_gpar.push_back(mcendtk1pos.Z());
	  mc_endtk1_gpar.push_back(mcendtk1mom.X());
	  mc_endtk1_gpar.push_back(mcendtk1mom.Y());
	  mc_endtk1_gpar.push_back(mcendtk1mom.Z());
	  Plane mcendtk1plane(mcendtk1pos,mcendtk1dir);
	  bool propok = true;
	  auto propState = prop.propagateToPlane(propok, endstate_tk1, mcendtk1plane, true, true, TrackStatePropagator::UNKNOWN);
	  if (propok) {
	    auto p5 = propState.parameters();
	    auto c5 = propState.covariance();
	    auto p6 = propState.parameters6D();
	    auto c6 = propState.covariance6D();
	    //cout << "prop pos=" << propState.position() << " mom=" << propState.momentum()<< " dir=" << propState.momentum().Unit() << endl;
	    tk1_endmc_gpar.assign( p6.begin(), p6.end() );
	    tk1_endmc_lpar.assign( p5.begin(), p5.end() );
	    tk1_endmc_gcov.assign( c6.begin(), c6.end() );
	    tk1_endmc_lcov.assign( c5.begin(), c5.end() );
	  }
	}
	if (istep_start2>=0) {
	  auto scecorr = SCE->GetPosOffsets( {mctk->at(istep_start2).X(),mctk->at(istep_start2).Y(),mctk->at(istep_start2).Z()} );
	  double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
	  double yOffset = scecorr.Y();
	  double zOffset = scecorr.Z();
	  Point_t mcstarttk2pos(mctk->at(istep_start2).X()+xOffset,mctk->at(istep_start2).Y()+yOffset,mctk->at(istep_start2).Z()+zOffset);
	  Vector_t mcstarttk2mom(mctk->at(istep_start2).Px(),mctk->at(istep_start2).Py(),mctk->at(istep_start2).Pz());
	  Vector_t mcstarttk2dir = mcstarttk2mom.Unit();
	  //cout << "istep_start2=" << istep_start2 << " istep_start2_dist=" << istep_start2_dist << " mcstarttk2pos=" << mcstarttk2pos << " mcstarttk2dir=" << mcstarttk2dir << endl;
	  mc_vtxtk2_gpar.push_back(mcstarttk2pos.X());
	  mc_vtxtk2_gpar.push_back(mcstarttk2pos.Y());
	  mc_vtxtk2_gpar.push_back(mcstarttk2pos.Z());
	  mc_vtxtk2_gpar.push_back(mcstarttk2mom.X());
	  mc_vtxtk2_gpar.push_back(mcstarttk2mom.Y());
	  mc_vtxtk2_gpar.push_back(mcstarttk2mom.Z());
	  Plane mcstarttk2plane(mcstarttk2pos,mcstarttk2dir);
	  //cout << "startstate_tk2 pos=" << startstate_tk2.position() << " dir=" << startstate_tk2.momentum().Unit() << endl;
	  bool propok = true;
	  auto propState = prop.propagateToPlane(propok, startstate_tk2, mcstarttk2plane, true, true, TrackStatePropagator::UNKNOWN);
	  //auto propState = prop.propagateToPlane(propok, startstate_tk2, mcstarttk2plane, false, true, TrackStatePropagator::UNKNOWN);
	  if (propok) {
	    auto p5 = propState.parameters();
	    auto c5 = propState.covariance();
	    auto p6 = propState.parameters6D();
	    auto c6 = propState.covariance6D();
	    //cout << "prop pos=" << propState.position() << " mom=" << propState.momentum()<< " dir=" << propState.momentum().Unit() << endl;
	    tk2_vtxmc_gpar.assign( p6.begin(), p6.end() );
	    tk2_vtxmc_lpar.assign( p5.begin(), p5.end() );
	    tk2_vtxmc_gcov.assign( c6.begin(), c6.end() );
	    tk2_vtxmc_lcov.assign( c5.begin(), c5.end() );
	  }
	}
	if (istep_end2>=0) {
	  auto scecorr = SCE->GetPosOffsets( {mctk->at(istep_end2).X(),mctk->at(istep_end2).Y(),mctk->at(istep_end2).Z()} );
	  double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
	  double yOffset = scecorr.Y();
	  double zOffset = scecorr.Z();
	  Point_t mcendtk2pos(mctk->at(istep_end2).X()+xOffset,mctk->at(istep_end2).Y()+yOffset,mctk->at(istep_end2).Z()+zOffset);
	  Vector_t mcendtk2mom(mctk->at(istep_end2).Px(),mctk->at(istep_end2).Py(),mctk->at(istep_end2).Pz());
	  Vector_t mcendtk2dir = mcendtk2mom.Unit();
	  //cout << "istep_end2=" << istep_end2 << " istep_end2_dist=" << istep_end2_dist << " mcendtk2pos=" << mcendtk2pos << " mcendtk2dir=" << mcendtk2dir << endl;
	  mc_endtk2_gpar.push_back(mcendtk2pos.X());
	  mc_endtk2_gpar.push_back(mcendtk2pos.Y());
	  mc_endtk2_gpar.push_back(mcendtk2pos.Z());
	  mc_endtk2_gpar.push_back(mcendtk2mom.X());
	  mc_endtk2_gpar.push_back(mcendtk2mom.Y());
	  mc_endtk2_gpar.push_back(mcendtk2mom.Z());
	  Plane mcendtk2plane(mcendtk2pos,mcendtk2dir);
	  bool propok = true;
	  auto propState = prop.propagateToPlane(propok, endstate_tk2, mcendtk2plane, true, true, TrackStatePropagator::UNKNOWN);
	  if (propok) {
	    auto p5 = propState.parameters();
	    auto c5 = propState.covariance();
	    auto p6 = propState.parameters6D();
	    auto c6 = propState.covariance6D();
	    //cout << "prop pos=" << propState.position() << " mom=" << propState.momentum()<< " dir=" << propState.momentum().Unit() << endl;
	    tk2_endmc_gpar.assign( p6.begin(), p6.end() );
	    tk2_endmc_lpar.assign( p5.begin(), p5.end() );
	    tk2_endmc_gcov.assign( c6.begin(), c6.end() );
	    tk2_endmc_lcov.assign( c5.begin(), c5.end() );
	  }
	}
      }//if (mctk)
      //
      tree->Fill();
      //
    }
    //
  }
  //
}

DEFINE_ART_MODULE(SplitTrackNtuplizer)
