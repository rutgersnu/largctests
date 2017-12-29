////////////////////////////////////////////////////////////////////////
// Class:       TrackValidation
// Plugin Type: analyzer (art v2_05_00)
// File:        TrackValidation_module.cc
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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcore/Geometry/Geometry.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingPlane.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"

#include "lardata/RecoObjects/PropXYZPlane.h"
#include "lardata/RecoObjects/SurfXYZPlane.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardata/RecoObjects/TrackingPlaneHelper.h"
#include "lardata/RecoObjects/TrackState.h"

#include "TH1.h"
#include "TH2.h"

class TrackValidation;

class TrackValidation : public art::EDAnalyzer {
public:
  explicit TrackValidation(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackValidation(TrackValidation const &) = delete;
  TrackValidation(TrackValidation &&) = delete;
  TrackValidation & operator = (TrackValidation const &) = delete;
  TrackValidation & operator = (TrackValidation &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  void beginJob() override;
  
private:

  std::string inputTracksLabel;
  unsigned int minHits;
  //
  TH1F* all_NTracks;
  //
  TH1F* all_NumberTrajectoryPoints;
  TH1F* all_CountValidPoints;
  TH1F* all_HasMomentum;
  TH1F* all_Length;
  TH1F* all_Chi2;
  TH1F* all_Chi2PerNdof;
  TH1F* all_Ndof;
  TH1F* all_ParticleId;
  TH1F* all_Theta;
  TH1F* all_Phi;
  TH1F* all_ZenithAngle;
  TH1F* all_AzimuthAngle;
  TH1F* all_NumberCovariance;
  //
  TH1F* all_hitPlane0Frac;
  TH1F* all_hitPlane1Frac;
  TH1F* all_hitPlane2Frac;
  //
  TH2F* all_XYStart;
  TH2F* all_XYEnd;
  TH2F* all_ZRStart;
  TH2F* all_ZREnd;
  //
  TH1F* NumberTrajectoryPoints;
  TH1F* CountValidPoints;
  TH1F* HasMomentum;
  TH1F* Length;
  TH1F* dLength;
  TH1F* dLengthRel;
  TH1F* Chi2;
  TH1F* Chi2PerNdof;
  TH1F* Ndof;
  TH1F* ParticleId;
  TH1F* Theta;
  TH1F* Phi;
  TH1F* ZenithAngle;
  TH1F* AzimuthAngle;
  TH1F* NumberCovariance;
  //
  TH1F* dx_assoc;
  TH1F* dy_assoc;
  TH1F* dz_assoc;
  TH1F* dx_prop_assoc;
  TH1F* dy_prop_assoc;
  TH1F* dz_prop_assoc;
  TH1F* dx_prop_assoc_rerr;
  TH1F* dy_prop_assoc_rerr;
  TH1F* dz_prop_assoc_rerr;
  TH1F* dx_pull_prop_assoc;
  TH1F* dy_pull_prop_assoc;
  TH1F* dz_pull_prop_assoc;
  //
  TH1F* duxr_assoc;
  TH1F* duyr_assoc;
  TH1F* duzr_assoc;
  TH1F* duxr_prop_assoc;
  TH1F* duyr_prop_assoc;
  TH1F* duzr_prop_assoc;
  TH1F* dux_prop_assoc_rerr;
  TH1F* duy_prop_assoc_rerr;
  TH1F* duz_prop_assoc_rerr;
  TH1F* dux_pull_prop_assoc;
  TH1F* duy_pull_prop_assoc;
  TH1F* duz_pull_prop_assoc;
  //
  TH1F* dp0_prop_assoc;
  TH1F* dp1_prop_assoc;
  TH1F* dp2_prop_assoc;
  TH1F* dp3_prop_assoc;
  TH1F* dp0_prop_assoc_rerr;
  TH1F* dp1_prop_assoc_rerr;
  TH1F* dp2_prop_assoc_rerr;
  TH1F* dp3_prop_assoc_rerr;
  TH1F* dp0_pull_prop_assoc;
  TH1F* dp1_pull_prop_assoc;
  TH1F* dp2_pull_prop_assoc;
  TH1F* dp3_pull_prop_assoc;
  //
  TH2F* x_vs_prop_pullX;
  TH2F* y_vs_prop_pullX;
  TH2F* z_vs_prop_pullX;
  TH2F* dirx_vs_prop_pullX;
  TH2F* diry_vs_prop_pullX;
  TH2F* dirz_vs_prop_pullX;
  TH2F* p_vs_prop_pullX;
  TH2F* theta_vs_prop_pullX;
  TH2F* phi_vs_prop_pullX;
  TH2F* zenith_vs_prop_pullX;
  TH2F* azimuth_vs_prop_pullX;
  //
  TH1F* hit_track_res;
  TH1F* hit_track_pull;
  TH1F* hit_track_terr;
  TH1F* hit_track_herr;
  TH1F* hit_track_plane0_res;
  TH1F* hit_track_plane0_pull;
  TH1F* hit_track_plane0_terr;
  TH1F* hit_track_plane0_herr;
  TH1F* hit_track_plane1_res;
  TH1F* hit_track_plane1_pull;
  TH1F* hit_track_plane1_terr;
  TH1F* hit_track_plane1_herr;
  TH1F* hit_track_plane2_res;
  TH1F* hit_track_plane2_pull;
  TH1F* hit_track_plane2_terr;
  TH1F* hit_track_plane2_herr;
  TH2F* hit_track_plane0_res_vs_wire;
  TH2F* hit_track_plane1_res_vs_wire;
  TH2F* hit_track_plane2_res_vs_wire;
  TH2F* hit_track_plane0_res_vs_wire_test;
  TH2F* hit_track_plane1_res_vs_wire_test;
  TH2F* hit_track_plane2_res_vs_wire_test;
  //
};

void TrackValidation::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  //
  all_NTracks = tfs->make<TH1F>("all_NTracks","all_NTracks", 100, 0, 100);
  //
  all_NumberTrajectoryPoints = tfs->make<TH1F>("all_NumberTrajectoryPoints","all_NumberTrajectoryPoints", 100, 0, 2500);
  all_CountValidPoints       = tfs->make<TH1F>("all_CountValidPoints      ","all_CountValidPoints      ", 100, 0, 2500);
  all_HasMomentum            = tfs->make<TH1F>("all_HasMomentum           ","all_HasMomentum           ",   2, 0,    2);
  all_Length                 = tfs->make<TH1F>("all_Length                ","all_Length                ", 100,  0, 500);
  all_Chi2                   = tfs->make<TH1F>("all_Chi2                  ","all_Chi2                  ", 100,  0, 100);
  all_Chi2PerNdof            = tfs->make<TH1F>("all_Chi2PerNdof           ","all_Chi2PerNdof           ", 100,  0,  10);
  all_Ndof                   = tfs->make<TH1F>("all_Ndof                  ","all_Ndof                  ", 100, 0, 2500);
  all_ParticleId             = tfs->make<TH1F>("all_ParticleId            ","all_ParticleId            ", 100,-50,  50);
  all_Theta                  = tfs->make<TH1F>("all_Theta                 ","all_Theta                 ",  66,  0, 3.3);
  all_Phi                    = tfs->make<TH1F>("all_Phi                   ","all_Phi                   ",  66,-3.3,3.3);
  all_ZenithAngle            = tfs->make<TH1F>("all_ZenithAngle           ","all_ZenithAngle           ",  66,  0, 3.3);
  all_AzimuthAngle           = tfs->make<TH1F>("all_AzimuthAngle          ","all_AzimuthAngle          ",  66,-3.3,3.3);
  all_NumberCovariance       = tfs->make<TH1F>("all_NumberCovariance      ","all_NumberCovariance      ",  10, 0,   10);
  //
  all_hitPlane0Frac = tfs->make<TH1F>("all_hitPlane0Frac","all_hitPlane0Frac",101,0,1.01);
  all_hitPlane1Frac = tfs->make<TH1F>("all_hitPlane1Frac","all_hitPlane1Frac",101,0,1.01);
  all_hitPlane2Frac = tfs->make<TH1F>("all_hitPlane2Frac","all_hitPlane2Frac",101,0,1.01);
  //
  all_XYStart = tfs->make<TH2F>("all_XYStart","all_XYStart",100,-50,350,  100,-250,250);
  all_XYEnd   = tfs->make<TH2F>("all_XYEnd  ","all_XYEnd  ",100,-50,350,  100,-250,250);
  all_ZRStart = tfs->make<TH2F>("all_ZRStart","all_ZRStart",120,-50,1050, 100,   0,350);
  all_ZREnd   = tfs->make<TH2F>("all_ZREnd  ","all_ZREnd  ",120,-50,1050, 100,   0,350);
  //
  NumberTrajectoryPoints = tfs->make<TH1F>("NumberTrajectoryPoints","NumberTrajectoryPoints", 100, 0, 2500);
  CountValidPoints       = tfs->make<TH1F>("CountValidPoints      ","CountValidPoints      ", 100, 0, 2500);
  HasMomentum            = tfs->make<TH1F>("HasMomentum           ","HasMomentum           ",   2, 0,    2);
  Length                 = tfs->make<TH1F>("Length                ","Length                ", 100,  0, 500);
  dLength                = tfs->make<TH1F>("dLength               ","dLength               ", 100,-100,100);
  dLengthRel             = tfs->make<TH1F>("dLengthRel            ","dLengthRel            ", 100, -1., 1.);
  Chi2                   = tfs->make<TH1F>("Chi2                  ","Chi2                  ", 100,  0, 100);
  Chi2PerNdof            = tfs->make<TH1F>("Chi2PerNdof           ","Chi2PerNdof           ", 100,  0,  10);
  Ndof                   = tfs->make<TH1F>("Ndof                  ","Ndof                  ", 100, 0, 2500);
  ParticleId             = tfs->make<TH1F>("ParticleId            ","ParticleId            ", 100,-50,  50);
  Theta                  = tfs->make<TH1F>("Theta                 ","Theta                 ",  66,  0, 3.3);
  Phi                    = tfs->make<TH1F>("Phi                   ","Phi                   ",  66,-3.3,3.3);
  ZenithAngle            = tfs->make<TH1F>("ZenithAngle           ","ZenithAngle           ",  66,  0, 3.3);
  AzimuthAngle           = tfs->make<TH1F>("AzimuthAngle          ","AzimuthAngle          ",  66,-3.3,3.3);
  NumberCovariance       = tfs->make<TH1F>("NumberCovariance      ","NumberCovariance      ",  10, 0,   10);
  //
  dx_assoc           = tfs->make<TH1F>("dx_assoc          ","dx_assoc          ",100,-2.5,2.5);
  dy_assoc           = tfs->make<TH1F>("dy_assoc          ","dy_assoc          ",100,-2.5,2.5);
  dz_assoc           = tfs->make<TH1F>("dz_assoc          ","dz_assoc          ",100,-2.5,2.5);
  dx_prop_assoc      = tfs->make<TH1F>("dx_prop_assoc     ","dx_prop_assoc     ",100,-1.0,1.0);
  dy_prop_assoc      = tfs->make<TH1F>("dy_prop_assoc     ","dy_prop_assoc     ",100,-1.0,1.0);
  dz_prop_assoc      = tfs->make<TH1F>("dz_prop_assoc     ","dz_prop_assoc     ",100,-1.0,1.0);
  dx_prop_assoc_rerr = tfs->make<TH1F>("dx_prop_assoc_rerr","dx_prop_assoc_rerr",100, 0.,1.);
  dy_prop_assoc_rerr = tfs->make<TH1F>("dy_prop_assoc_rerr","dy_prop_assoc_rerr",100, 0.,1.);
  dz_prop_assoc_rerr = tfs->make<TH1F>("dz_prop_assoc_rerr","dz_prop_assoc_rerr",100, 0.,1.);
  dx_pull_prop_assoc = tfs->make<TH1F>("dx_pull_prop_assoc","dx_pull_prop_assoc",100,-10,10);
  dy_pull_prop_assoc = tfs->make<TH1F>("dy_pull_prop_assoc","dy_pull_prop_assoc",100,-10,10);
  dz_pull_prop_assoc = tfs->make<TH1F>("dz_pull_prop_assoc","dz_pull_prop_assoc",100,-10,10);
  //
  duxr_assoc          = tfs->make<TH1F>("duxr_assoc         ","duxr_assoc         ",100,-0.25,0.25);
  duyr_assoc          = tfs->make<TH1F>("duyr_assoc         ","duyr_assoc         ",100,-0.25,0.25);
  duzr_assoc          = tfs->make<TH1F>("duzr_assoc         ","duzr_assoc         ",100,-0.25,0.25);
  duxr_prop_assoc     = tfs->make<TH1F>("duxr_prop_assoc    ","duxr_prop_assoc    ",100,-0.25,0.25);
  duyr_prop_assoc     = tfs->make<TH1F>("duyr_prop_assoc    ","duyr_prop_assoc    ",100,-0.25,0.25);
  duzr_prop_assoc     = tfs->make<TH1F>("duzr_prop_assoc    ","duzr_prop_assoc    ",100,-0.25,0.25);
  dux_prop_assoc_rerr = tfs->make<TH1F>("dux_prop_assoc_rerr","dux_prop_assoc_rerr",100, 0.,1.);
  duy_prop_assoc_rerr = tfs->make<TH1F>("duy_prop_assoc_rerr","duy_prop_assoc_rerr",100, 0.,1.);
  duz_prop_assoc_rerr = tfs->make<TH1F>("duz_prop_assoc_rerr","duz_prop_assoc_rerr",100, 0.,1.);
  dux_pull_prop_assoc = tfs->make<TH1F>("dux_pull_prop_assoc","dux_pull_prop_assoc",100,-10,10);
  duy_pull_prop_assoc = tfs->make<TH1F>("duy_pull_prop_assoc","duy_pull_prop_assoc",100,-10,10);
  duz_pull_prop_assoc = tfs->make<TH1F>("duz_pull_prop_assoc","duz_pull_prop_assoc",100,-10,10);
  //
  dp0_prop_assoc      = tfs->make<TH1F>("dp0_prop_assoc     ","dp0_prop_assoc     ",100,-5,5);
  dp1_prop_assoc      = tfs->make<TH1F>("dp1_prop_assoc     ","dp1_prop_assoc     ",100,-10,10);
  dp2_prop_assoc      = tfs->make<TH1F>("dp2_prop_assoc     ","dp2_prop_assoc     ",100,-0.25,0.25);
  dp3_prop_assoc      = tfs->make<TH1F>("dp3_prop_assoc     ","dp3_prop_assoc     ",100,-0.25,0.25);
  dp0_prop_assoc_rerr = tfs->make<TH1F>("dp0_prop_assoc_rerr","dp0_prop_assoc_rerr",100, 0.,1.);
  dp1_prop_assoc_rerr = tfs->make<TH1F>("dp1_prop_assoc_rerr","dp1_prop_assoc_rerr",100, 0.,1.);
  dp2_prop_assoc_rerr = tfs->make<TH1F>("dp2_prop_assoc_rerr","dp2_prop_assoc_rerr",100, 0.,1.);
  dp3_prop_assoc_rerr = tfs->make<TH1F>("dp3_prop_assoc_rerr","dp3_prop_assoc_rerr",100, 0.,1.);
  dp0_pull_prop_assoc = tfs->make<TH1F>("dp0_pull_prop_assoc","dp0_pull_prop_assoc",100,-10,10);
  dp1_pull_prop_assoc = tfs->make<TH1F>("dp1_pull_prop_assoc","dp1_pull_prop_assoc",100,-10,10);
  dp2_pull_prop_assoc = tfs->make<TH1F>("dp2_pull_prop_assoc","dp2_pull_prop_assoc",100,-10,10);
  dp3_pull_prop_assoc = tfs->make<TH1F>("dp3_pull_prop_assoc","dp3_pull_prop_assoc",100,-10,10);
  //
  x_vs_prop_pullX       = tfs->make<TH2F>("x_vs_prop_pullX      ", "x_vs_prop_pullX      ", 100, -10, 10, 100, 0, 1000);
  y_vs_prop_pullX       = tfs->make<TH2F>("y_vs_prop_pullX      ", "y_vs_prop_pullX      ", 100, -10, 10, 100, -500, 500);
  z_vs_prop_pullX       = tfs->make<TH2F>("z_vs_prop_pullX      ", "z_vs_prop_pullX      ", 100, -10, 10, 100, 0, 1000);
  dirx_vs_prop_pullX    = tfs->make<TH2F>("dirx_vs_prop_pullX   ", "dirx_vs_prop_pullX   ", 100, -10, 10, 10, -1, 1);
  diry_vs_prop_pullX    = tfs->make<TH2F>("diry_vs_prop_pullX   ", "diry_vs_prop_pullX   ", 100, -10, 10, 10, -1, 1);
  dirz_vs_prop_pullX    = tfs->make<TH2F>("dirz_vs_prop_pullX   ", "dirz_vs_prop_pullX   ", 100, -10, 10, 10, -1, 1);
  p_vs_prop_pullX       = tfs->make<TH2F>("p_vs_prop_pullX      ", "p_vs_prop_pullX      ", 100, -10, 10, 10, 0, 5);
  theta_vs_prop_pullX   = tfs->make<TH2F>("theta_vs_prop_pullX  ", "theta_vs_prop_pullX  ", 100, -10, 10, 35, 0, 3.5);
  phi_vs_prop_pullX     = tfs->make<TH2F>("phi_vs_prop_pullX    ", "phi_vs_prop_pullX    ", 100, -10, 10, 35, -3.5, 3.5);
  zenith_vs_prop_pullX  = tfs->make<TH2F>("zenith_vs_prop_pullX ", "zenith_vs_prop_pullX ", 100, -10, 10, 35, 0, 3.5);
  azimuth_vs_prop_pullX = tfs->make<TH2F>("azimuth_vs_prop_pullX", "azimuth_vs_prop_pullX", 100, -10, 10, 35, -3.5, 3.5);
  //
  hit_track_res  = tfs->make<TH1F>("hit_track_res ","hit_track_res ",100,-0.5,0.5);
  hit_track_pull = tfs->make<TH1F>("hit_track_pull","hit_track_pull",100,-5,5);
  hit_track_terr = tfs->make<TH1F>("hit_track_terr","hit_track_terr",200,0,0.5);
  hit_track_herr = tfs->make<TH1F>("hit_track_herr","hit_track_herr",200,0,0.5);
  hit_track_plane0_res  = tfs->make<TH1F>("hit_track_plane0_res ","hit_track_plane0_res ",100,-0.5,0.5);
  hit_track_plane0_pull = tfs->make<TH1F>("hit_track_plane0_pull","hit_track_plane0_pull",100,-5,5);
  hit_track_plane0_terr = tfs->make<TH1F>("hit_track_plane0_terr","hit_track_plane0_terr",200,0,0.5);
  hit_track_plane0_herr = tfs->make<TH1F>("hit_track_plane0_herr","hit_track_plane0_herr",200,0,0.5);
  hit_track_plane1_res  = tfs->make<TH1F>("hit_track_plane1_res ","hit_track_plane1_res ",100,-0.5,0.5);
  hit_track_plane1_pull = tfs->make<TH1F>("hit_track_plane1_pull","hit_track_plane1_pull",100,-5,5);
  hit_track_plane1_terr = tfs->make<TH1F>("hit_track_plane1_terr","hit_track_plane1_terr",200,0,0.5);
  hit_track_plane1_herr = tfs->make<TH1F>("hit_track_plane1_herr","hit_track_plane1_herr",200,0,0.5);
  hit_track_plane2_res  = tfs->make<TH1F>("hit_track_plane2_res ","hit_track_plane2_res ",100,-0.5,0.5);
  hit_track_plane2_pull = tfs->make<TH1F>("hit_track_plane2_pull","hit_track_plane2_pull",100,-5,5);
  hit_track_plane2_terr = tfs->make<TH1F>("hit_track_plane2_terr","hit_track_plane2_terr",200,0,0.5);
  hit_track_plane2_herr = tfs->make<TH1F>("hit_track_plane2_herr","hit_track_plane2_herr",200,0,0.5);
  hit_track_plane0_res_vs_wire = tfs->make<TH2F>("hit_track_plane0_res_vs_wire","hit_track_plane0_res_vs_wire",2500,0,2500,100,-0.5,0.5);
  hit_track_plane1_res_vs_wire = tfs->make<TH2F>("hit_track_plane1_res_vs_wire","hit_track_plane1_res_vs_wire",2500,0,2500,100,-0.5,0.5);
  hit_track_plane2_res_vs_wire = tfs->make<TH2F>("hit_track_plane2_res_vs_wire","hit_track_plane2_res_vs_wire",3500,0,3500,100,-0.5,0.5);
  hit_track_plane0_res_vs_wire_test = tfs->make<TH2F>("hit_track_plane0_res_vs_wire_test","hit_track_plane0_res_vs_wire_test",2500,0,2500,100,-0.5,0.5);
  hit_track_plane1_res_vs_wire_test = tfs->make<TH2F>("hit_track_plane1_res_vs_wire_test","hit_track_plane1_res_vs_wire_test",2500,0,2500,100,-0.5,0.5);
  hit_track_plane2_res_vs_wire_test = tfs->make<TH2F>("hit_track_plane2_res_vs_wire_test","hit_track_plane2_res_vs_wire_test",3500,0,3500,100,-0.5,0.5);
  //
}

TrackValidation::TrackValidation(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  inputTracksLabel(p.get<std::string>("inputTracksLabel")),
  minHits(p.get<unsigned int>("minHits"))
{}

void TrackValidation::analyze(art::Event const & e)
{

  using namespace std;
  using namespace trkf;
  using namespace recob::tracking;

  const sim::MCTrack* mctk = 0;
  const recob::Track* retk = 0;

  detinfo::DetectorProperties const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
  detinfo::DetectorClocks const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  trkf::PropXYZPlane prop_(10., true);

  if (e.isRealData()==0) {
    art::InputTag SimTrackInputTag("mcreco");
    art::ValidHandle<std::vector<sim::MCTrack> > simTracks = e.getValidHandle<std::vector<sim::MCTrack> >(SimTrackInputTag);
    for (unsigned int iMC = 0; iMC < simTracks->size(); ++iMC) {
      const sim::MCTrack& mctrack = simTracks->at(iMC);
      if (mctrack.PdgCode()!=13)        continue;
      if (mctrack.Process()!="primary") continue;
      mctk = &mctrack;
      break;
    }
  }
  
  art::InputTag TrackInputTag(inputTracksLabel);
  art::ValidHandle<std::vector<recob::Track> > Tracks = e.getValidHandle<std::vector<recob::Track> >(TrackInputTag);
  // std::unique_ptr<art::FindManyP<recob::Hit, recob::TrackHitMeta> > pHits(new art::FindManyP<recob::Hit, recob::TrackHitMeta>(Tracks, e, TrackInputTag));
  // auto const& tkHitsAssn = *e.getValidHandle<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >(TrackInputTag);
  std::unique_ptr<art::FindManyP<recob::Hit> > pHits(new art::FindManyP<recob::Hit>(Tracks, e, TrackInputTag));
  auto const& tkHitsAssn = *e.getValidHandle<art::Assns<recob::Track, recob::Hit> >(TrackInputTag);

  art::ValidHandle<std::vector<std::vector<recob::TrackFitHitInfo> > > AllResiduals =
    e.getValidHandle<std::vector<std::vector<recob::TrackFitHitInfo> > >(TrackInputTag);
  
  //cout << "Tracks->size()=" << Tracks->size() << endl;
  all_NTracks->Fill(Tracks->size());

  unsigned int maxValidPoints = 0;
  for (unsigned int iTrack = 0; iTrack < Tracks->size(); ++iTrack) {
    const recob::Track& track = Tracks->at(iTrack);
    if (track.NumberTrajectoryPoints()<minHits) continue;
    all_NumberTrajectoryPoints->Fill(track.NumberTrajectoryPoints());
    all_CountValidPoints      ->Fill(track.CountValidPoints      ());
    all_HasMomentum           ->Fill(track.HasMomentum           ());
    all_Length                ->Fill(track.Length                ());
    all_Chi2                  ->Fill(track.Chi2                  ());
    all_Chi2PerNdof           ->Fill(track.Chi2PerNdof           ());
    all_Ndof                  ->Fill(track.Ndof                  ());
    all_ParticleId            ->Fill(track.ParticleId            ());
    all_Theta                 ->Fill(track.Theta                 ());
    all_Phi                   ->Fill(track.Phi                   ());
    all_ZenithAngle           ->Fill(track.ZenithAngle           ());
    all_AzimuthAngle          ->Fill(track.AzimuthAngle          ());
    all_NumberCovariance      ->Fill(track.NumberCovariance      ());
    //
    const auto start = track.Trajectory().Start();
    const auto end   = track.Trajectory().End();
    all_XYStart->Fill(start.X(),start.Y());
    all_XYEnd  ->Fill(end.X()  ,end.Y()  );
    all_ZRStart->Fill(start.Z(),std::sqrt(start.Perp2()));
    all_ZREnd  ->Fill(end.Z()  ,std::sqrt(end.Perp2()));
    //
    auto trackhits = pHits->at(iTrack);
    std::vector<recob::TrackFitHitInfo> residuals = AllResiduals->at(iTrack);
    assert(trackhits.size()==residuals.size());
    float nhits = trackhits.size();
    float nhits0 = 0.;
    float nhits1 = 0.;
    float nhits2 = 0.;
    // trick to preserve the order
    art::Ptr<recob::Track> ptrack(Tracks, iTrack);
    std::vector<art::Ptr<recob::Hit> > inHits;
    for (auto it = tkHitsAssn.begin(); it!=tkHitsAssn.end(); ++it) {
      if (it->first == ptrack) inHits.push_back(it->second);
      else if (inHits.size()>0) break;
    }
    //
    for (int ihit=0;ihit<nhits;ihit++) {
      auto hit = inHits[ihit];
      auto res = residuals[ihit];
      if (track.FlagsAtPoint(ihit).isExcludedFromFit()) continue;
      //std::cout << hit->WireID().Plane << ", " << hit->WireID().Wire << " --- " << res.WireId().Plane << ", " << res.WireId().Wire << std::endl;
      assert(hit->WireID().Plane==res.WireId().Plane);
      assert(hit->WireID().Wire==res.WireId().Wire);
      if (hit->WireID().Plane==0) {
	nhits0++;
	hit_track_plane0_res->Fill(res.trackStatePar()(0)-res.hitMeas());
	hit_track_plane0_pull->Fill( (res.trackStatePar()(0)-res.hitMeas())/sqrt(res.trackStateCov()(0,0)+res.hitMeasErr2()) );
	hit_track_plane0_res_vs_wire->Fill(res.WireId().Wire, res.trackStatePar()(0)-res.hitMeas());
	hit_track_plane0_terr->Fill( sqrt(res.trackStateCov()(0,0)) );
	hit_track_plane0_herr->Fill( sqrt(res.hitMeasErr2()) );
      }
      if (hit->WireID().Plane==1) {
	nhits1++;
	hit_track_plane1_res->Fill(res.trackStatePar()(0)-res.hitMeas());
	hit_track_plane1_pull->Fill( (res.trackStatePar()(0)-res.hitMeas())/sqrt(res.trackStateCov()(0,0)+res.hitMeasErr2()) );
	hit_track_plane1_res_vs_wire->Fill(res.WireId().Wire, res.trackStatePar()(0)-res.hitMeas());
	hit_track_plane1_terr->Fill( sqrt(res.trackStateCov()(0,0)) );
	hit_track_plane1_herr->Fill( sqrt(res.hitMeasErr2()) );
      }
      if (hit->WireID().Plane==2) {
	nhits2++;
	hit_track_plane2_res->Fill(res.trackStatePar()(0)-res.hitMeas());
	hit_track_plane2_pull->Fill( (res.trackStatePar()(0)-res.hitMeas())/sqrt(res.trackStateCov()(0,0)+res.hitMeasErr2()) );
	hit_track_plane2_res_vs_wire->Fill(res.WireId().Wire, res.trackStatePar()(0)-res.hitMeas());
	hit_track_plane2_terr->Fill( sqrt(res.trackStateCov()(0,0)) );
	hit_track_plane2_herr->Fill( sqrt(res.hitMeasErr2()) );
      }
      art::ServiceHandle<geo::Geometry> geom;
      recob::tracking::Plane plane = makePlane(geom->WireIDToWireGeo(res.WireId()));
      trkf::TrackState tkstate(res.trackStatePar(),res.trackStateCov(),plane,track.Trajectory().DirectionAtPoint(ihit).Dot(plane.direction())>0,13);
      // if (start.X()<50 && end.X()>200. && start.Z()<200 && end.Z()>800. ) {
      // if (tkstate.position().Z()>200. && tkstate.position().Z()<800. && fabs(tkstate.position().Y())<50. && tkstate.position().X()>230. && tkstate.momentum().Unit().X()<-0.5 && track.Length()>40.) {
      // if (tkstate.position().X()>100. && tkstate.position().X()<200. && fabs(tkstate.position().Y())<50. && track.Length()>40. && tkstate.momentum().Unit().X()>0. &&  tkstate.momentum().Unit().Z()>0.) {
      // if (start.X()>10 && start.X()<50 && end.X()>250.) {
      if (ihit==0) {
	if (hit->WireID().Plane==0) hit_track_plane0_res_vs_wire_test->Fill(res.WireId().Wire, res.trackStatePar()(0)-res.hitMeas());
	if (hit->WireID().Plane==1) hit_track_plane1_res_vs_wire_test->Fill(res.WireId().Wire, res.trackStatePar()(0)-res.hitMeas());
	if (hit->WireID().Plane==2) hit_track_plane2_res_vs_wire_test->Fill(res.WireId().Wire, res.trackStatePar()(0)-res.hitMeas());
      }
      // //std::cout << (res.trackStatePar()(0)-res.hitMeas())/(0.5*res.trackStatePar()(0)+0.5*res.hitMeas()) << std::endl;
      hit_track_res->Fill(res.trackStatePar()(0)-res.hitMeas());
      hit_track_pull->Fill( (res.trackStatePar()(0)-res.hitMeas())/sqrt(res.trackStateCov()(0,0)+res.hitMeasErr2()) );
      hit_track_terr->Fill( sqrt(res.trackStateCov()(0,0)) );
      hit_track_herr->Fill( sqrt(res.hitMeasErr2()) );
      //
      //test jacobians
      // auto par5 = res.trackStatePar();
      // auto cov5 = res.trackStateCov();
      // // for (int i=0;i<5;++i) for (int j=0;j<5;++j) cov5(i,j) = 1+i/3.+j/2.;
      // art::ServiceHandle<geo::Geometry> geom;
      // recob::tracking::Plane plane = makePlane(geom->WireIDToWireGeo(res.WireId()));
      // auto cov6newT = plane.Local5DToGlobal6DCovariance(cov5,true ,track.Trajectory().MomentumVectorAtPoint(ihit));
      // // auto cov6newF = plane.Local5DToGlobal6DCovariance(cov5,false,track.Trajectory().MomentumVectorAtPoint(ihit));
      // auto cov5newT = plane.Global6DToLocal5DCovariance(cov6newT,true ,track.Trajectory().MomentumVectorAtPoint(ihit));
      // // auto cov5newF = plane.Global6DToLocal5DCovariance(cov6newF,false,track.Trajectory().MomentumVectorAtPoint(ihit));
      // std::cout << "local cov=\n" << cov5 << std::endl;
      // std::cout << "global covnewT=\n" << cov6newT << std::endl;
      // std::cout << "local covnewT=\n" << cov5newT << std::endl;
      // // std::cout << "local covnewF =\n" << cov5newF  << std::endl;
      // bool alongdir = track.Trajectory().DirectionAtPoint(ihit).Dot(plane.direction())>0;
      // auto par6 = plane.Local5DToGlobal6DParameters(par5, alongdir);
      // auto par5new = plane.Global6DToLocal5DParameters(par6);
      // std::cout << "local par=" << par5 << std::endl;
      // std::cout << "global parnewT=" << par6 << std::endl;
      // std::cout << "local parnew=" << par5new << std::endl;
      // std::cout << "alongdir=" << alongdir << std::endl;
      //
    }
    all_hitPlane0Frac->Fill(nhits0/nhits);
    all_hitPlane1Frac->Fill(nhits1/nhits);
    all_hitPlane2Frac->Fill(nhits2/nhits);
    //
    if (track.CountValidPoints()>maxValidPoints) {
      maxValidPoints = track.CountValidPoints();
      retk = &track;
    }
  }

  if (mctk==0) cout << "no sim::MCTrack found" << endl;
  if (retk==0) cout << "no recob::Track found" << endl;
  if (mctk==0 || retk==0) return;

  double mclen = 0.;
  for (unsigned int imc=0; imc<mctk->size(); ++imc) {
    if (imc>0) {
      mclen+=sqrt( ((*mctk)[imc].X()-(*mctk)[imc-1].X())*((*mctk)[imc].X()-(*mctk)[imc-1].X()) +
		   ((*mctk)[imc].Y()-(*mctk)[imc-1].Y())*((*mctk)[imc].Y()-(*mctk)[imc-1].Y()) +
		   ((*mctk)[imc].Z()-(*mctk)[imc-1].Z())*((*mctk)[imc].Z()-(*mctk)[imc-1].Z()) );
    }
  }
    
  float dotvtx = retk->VertexDirection().X()*mctk->Start().Momentum().X()+retk->VertexDirection().Y()*mctk->Start().Momentum().Y()+retk->VertexDirection().Z()*mctk->Start().Momentum().Z();
  bool isSwapEndVtx = dotvtx<0;
  
  double g4Ticks = detClocks->TPCG4Time2Tick(mctk->Start().T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
  double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0);
  Point_t  mcpos(mctk->Start().Position().X()+xOffset,mctk->Start().Position().Y(),mctk->Start().Position().Z());
  Vector_t mcmom(mctk->Start().Momentum().X()*0.001,mctk->Start().Momentum().Y()*0.001,mctk->Start().Momentum().Z()*0.001);
  Vector_t mcdir = mcmom.Unit();

  SVector6 mcpar6(mcpos.X(), mcpos.Y(), mcpos.Z(), mcdir.X(), mcdir.Y(), mcdir.Z());
  SVector5 mcpar = Plane::Global6DToLocal5DParameters(mcpar6, mcpos, mcdir);
  const std::shared_ptr<const SurfXYZPlane> mcvtxsurf(new SurfXYZPlane(mcpos.X(), mcpos.Y(), mcpos.Z(), mcdir.X(), mcdir.Y(), mcdir.Z()));

  Point_t      repos = (isSwapEndVtx ?  retk->Trajectory().End()  : retk->Trajectory().Vertex() );
  Vector_t     redir = (isSwapEndVtx ? -retk->Trajectory().EndDirection() : retk->Trajectory().VertexDirection());
  SVector5     repar = (isSwapEndVtx ? retk->EndParametersLocal5D() : retk->VertexParametersLocal5D() );
  SMatrixSym55 recov = (isSwapEndVtx ? retk->EndCovarianceLocal5D() : retk->VertexCovarianceLocal5D() );
  const std::shared_ptr<const SurfXYZPlane> vtxsurf(new SurfXYZPlane(repos.X(), repos.Y(), repos.Z(), redir.X(), redir.Y(), redir.Z()));
  TrackVector vtxvec(5);
  for (int i=0;i<5;++i) vtxvec(i) = repar(i);
  vtxvec(4)=1./mcmom.R();
  KTrack vtxTrack(vtxsurf,vtxvec,Surface::FORWARD, 13);
  TrackError vtxerr; vtxerr.resize(5, false);
  for (int i=0;i<5;++i) for (int j=0;j<5;++j) vtxerr(i,j) = recov(i,j);
  KETrack tre(vtxTrack, vtxerr);
  prop_.noise_prop(tre,mcvtxsurf,Propagator::UNKNOWN,true);
  double  prop_pos[3];
  tre.getPosition(prop_pos);
  double  prop_mom[3];
  tre.getMomentum(prop_mom);
  Point_t      prpos(prop_pos[0],prop_pos[1],prop_pos[2]);
  Vector_t     prmom(prop_mom[0],prop_mom[1],prop_mom[2]);
  Vector_t     prdir = prmom.Unit();
  SVector5     prpar;
  SMatrixSym55 prcov;
  for (int i=0;i<5;++i) prpar(i) = tre.getVector()(i);
  for (int i=0;i<5;++i) for (int j=0;j<5;++j) prcov(i,j) = tre.getError()(i,j);

  if ( (prpos-mcpos).R()>10 ) {
    cout << "no associated tracks" << endl;
    return;
  }

  NumberTrajectoryPoints->Fill(retk->NumberTrajectoryPoints());
  CountValidPoints      ->Fill(retk->CountValidPoints      ());
  HasMomentum           ->Fill(retk->HasMomentum           ());
  Length                ->Fill(retk->Length                ());
  Chi2                  ->Fill(retk->Chi2                  ());
  Chi2PerNdof           ->Fill(retk->Chi2PerNdof           ());
  Ndof                  ->Fill(retk->Ndof                  ());
  ParticleId            ->Fill(retk->ParticleId            ());
  Theta                 ->Fill(retk->Theta                 ());
  Phi                   ->Fill(retk->Phi                   ());
  ZenithAngle           ->Fill(retk->ZenithAngle           ());
  AzimuthAngle          ->Fill(retk->AzimuthAngle          ());
  NumberCovariance      ->Fill(retk->NumberCovariance      ());
  //
  dp0_prop_assoc->Fill( prpar(0)-mcpar(0) );
  dp1_prop_assoc->Fill( prpar(1)-mcpar(1) );
  dp2_prop_assoc->Fill( prpar(2)-mcpar(2) );
  dp3_prop_assoc->Fill( prpar(3)-mcpar(3) );
  //
  dx_assoc->Fill(repos.X()-mcpos.X());
  dy_assoc->Fill(repos.Y()-mcpos.Y());
  dz_assoc->Fill(repos.Z()-mcpos.Z());
  duxr_assoc->Fill( (redir.X()-mcdir.X())/mcdir.X() );
  duyr_assoc->Fill( (redir.Y()-mcdir.Y())/mcdir.Y() );
  duzr_assoc->Fill( (redir.Z()-mcdir.Z())/mcdir.Z() );
  //
  dx_prop_assoc->Fill(prpos.X()-mcpos.X());
  dy_prop_assoc->Fill(prpos.Y()-mcpos.Y());
  dz_prop_assoc->Fill(prpos.Z()-mcpos.Z());
  duxr_prop_assoc->Fill( (prdir.X()-mcdir.X())/mcdir.X() );
  duyr_prop_assoc->Fill( (prdir.Y()-mcdir.Y())/mcdir.Y() );
  duzr_prop_assoc->Fill( (prdir.Z()-mcdir.Z())/mcdir.Z() );
  //
  if (retk->NumberCovariance()>1) {
    dp0_pull_prop_assoc->Fill( ( prpar(0)-mcpar(0) )/sqrt(prcov(0,0)) );
    dp1_pull_prop_assoc->Fill( ( prpar(1)-mcpar(1) )/sqrt(prcov(1,1)) );
    dp2_pull_prop_assoc->Fill( ( prpar(2)-mcpar(2) )/sqrt(prcov(2,2)) );
    dp3_pull_prop_assoc->Fill( ( prpar(3)-mcpar(3) )/sqrt(prcov(3,3)) );
    //
    dp0_prop_assoc_rerr->Fill( sqrt(prcov(0,0))/fabs( prpar(0) ) );
    dp1_prop_assoc_rerr->Fill( sqrt(prcov(1,1))/fabs( prpar(1) ) );
    dp2_prop_assoc_rerr->Fill( sqrt(prcov(2,2))/fabs( prpar(2) ) );
    dp3_prop_assoc_rerr->Fill( sqrt(prcov(3,3))/fabs( prpar(3) ) );
    //
    SMatrix65     j = Plane::Local5DToGlobal6DJacobian(false,prdir,mcdir);
    SMatrixSym66 c6 = ROOT::Math::Similarity(j,prcov);
    dx_pull_prop_assoc->Fill( (prpos.X()-mcpos.X())/sqrt(c6(0,0)) );
    dy_pull_prop_assoc->Fill( (prpos.Y()-mcpos.Y())/sqrt(c6(1,1)) );
    dz_pull_prop_assoc->Fill( (prpos.Z()-mcpos.Z())/sqrt(c6(2,2)) );
    dx_prop_assoc_rerr->Fill( sqrt(c6(0,0))/fabs(prpos.X()) );
    dy_prop_assoc_rerr->Fill( sqrt(c6(1,1))/fabs(prpos.Y()) );
    dz_prop_assoc_rerr->Fill( sqrt(c6(2,2))/fabs(prpos.Z()) );
    dux_pull_prop_assoc->Fill( (prdir.X()-mcdir.X())/sqrt(c6(3,3)) );
    duy_pull_prop_assoc->Fill( (prdir.Y()-mcdir.Y())/sqrt(c6(4,4)) );
    duz_pull_prop_assoc->Fill( (prdir.Z()-mcdir.Z())/sqrt(c6(5,5)) );
    dux_prop_assoc_rerr->Fill( sqrt(c6(3,3))/fabs(prdir.X()) );
    duy_prop_assoc_rerr->Fill( sqrt(c6(4,4))/fabs(prdir.Y()) );
    duz_prop_assoc_rerr->Fill( sqrt(c6(5,5))/fabs(prdir.Z()) );
    //
    double xpull = (prpos.X()-mcpos.X())/sqrt(c6(0,0));
    x_vs_prop_pullX->Fill( xpull, mcpos.X() );
    y_vs_prop_pullX->Fill( xpull, mcpos.Y() );
    z_vs_prop_pullX->Fill( xpull, mcpos.Z() );
    dirx_vs_prop_pullX->Fill( xpull, mcdir.X() );
    diry_vs_prop_pullX->Fill( xpull, mcdir.Y() );
    dirz_vs_prop_pullX->Fill( xpull, mcdir.Z() );
    p_vs_prop_pullX->Fill( xpull, mcmom.R() );
    theta_vs_prop_pullX->Fill( xpull, retk->Theta() );
    phi_vs_prop_pullX->Fill( xpull, retk->Phi() );
    zenith_vs_prop_pullX->Fill( xpull, retk->ZenithAngle() );
    azimuth_vs_prop_pullX->Fill( xpull, retk->AzimuthAngle() );
    //
  }
  //
  dLength->Fill(retk->Length()-mclen);
  dLengthRel->Fill( (retk->Length()-mclen)/mclen );
}

DEFINE_ART_MODULE(TrackValidation)
