////////////////////////////////////////////////////////////////////////
// Class:       TrackingPerformance
// Plugin Type: analyzer (art v2_08_03)
// File:        TrackingPerformance_module.cc
//
// Generated at Thu Sep 28 21:54:55 2017 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_01_01.
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
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardata/RecoObjects/TrackStatePropagator.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"

#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "TH1.h"
#include "TH2.h"

class TrackingPerformance : public art::EDAnalyzer {
public:
  explicit TrackingPerformance(fhicl::ParameterSet const & p);

  TrackingPerformance(TrackingPerformance const &) = delete;
  TrackingPerformance(TrackingPerformance &&) = delete;
  TrackingPerformance & operator = (TrackingPerformance const &) = delete;
  TrackingPerformance & operator = (TrackingPerformance &&) = delete;

  void analyze(art::Event const & e) override;

  void beginJob() override;

  const simb::MCParticle* getAssocMCParticle(const std::vector<art::Ptr<recob::Hit> >&,
                                             const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& hittruth);

private:
  //
  std::string inputTracksLabel;
  unsigned int minHits;
  int selectPdgCode;
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
  TH1F* dux_assoc;
  TH1F* duy_assoc;
  TH1F* duz_assoc;
  TH1F* dux_prop_assoc;
  TH1F* duy_prop_assoc;
  TH1F* duz_prop_assoc;
  TH1F* dux_prop_assoc_rerr;
  TH1F* duy_prop_assoc_rerr;
  TH1F* duz_prop_assoc_rerr;
  TH1F* dux_pull_prop_assoc;
  TH1F* duy_pull_prop_assoc;
  TH1F* duz_pull_prop_assoc;
  //
  // TH1F* dp0_prop_assoc;
  // TH1F* dp1_prop_assoc;
  // TH1F* dp2_prop_assoc;
  // TH1F* dp3_prop_assoc;
  // TH1F* dp0_prop_assoc_rerr;
  // TH1F* dp1_prop_assoc_rerr;
  // TH1F* dp2_prop_assoc_rerr;
  // TH1F* dp3_prop_assoc_rerr;
  // TH1F* dp0_pull_prop_assoc;
  // TH1F* dp1_pull_prop_assoc;
  // TH1F* dp2_pull_prop_assoc;
  // TH1F* dp3_pull_prop_assoc;
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
  TH2F* mu_mcsmom_vs_truemom;
  TH2F* mu_mcsmom_vs_truemom_contained;
  //
  TH1F* ptypemc;
  TH1F* deltaP;
  TH1F* deltaPrel;
  TH1F* dRVtxMC;
  TH1F* nHits;
  TH1F* dotpdir;
  //
  TH2F* mom_vs_truemom;
  TH2F* mom_vs_truemom_okid;
  //
  TH2F* ptype_rec_vs_mc;
  TH1F* deltaP_okid;
  TH1F* deltaPrel_okid;
  //
};


TrackingPerformance::TrackingPerformance(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
  , inputTracksLabel(p.get<std::string>("inputTracksLabel"))
  , minHits(p.get<unsigned int>("minHits"))
  , selectPdgCode(p.get<int>("selectPdgCode"))
{}

void TrackingPerformance::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
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
  ParticleId             = tfs->make<TH1F>("ParticleId            ","ParticleId            ",   6, 0,    6);
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
  dx_pull_prop_assoc = tfs->make<TH1F>("dx_pull_prop_assoc","dx_pull_prop_assoc",100,-5,5);
  dy_pull_prop_assoc = tfs->make<TH1F>("dy_pull_prop_assoc","dy_pull_prop_assoc",100,-5,5);
  dz_pull_prop_assoc = tfs->make<TH1F>("dz_pull_prop_assoc","dz_pull_prop_assoc",100,-5,5);
  //
  dux_assoc          = tfs->make<TH1F>("dux_assoc         ","dux_assoc         ",100,-0.25,0.25);
  duy_assoc          = tfs->make<TH1F>("duy_assoc         ","duy_assoc         ",100,-0.25,0.25);
  duz_assoc          = tfs->make<TH1F>("duz_assoc         ","duz_assoc         ",100,-0.25,0.25);
  dux_prop_assoc     = tfs->make<TH1F>("dux_prop_assoc    ","dux_prop_assoc    ",100,-0.25,0.25);
  duy_prop_assoc     = tfs->make<TH1F>("duy_prop_assoc    ","duy_prop_assoc    ",100,-0.25,0.25);
  duz_prop_assoc     = tfs->make<TH1F>("duz_prop_assoc    ","duz_prop_assoc    ",100,-0.25,0.25);
  dux_prop_assoc_rerr = tfs->make<TH1F>("dux_prop_assoc_rerr","dux_prop_assoc_rerr",100, 0.,1.);
  duy_prop_assoc_rerr = tfs->make<TH1F>("duy_prop_assoc_rerr","duy_prop_assoc_rerr",100, 0.,1.);
  duz_prop_assoc_rerr = tfs->make<TH1F>("duz_prop_assoc_rerr","duz_prop_assoc_rerr",100, 0.,1.);
  dux_pull_prop_assoc = tfs->make<TH1F>("dux_pull_prop_assoc","dux_pull_prop_assoc",100,-5,5);
  duy_pull_prop_assoc = tfs->make<TH1F>("duy_pull_prop_assoc","duy_pull_prop_assoc",100,-5,5);
  duz_pull_prop_assoc = tfs->make<TH1F>("duz_pull_prop_assoc","duz_pull_prop_assoc",100,-5,5);
  //
  // dp0_prop_assoc      = tfs->make<TH1F>("dp0_prop_assoc     ","dp0_prop_assoc     ",100,-5,5);
  // dp1_prop_assoc      = tfs->make<TH1F>("dp1_prop_assoc     ","dp1_prop_assoc     ",100,-10,10);
  // dp2_prop_assoc      = tfs->make<TH1F>("dp2_prop_assoc     ","dp2_prop_assoc     ",100,-0.25,0.25);
  // dp3_prop_assoc      = tfs->make<TH1F>("dp3_prop_assoc     ","dp3_prop_assoc     ",100,-0.25,0.25);
  // dp0_prop_assoc_rerr = tfs->make<TH1F>("dp0_prop_assoc_rerr","dp0_prop_assoc_rerr",100, 0.,1.);
  // dp1_prop_assoc_rerr = tfs->make<TH1F>("dp1_prop_assoc_rerr","dp1_prop_assoc_rerr",100, 0.,1.);
  // dp2_prop_assoc_rerr = tfs->make<TH1F>("dp2_prop_assoc_rerr","dp2_prop_assoc_rerr",100, 0.,1.);
  // dp3_prop_assoc_rerr = tfs->make<TH1F>("dp3_prop_assoc_rerr","dp3_prop_assoc_rerr",100, 0.,1.);
  // dp0_pull_prop_assoc = tfs->make<TH1F>("dp0_pull_prop_assoc","dp0_pull_prop_assoc",100,-10,10);
  // dp1_pull_prop_assoc = tfs->make<TH1F>("dp1_pull_prop_assoc","dp1_pull_prop_assoc",100,-10,10);
  // dp2_pull_prop_assoc = tfs->make<TH1F>("dp2_pull_prop_assoc","dp2_pull_prop_assoc",100,-10,10);
  // dp3_pull_prop_assoc = tfs->make<TH1F>("dp3_pull_prop_assoc","dp3_pull_prop_assoc",100,-10,10);
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
  mu_mcsmom_vs_truemom = tfs->make<TH2F>("mu_mcsmom_vs_truemom", "mu_mcsmom_vs_truemom", 50, 0., 5., 50, 0., 5.);
  mu_mcsmom_vs_truemom_contained = tfs->make<TH2F>("mu_mcsmom_vs_truemom_contained", "mu_mcsmom_vs_truemom_contained", 50, 0., 5., 50, 0., 5.);
  //
  ptypemc = tfs->make<TH1F>("ptypemc","ptypemc",6,0,6);
  deltaP = tfs->make<TH1F>("deltaP","deltaP",20,-1,1);
  deltaPrel = tfs->make<TH1F>("deltaPrel","deltaPrel",20,-1,1);
  dRVtxMC = tfs->make<TH1F>("dRVtxMC","dRVtxMC",20,0,10);
  nHits = tfs->make<TH1F>("nHits","nHits",20,0,100);
  dotpdir = tfs->make<TH1F>("dotpdir","dotpdir",20,-1,1);
  //
  mom_vs_truemom = tfs->make<TH2F>("mom_vs_truemom", "mom_vs_truemom", 50, 0., 5., 50, 0., 5.);
  mom_vs_truemom_okid = tfs->make<TH2F>("mom_vs_truemom_okid", "mom_vs_truemom_okid", 50, 0., 5., 50, 0., 5.);
  //
  ptype_rec_vs_mc = tfs->make<TH2F>("ptype_rec_vs_mc","ptype_rec_vs_mc",6,0,6,6,0,6);
  deltaP_okid  = tfs->make<TH1F>("deltaP_okid","deltaP_okid",50,-1,1);
  deltaPrel_okid  = tfs->make<TH1F>("deltaPrel_okid","deltaPrel_okid",50,-1,1);
}

void TrackingPerformance::analyze(art::Event const & e)
{
  //
  using namespace std;
  using namespace recob;
  using namespace recob::tracking;
  //
  detinfo::DetectorProperties const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
  detinfo::DetectorClocks const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  //
  const auto& inputPFParticle = e.getValidHandle<vector<PFParticle> >("pandora");
  auto assocTracks = unique_ptr<art::FindManyP<Track> >(new art::FindManyP<Track>(inputPFParticle, e, inputTracksLabel));
  const auto& inputTracks = e.getValidHandle<vector<Track> >(inputTracksLabel);
  auto assocHits = unique_ptr<art::FindManyP<Hit> >(new art::FindManyP<Hit>(inputTracks, e, inputTracksLabel));
  const auto& inputHits = e.getValidHandle<vector<Hit> >("gaushit");
  const auto& hittruth = unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>(inputHits,e,"gaushitTruthMatch"));
  //
  auto const& mctruth = *e.getValidHandle<std::vector<simb::MCTruth> >("generator");
  //
  const auto& mcsmom = *e.getValidHandle<vector<MCSFitResult> >("pandoraMCSMu");
  const std::vector<anab::CosmicTag>* cont = e.getValidHandle<std::vector<anab::CosmicTag> >("pandoraContTag").product();

  Point_t nuvtx(mctruth[0].GetNeutrino().Nu().Position().X(),mctruth[0].GetNeutrino().Nu().Position().Y(),mctruth[0].GetNeutrino().Nu().Position().Z());
  if (0) std::cout << "nu vtx=" << nuvtx << " with daughters=" << mctruth[0].GetNeutrino().Nu().NumberDaughters() << std::endl;
  for (int i=0; i<mctruth[0].NParticles(); ++i) {
    if (mctruth[0].GetParticle(i).StatusCode()!=1) continue;
    if (0) cout << "part pdgid=" << mctruth[0].GetParticle(i).PdgCode() << " pos=(" << mctruth[0].GetParticle(i).Vx() << "," << mctruth[0].GetParticle(i).Vy() << "," << mctruth[0].GetParticle(i).Vz() << ") dir=(" << mctruth[0].GetParticle(i).Px()/mctruth[0].GetParticle(i).P() << "," << mctruth[0].GetParticle(i).Py()/mctruth[0].GetParticle(i).P() << "," << mctruth[0].GetParticle(i).Pz()/mctruth[0].GetParticle(i).P() << ") p=" << mctruth[0].GetParticle(i).P() << " status=" << mctruth[0].GetParticle(i).StatusCode() << " process=" << mctruth[0].GetParticle(i).Process() << endl;
  }

  for (size_t iPF = 0; iPF < inputPFParticle->size(); ++iPF) {
    art::Ptr<PFParticle> pfp(inputPFParticle, iPF);
    if (pfp->IsPrimary()==false || pfp->NumDaughters()<2) continue;
    if (0) std::cout << "pfp#" << iPF << " PdgCode=" << pfp->PdgCode()
	      << " IsPrimary=" << pfp->IsPrimary()
	      << " NumDaughters=" << pfp->NumDaughters()
	      << std::endl;
    auto& pfd = pfp->Daughters();
    for (auto ipfd : pfd) {
      vector< art::Ptr<Track> > pftracks = assocTracks->at(ipfd);
      art::Ptr<PFParticle> pfpd(inputPFParticle, ipfd);
      if (0) cout << "pfp id=" << ipfd << " pdg=" << pfpd->PdgCode() << endl;
      for (auto t : pftracks) {
	//
	if (selectPdgCode!=-1 && std::abs(t->ParticleId())!=selectPdgCode) {
	  if (0) cout << "track has not the selected pdgCode but=" << t->ParticleId() << endl;
	  continue;
	}
	nHits->Fill(t->CountValidPoints());
	if (t->CountValidPoints()<minHits) {
	  if (0) cout << "track has only npoints=" << t->CountValidPoints() << endl;
	  continue;
	}
	//
	const simb::MCParticle* mcp = getAssocMCParticle(assocHits->at(t.key()),hittruth);
	if (mcp==nullptr) continue;
	dRVtxMC->Fill((nuvtx-Point_t(mcp->Vx(),mcp->Vy(),mcp->Vz())).R());
	if (0) cout << "MCParticle status=" << mcp->StatusCode() << " Mother=" << mcp->Mother() << " Process=" << mcp->Process() /*<< " MotherPdgCode=" << (mcp->Mother()>=0 ? mctruth[0].GetParticle(mcp->Mother()).PdgCode() : 0)*/ << endl;
	if ((nuvtx-Point_t(mcp->Vx(),mcp->Vy(),mcp->Vz())).R()>0.5) {
	  if (0) cout << "track has matched particle not from neutrino" << endl;
	  continue;
	}
	//
	auto repos = t->Start();
	auto redir = t->StartDirection();
	//
	auto scecorr = SCE->GetPosOffsets(geo::Point_t(mcp->Vx(),mcp->Vy(),mcp->Vz()));
	double g4Ticks = detClocks->TPCG4Time2Tick(mcp->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
	double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
	double yOffset = scecorr.Y();
	double zOffset = scecorr.Z();
	Point_t mcpos(mcp->Vx()+xOffset,mcp->Vy()+yOffset,mcp->Vz()+zOffset);
	Vector_t mcdir(mcp->Px()/mcp->P(),mcp->Py()/mcp->P(),mcp->Pz()/mcp->P());
	//
	dotpdir->Fill(redir.Dot(mcdir));
	if (redir.Dot(mcdir)<0.5) {
	  if (0) cout << "track not in the right direction..." << endl;
	  continue;//fixme, if too frequent impose direction exiting from vertex
	}
	//
	int pidmc = 0;
	if (std::abs(mcp->PdgCode())==13)        pidmc = 1;
	else if (std::abs(mcp->PdgCode())==211)  pidmc = 2;
	else if (std::abs(mcp->PdgCode())==321)  pidmc = 3;
	else if (std::abs(mcp->PdgCode())==2212) pidmc = 4;
	else if (std::abs(mcp->PdgCode())==11)   pidmc = 5;
	ptypemc->Fill(pidmc);
	int pidr = 0;
	if (std::abs(t->ParticleId())==13)        pidr = 1;
	else if (std::abs(t->ParticleId())==211)  pidr = 2;
	else if (std::abs(t->ParticleId())==321)  pidr = 3;
	else if (std::abs(t->ParticleId())==2212) pidr = 4;
	else if (std::abs(t->ParticleId())==11)   pidr = 5;
	ptype_rec_vs_mc->Fill(pidmc,pidr);
	if (selectPdgCode!=-1 && std::abs(t->ParticleId())!=std::abs(mcp->PdgCode())) {
	  if (0) cout << "track has not the same pdgCode as MC:" << t->ParticleId() << " vs " << mcp->PdgCode() << endl;
	  continue;
	}
	//
	if (0) cout << "trk pos=" << t->Start() << " dir=" << t->StartDirection() << " p=" << t->VertexMomentum() << " pdg=" << t->ParticleId() << endl;
	if (0) cout << "mcp pos=(" << mcp->Vx() << "," << mcp->Vy() << "," << mcp->Vz() << ") dir=(" << mcp->Px()/mcp->P() << "," << mcp->Py()/mcp->P() << "," << mcp->Pz()/mcp->P() << ") p=" << mcp->P() << " pdg=" << mcp->PdgCode() << " mom=" << mcp->Mother() << " id=" << mcp->TrackId() << " proc=" << mcp->Process() << endl;
	//
	NumberTrajectoryPoints->Fill(t->NumberTrajectoryPoints());
	CountValidPoints->Fill(t->CountValidPoints());
	HasMomentum->Fill(t->HasMomentum());
	auto rl = t->Length();
	auto mcl = mcp->Trajectory().TotalLength();
	Length->Fill(rl);
	dLength->Fill(rl-mcl);
	dLengthRel->Fill( (rl-mcl)/mcl );
	Chi2->Fill(t->Chi2());
	Chi2PerNdof->Fill(t->Chi2PerNdof());
	Ndof->Fill(t->Ndof());
	if (std::abs(t->ParticleId())==13)        ParticleId->Fill(1);
	else if (std::abs(t->ParticleId())==211)  ParticleId->Fill(2);
	else if (std::abs(t->ParticleId())==321)  ParticleId->Fill(3);
	else if (std::abs(t->ParticleId())==2212) ParticleId->Fill(4);
	else if (std::abs(t->ParticleId())==11)   ParticleId->Fill(5);
	else ParticleId->Fill(0);
	Theta->Fill(t->Theta());
	Phi->Fill(t->Phi());
	ZenithAngle->Fill(t->ZenithAngle());
	AzimuthAngle->Fill(t->AzimuthAngle());
	NumberCovariance->Fill(t->NumberCovariance());
	deltaP->Fill( t->StartMomentum()-mcp->P() );
	deltaPrel->Fill( (t->StartMomentum()-mcp->P())/mcp->P() );
	mom_vs_truemom->Fill(mcp->P(),t->StartMomentum());
	if (std::abs(t->ParticleId())==std::abs(mcp->PdgCode())) {
	  deltaP_okid->Fill( t->StartMomentum()-mcp->P() );
	  deltaPrel_okid->Fill( (t->StartMomentum()-mcp->P())/mcp->P() );
	  mom_vs_truemom_okid->Fill(mcp->P(),t->StartMomentum());
	}
	//
	if (0) cout << "repos=" << repos << " endpos=" << t->Trajectory().End() << " mcpos=" << mcpos << " mcposorig=" << Point_t(mcp->Vx(),mcp->Vy(),mcp->Vz()) << endl;
	if (0) cout << "scecorr=" << scecorr << endl;
	if (0) cout << "redir=" << redir << " mcdir=" << mcdir << endl;
	if (0) cout << "dot=" << redir.Dot(mcdir) << endl;
	//
	dx_assoc->Fill(repos.X()-mcpos.X());
	dy_assoc->Fill(repos.Y()-mcpos.Y());
	dz_assoc->Fill(repos.Z()-mcpos.Z());
	dux_assoc->Fill(redir.X()-mcdir.X());
	duy_assoc->Fill(redir.Y()-mcdir.Y());
	duz_assoc->Fill(redir.Z()-mcdir.Z());
	//
	trkf::TrackStatePropagator prop(1.0, 0.1, 10, 10., 0.01, false);
	recob::tracking::Plane mcplane(mcpos,mcdir);
	recob::tracking::Plane tkplane(repos,redir);
	trkf::TrackState tkstart(t->VertexParametersLocal5D(), t->VertexCovarianceLocal5D(), tkplane, true, t->ParticleId());
	bool propok = true;
	trkf::TrackState tkatmc = prop.propagateToPlane(propok, tkstart, mcplane, true, true, trkf::TrackStatePropagator::UNKNOWN);
	//
	if (!propok) continue;
	//
	auto prpos = tkatmc.position();
	if (0) cout << "prpos=" << prpos << endl;
	auto prdir = tkatmc.momentum().Unit();
	auto c = tkatmc.plane().Local5DToGlobal6DCovariance(tkatmc.covariance(),false,tkatmc.momentum());//consider direction not momentum
	// auto cv = pfp.Track.VertexCovarianceLocal5D();
	//
	dx_prop_assoc->Fill(prpos.X()-mcpos.X());
	dy_prop_assoc->Fill(prpos.Y()-mcpos.Y());
	dz_prop_assoc->Fill(prpos.Z()-mcpos.Z());
	dux_prop_assoc->Fill(prdir.X()-mcdir.X());
	duy_prop_assoc->Fill(prdir.Y()-mcdir.Y());
	duz_prop_assoc->Fill(prdir.Z()-mcdir.Z());
	//
	dx_prop_assoc_rerr->Fill( sqrt(c(0,0))/fabs(prpos.X()) );
	dy_prop_assoc_rerr->Fill( sqrt(c(1,1))/fabs(prpos.Y()) );
	dz_prop_assoc_rerr->Fill( sqrt(c(2,2))/fabs(prpos.Z()) );
	dux_prop_assoc_rerr->Fill( sqrt(c(3,3))/fabs(prdir.X()) );
	duy_prop_assoc_rerr->Fill( sqrt(c(4,4))/fabs(prdir.Y()) );
	duz_prop_assoc_rerr->Fill( sqrt(c(5,5))/fabs(prdir.Z()) );
	//
	dx_pull_prop_assoc->Fill( (prpos.X()-mcpos.X())/sqrt(c(0,0)) );
	dy_pull_prop_assoc->Fill( (prpos.Y()-mcpos.Y())/sqrt(c(1,1)) );
	dz_pull_prop_assoc->Fill( (prpos.Z()-mcpos.Z())/sqrt(c(2,2)) );
	dux_pull_prop_assoc->Fill( (prdir.X()-mcdir.X())/sqrt(c(3,3)) );
	duy_pull_prop_assoc->Fill( (prdir.Y()-mcdir.Y())/sqrt(c(4,4)) );
	duz_pull_prop_assoc->Fill( (prdir.Z()-mcdir.Z())/sqrt(c(5,5)) );
	//
	auto pullX = (prpos.X()-mcpos.X())/sqrt(c(0,0));
	x_vs_prop_pullX->Fill(pullX,mcpos.X());
	y_vs_prop_pullX->Fill(pullX,mcpos.Y());
	z_vs_prop_pullX->Fill(pullX,mcpos.Z());
	dirx_vs_prop_pullX->Fill(pullX,mcdir.X());
	diry_vs_prop_pullX->Fill(pullX,mcdir.Y());
	dirz_vs_prop_pullX->Fill(pullX,mcdir.Z());
	p_vs_prop_pullX->Fill(pullX,mcp->P());
	theta_vs_prop_pullX->Fill(pullX,t->Theta());
	phi_vs_prop_pullX->Fill(pullX,t->Phi());
	zenith_vs_prop_pullX->Fill(pullX,t->ZenithAngle());
	azimuth_vs_prop_pullX->Fill(pullX,t->AzimuthAngle());
	//
        if (std::abs(mcp->PdgCode())==13) {
          //a muon, check mcs mom
          auto mcs = mcsmom[t.key()];
          if (0) std::cout << "muon with mc mom=" << mcp->P() << " mcs=" << mcs.bestMomentum() << std::endl;
          mu_mcsmom_vs_truemom->Fill(mcp->P(),mcs.bestMomentum());
	  //
	  bool isContained = cont->at(t->ID()).CosmicType()==anab::kNotTagged;
          if (isContained) mu_mcsmom_vs_truemom_contained->Fill(mcp->P(),mcs.bestMomentum());
        }
      }
    }
  }

}

const simb::MCParticle* TrackingPerformance::getAssocMCParticle(const std::vector<art::Ptr<recob::Hit> >& hits,
                                                                const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& hittruth) {
  //credit: Wes Ketchum
  std::unordered_map<int,double> trkide;
  double maxe=-1, tote=0;
  simb::MCParticle const* maxp_me = nullptr; //pointer for the particle match we will calculate
  for (auto h : hits) {
    std::vector<art::Ptr<simb::MCParticle> > particle_vec = hittruth->at(h.key());
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec = hittruth->data(h.key());;
    //loop over particles
    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
      trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
      tote += match_vec[i_p]->energy; //calculate total energy deposited
      if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
	maxe = trkide[ particle_vec[i_p]->TrackId() ];
	maxp_me = particle_vec[i_p].get();
      }
    }//end loop over particles per hit
  }
  return maxp_me;
}

DEFINE_ART_MODULE(TrackingPerformance)
