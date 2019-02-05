////////////////////////////////////////////////////////////////////////
// Class:       ShowerPerformance
// Plugin Type: analyzer (art v3_00_00)
// File:        ShowerPerformance_module.cc
//
// Generated at Wed Jan  9 10:02:44 2019 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_04_00.
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
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardata/RecoObjects/TrackStatePropagator.h"

#include "lardata/RecoBaseProxy/Track.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/GeometryUtilities.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/Utilities/ForEachAssociatedGroup.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

class ShowerPerformance;

class ShowerPerformance : public art::EDAnalyzer {
public:
  explicit ShowerPerformance(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShowerPerformance(ShowerPerformance const&) = delete;
  ShowerPerformance(ShowerPerformance&&) = delete;
  ShowerPerformance& operator=(ShowerPerformance const&) = delete;
  ShowerPerformance& operator=(ShowerPerformance&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginJob() override;
  void resetTree();

  art::Ptr<simb::MCParticle> getAssocMCParticle(const std::vector<art::Ptr<recob::Hit> >&,
						const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& hittruth) const;

  int nHitsFromMCParticle(size_t mcid,
			  const std::vector<art::Ptr<recob::Hit> >&,
			  const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& hittruth) const;

  template <typename T> std::vector<float> dEdx(const recob::Shower& shower, const std::vector<art::Ptr<recob::Cluster> >& clusters, const T& clusterHitsGroups, float dtrunk);

  template <typename T> std::vector<float> dEdx(const recob::Track& track, const T& hits, float range);

private:
  //
  std::string inputShowerLabel;
  std::string inputTrackLabel;

  TH1F* Length;
  //
  TH1F* dx_assoc;
  TH1F* dy_assoc;
  TH1F* dz_assoc;
  //
  TH1F* dux_assoc;
  TH1F* duy_assoc;
  TH1F* duz_assoc;
  //
  TH1F* dotpdir;
  //
  TH1F* purity;
  TH1F* completeness;
  //
  TH1F* dedxpl2rng4;

  TH1F* tk_Length;
  //
  TH1F* tk_dx_assoc;
  TH1F* tk_dy_assoc;
  TH1F* tk_dz_assoc;
  //
  TH1F* tk_dux_assoc;
  TH1F* tk_duy_assoc;
  TH1F* tk_duz_assoc;
  //
  TH1F* tk_dx_assoc_pull;
  TH1F* tk_dy_assoc_pull;
  TH1F* tk_dz_assoc_pull;
  //
  TH1F* tk_dux_assoc_pull;
  TH1F* tk_duy_assoc_pull;
  TH1F* tk_duz_assoc_pull;
  //
  TH1F* tk_dotpdir;
  //
  TH1F* tk_purity;
  TH1F* tk_completeness;
  //
  TH1F* tk_dedxpl2rng4;

  //
  TTree* tree;
  //
  int    run, subrun, eventid;
  //
  int mce_ispr;
  int mce_key;
  int mce_nh;
  int mce_isconv;
  int mce_pdg;
  int mce_mompdg;
  int mce_momispr;
  float mce_x, mce_y, mce_z, mce_ux, mce_uy, mce_uz;
  float mce_e/*, mce_dedx_r4p0, mce_dedx_r4p1, mce_dedx_r4p2, mce_dedx_r2p0, mce_dedx_r2p1, mce_dedx_r2p2*/;
  //
  float shw_len;
  int shw_nh;
  float shw_x, shw_y, shw_z, shw_ux, shw_uy, shw_uz;
  float shw_e_p0, shw_e_p1, shw_e_p2, shw_dedx_r4p0, shw_dedx_r4p1, shw_dedx_r4p2, shw_dedx_r2p0, shw_dedx_r2p1, shw_dedx_r2p2;
  float shw_dedx_defp0, shw_dedx_defp1, shw_dedx_defp2;
  float shw_pur, shw_cmpl, shw_mce_dotpdir;
  int shw_mce_key, shw_mce_ispr;
  int shw_mce_isconv;
  int shw_mce_pdg;
  int shw_mce_mompdg;
  int shw_mce_momispr;
  //
  float trk_len;
  int trk_nh, trk_nvh, trk_nvh_p0, trk_nvh_p1, trk_nvh_p2;
  int trk_nvh4, trk_nvh4_p0, trk_nvh4_p1, trk_nvh4_p2;
  float trk_x, trk_y, trk_z, trk_ux, trk_uy, trk_uz;
  float trk_ex, trk_ey, trk_ez, trk_eux, trk_euy, trk_euz;
  float trk_chi2, trk_ndof;
  float trk_dedx_r4p0, trk_dedx_r4p1, trk_dedx_r4p2, trk_dedx_r2p0, trk_dedx_r2p1, trk_dedx_r2p2;
  float trk_pur, trk_cmpl, trk_mce_dotpdir;
  int trk_mce_key, trk_mce_ispr;
  int trk_mce_isconv;
  int trk_mce_pdg;
  int trk_mce_mompdg;
  int trk_mce_momispr;
  //
};


ShowerPerformance::ShowerPerformance(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , inputShowerLabel(p.get<std::string>("inputShowerLabel"))
  , inputTrackLabel(p.get<std::string>("inputTrackLabel"))
  // More initializers here.
{}

void ShowerPerformance::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  //
  Length = tfs->make<TH1F>("Length","Length", 100,  0, 500);
  //
  dx_assoc = tfs->make<TH1F>("dx_assoc","dx_assoc",100,-2.5,2.5);
  dy_assoc = tfs->make<TH1F>("dy_assoc","dy_assoc",100,-2.5,2.5);
  dz_assoc = tfs->make<TH1F>("dz_assoc","dz_assoc",100,-2.5,2.5);
  //
  dux_assoc = tfs->make<TH1F>("dux_assoc","dux_assoc",100,-0.25,0.25);
  duy_assoc = tfs->make<TH1F>("duy_assoc","duy_assoc",100,-0.25,0.25);
  duz_assoc = tfs->make<TH1F>("duz_assoc","duz_assoc",100,-0.25,0.25);
  //
  dotpdir = tfs->make<TH1F>("dotpdir","dotpdir",20,-1,1);
  //
  purity = tfs->make<TH1F>("purity","purity",21,0,1.05);
  completeness = tfs->make<TH1F>("completeness","completeness",21,0,1.05);
  //
  dedxpl2rng4 = tfs->make<TH1F>("dedxpl2rng4","dedxpl2rng4",100,0,20);

  //
  tk_Length = tfs->make<TH1F>("tk_Length","tk_Length", 100,  0, 500);
  //
  tk_dx_assoc = tfs->make<TH1F>("tk_dx_assoc","tk_dx_assoc",100,-2.5,2.5);
  tk_dy_assoc = tfs->make<TH1F>("tk_dy_assoc","tk_dy_assoc",100,-2.5,2.5);
  tk_dz_assoc = tfs->make<TH1F>("tk_dz_assoc","tk_dz_assoc",100,-2.5,2.5);
  //
  tk_dux_assoc = tfs->make<TH1F>("tk_dux_assoc","tk_dux_assoc",100,-0.25,0.25);
  tk_duy_assoc = tfs->make<TH1F>("tk_duy_assoc","tk_duy_assoc",100,-0.25,0.25);
  tk_duz_assoc = tfs->make<TH1F>("tk_duz_assoc","tk_duz_assoc",100,-0.25,0.25);
  //
  tk_dx_assoc_pull = tfs->make<TH1F>("tk_dx_assoc_pull","tk_dx_assoc_pull",100,-5,5);
  tk_dy_assoc_pull = tfs->make<TH1F>("tk_dy_assoc_pull","tk_dy_assoc_pull",100,-5,5);
  tk_dz_assoc_pull = tfs->make<TH1F>("tk_dz_assoc_pull","tk_dz_assoc_pull",100,-5,5);
  //
  tk_dux_assoc_pull = tfs->make<TH1F>("tk_dux_assoc_pull","tk_dux_assoc_pull",100,-5,5);
  tk_duy_assoc_pull = tfs->make<TH1F>("tk_duy_assoc_pull","tk_duy_assoc_pull",100,-5,5);
  tk_duz_assoc_pull = tfs->make<TH1F>("tk_duz_assoc_pull","tk_duz_assoc_pull",100,-5,5);
  //
  tk_dotpdir = tfs->make<TH1F>("tk_dotpdir","tk_dotpdir",20,-1,1);
  //
  tk_purity = tfs->make<TH1F>("tk_purity","tk_purity",21,0,1.05);
  tk_completeness = tfs->make<TH1F>("tk_completeness","tk_completeness",21,0,1.05);
  //
  tk_dedxpl2rng4 = tfs->make<TH1F>("tk_dedxpl2rng4","tk_dedxpl2rng4",100,0,20);

  //
  tree = tfs->make<TTree>("tree", "tree");
  //
  tree->Branch("run", &run, "run/I");
  tree->Branch("subrun", &subrun, "subrun/I");
  tree->Branch("eventid", &eventid, "eventid/I");
  //
  tree->Branch("mce_ispr", &mce_ispr, "mce_ispr/I");
  tree->Branch("mce_key", &mce_key, "mce_key/I");
  tree->Branch("mce_nh", &mce_nh, "mce_nh/I");
  tree->Branch("mce_isconv", &mce_isconv, "mce_isconv/I");
  tree->Branch("mce_pdg", &mce_pdg, "mce_pdg/I");
  tree->Branch("mce_mompdg", &mce_mompdg, "mce_mompdg/I");
  tree->Branch("mce_momispr", &mce_momispr, "mce_momispr/I");
  tree->Branch("mce_x", &mce_x, "mce_x/F");
  tree->Branch("mce_y", &mce_y, "mce_y/F");
  tree->Branch("mce_z", &mce_z," mce_z/F");
  tree->Branch("mce_ux", &mce_ux, "mce_ux/F");
  tree->Branch("mce_uy", &mce_uy, "mce_uy/F");
  tree->Branch("mce_uz", &mce_uz, "mce_uz/F");
  tree->Branch("mce_e", &mce_e, "mce_e/F");
  //
  tree->Branch("shw_len", &shw_len, "shw_len/F");
  tree->Branch("shw_nh", &shw_nh, "shw_nh/I");
  tree->Branch("shw_x", &shw_x, "shw_x/F");
  tree->Branch("shw_y", &shw_y, "shw_y/F");
  tree->Branch("shw_z", &shw_z, "shw_z/F");
  tree->Branch("shw_ux", &shw_ux, "shw_ux/F");
  tree->Branch("shw_uy", &shw_uy, "shw_uy/F");
  tree->Branch("shw_uz", &shw_uz, "shw_uz/F");
  tree->Branch("shw_e_p0", &shw_e_p0, "shw_e_p0/F");
  tree->Branch("shw_e_p1", &shw_e_p1, "shw_e_p1/F");
  tree->Branch("shw_e_p2", &shw_e_p2, "shw_e_p2/F");
  tree->Branch("shw_dedx_r4p0", &shw_dedx_r4p0, "shw_shw_dedx_r4p0/F");
  tree->Branch("shw_dedx_r4p1", &shw_dedx_r4p1, "shw_shw_dedx_r4p1/F");
  tree->Branch("shw_dedx_r4p2", &shw_dedx_r4p2, "shw_shw_dedx_r4p2/F");
  tree->Branch("shw_dedx_r2p0", &shw_dedx_r2p0, "shw_shw_dedx_r2p0/F");
  tree->Branch("shw_dedx_r2p1", &shw_dedx_r2p1, "shw_shw_dedx_r2p1/F");
  tree->Branch("shw_dedx_r2p2", &shw_dedx_r2p2, "shw_shw_dedx_r2p2/F");
  tree->Branch("shw_dedx_defp0", &shw_dedx_defp0, "shw_shw_dedx_defp0/F");
  tree->Branch("shw_dedx_defp1", &shw_dedx_defp1, "shw_shw_dedx_defp1/F");
  tree->Branch("shw_dedx_defp2", &shw_dedx_defp2, "shw_shw_dedx_defp2/F");
  tree->Branch("shw_pur", &shw_pur, "shw_pu/F");
  tree->Branch("shw_cmpl", &shw_cmpl, "shw_cmpl/F");
  tree->Branch("shw_mce_dotpdir", &shw_mce_dotpdir, "shw_mce_dotpdir/F");
  tree->Branch("shw_mce_key", &shw_mce_key, "shw_mce_key/I");
  tree->Branch("shw_mce_ispr", &shw_mce_ispr, "shw_mce_ispr/I");
  tree->Branch("shw_mce_isconv", &shw_mce_isconv, "shw_mce_isconv/I");
  tree->Branch("shw_mce_pdg", &shw_mce_pdg, "shw_mce_pdg/I");
  tree->Branch("shw_mce_mompdg", &shw_mce_mompdg, "shw_mce_mompdg/I");
  tree->Branch("shw_mce_momispr", &shw_mce_momispr, "shw_mce_momispr/I");
  //
  tree->Branch("trk_len", &trk_len, "trk_len/F");
  tree->Branch("trk_nh", &trk_nh, "trk_nh/I");
  tree->Branch("trk_nvh", &trk_nvh, "trk_nvh/I");
  tree->Branch("trk_nvh_p0", &trk_nvh_p0, "trk_nvh_p0/I");
  tree->Branch("trk_nvh_p1", &trk_nvh_p1, "trk_nvh_p1/I");
  tree->Branch("trk_nvh_p2", &trk_nvh_p2, "trk_nvh_p2/I");
  tree->Branch("trk_nvh4", &trk_nvh4, "trk_nvh4/I");
  tree->Branch("trk_nvh4_p0", &trk_nvh4_p0, "trk_nvh4_p0/I");
  tree->Branch("trk_nvh4_p1", &trk_nvh4_p1, "trk_nvh4_p1/I");
  tree->Branch("trk_nvh4_p2", &trk_nvh4_p2, "trk_nvh4_p2/I");
  tree->Branch("trk_x", &trk_x, "trk_x/F");
  tree->Branch("trk_y", &trk_y, "trk_y/F");
  tree->Branch("trk_z", &trk_z, "trk_z/F");
  tree->Branch("trk_ux", &trk_ux, "trk_ux/F");
  tree->Branch("trk_uy", &trk_uy, "trk_uy/F");
  tree->Branch("trk_uz", &trk_uz, "trk_uz/F");
  tree->Branch("trk_ex", &trk_ex, "trk_ex/F");
  tree->Branch("trk_ey", &trk_ey, "trk_ey/F");
  tree->Branch("trk_ez", &trk_ez, "trk_ez/F");
  tree->Branch("trk_eux", &trk_eux, "trk_eux/F");
  tree->Branch("trk_euy", &trk_euy, "trk_euy/F");
  tree->Branch("trk_euz", &trk_euz, "trk_euz/F");
  tree->Branch("trk_chi2", &trk_chi2, "trk_chi2/F");
  tree->Branch("trk_ndof", &trk_ndof, "trk_ndof/F");
  tree->Branch("trk_dedx_r4p0", &trk_dedx_r4p0, "trk_trk_dedx_r4p0/F");
  tree->Branch("trk_dedx_r4p1", &trk_dedx_r4p1, "trk_trk_dedx_r4p1/F");
  tree->Branch("trk_dedx_r4p2", &trk_dedx_r4p2, "trk_trk_dedx_r4p2/F");
  tree->Branch("trk_dedx_r2p0", &trk_dedx_r2p0, "trk_trk_dedx_r2p0/F");
  tree->Branch("trk_dedx_r2p1", &trk_dedx_r2p1, "trk_trk_dedx_r2p1/F");
  tree->Branch("trk_dedx_r2p2", &trk_dedx_r2p2, "trk_trk_dedx_r2p2/F");
  tree->Branch("trk_pur", &trk_pur, "trk_pu/F");
  tree->Branch("trk_cmpl", &trk_cmpl, "trk_cmpl/F");
  tree->Branch("trk_mce_dotpdir", &trk_mce_dotpdir, "trk_mce_dotpdir/F");
  tree->Branch("trk_mce_key", &trk_mce_key, "trk_mce_key/I");
  tree->Branch("trk_mce_ispr", &trk_mce_ispr, "trk_mce_ispr/I");
  tree->Branch("trk_mce_isconv", &trk_mce_isconv, "trk_mce_isconv/I");
  tree->Branch("trk_mce_pdg", &trk_mce_pdg, "trk_mce_pdg/I");
  tree->Branch("trk_mce_mompdg", &trk_mce_mompdg, "trk_mce_mompdg/I");
  tree->Branch("trk_mce_momispr", &trk_mce_momispr, "trk_mce_momispr/I");
}

void ShowerPerformance::resetTree() {
  //
  run = -999;
  subrun = -999;
  eventid = -999;
  //
  mce_ispr = -999;
  mce_key = -999;
  mce_isconv = -999;
  mce_pdg = -999;
  mce_mompdg = -999;
  mce_momispr = -999;
  mce_nh = -999;
  mce_x = -999;
  mce_y = -999;
  mce_z = -999;
  mce_ux = -999;
  mce_uy = -999;
  mce_uz = -999;
  mce_e = -999;
  //
  shw_len = -999;
  shw_nh = -999;
  shw_x = -999;
  shw_y = -999;
  shw_z = -999;
  shw_ux = -999;
  shw_uy = -999;
  shw_uz = -999;
  shw_e_p0 = -999;
  shw_e_p1 = -999;
  shw_e_p2 = -999;
  shw_dedx_r4p0 = -999;
  shw_dedx_r4p1 = -999;
  shw_dedx_r4p2 = -999;
  shw_dedx_r2p0 = -999;
  shw_dedx_r2p1 = -999;
  shw_dedx_r2p2 = -999;
  shw_dedx_defp0 = -999;
  shw_dedx_defp1 = -999;
  shw_dedx_defp2 = -999;
  shw_pur = -999;
  shw_cmpl = -999;
  shw_mce_dotpdir = -999;
  shw_mce_key = -999;
  shw_mce_ispr = -999;
  shw_mce_isconv = -999;
  shw_mce_pdg = -999;
  shw_mce_mompdg = -999;
  shw_mce_momispr = -999;
  //
  trk_len = -999;
  trk_nh = -999;
  trk_nvh = -999;
  trk_nvh_p0 = -999;
  trk_nvh_p1 = -999;
  trk_nvh_p2 = -999;
  trk_nvh4 = -999;
  trk_nvh4_p0 = -999;
  trk_nvh4_p1 = -999;
  trk_nvh4_p2 = -999;
  trk_x = -999;
  trk_y = -999;
  trk_z = -999;
  trk_ux = -999;
  trk_uy = -999;
  trk_uz = -999;
  trk_ex = -999;
  trk_ey = -999;
  trk_ez = -999;
  trk_eux = -999;
  trk_euy = -999;
  trk_euz = -999;
  trk_chi2 = -999;
  trk_ndof = -999;
  trk_dedx_r4p0 = -999;
  trk_dedx_r4p1 = -999;
  trk_dedx_r4p2 = -999;
  trk_dedx_r2p0 = -999;
  trk_dedx_r2p1 = -999;
  trk_dedx_r2p2 = -999;
  trk_pur = -999;
  trk_cmpl = -999;
  trk_mce_dotpdir = -999;
  trk_mce_key = -999;
  trk_mce_ispr = -999;
  trk_mce_isconv = -999;
  trk_mce_pdg = -999;
  trk_mce_mompdg = -999;
  trk_mce_momispr = -999;
}

void ShowerPerformance::analyze(art::Event const& e)
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
  auto assocShowers = unique_ptr<art::FindManyP<Shower> >(new art::FindManyP<Shower>(inputPFParticle, e, inputShowerLabel));
  auto const& inputTracks = proxy::getCollection<proxy::Tracks>(e,inputTrackLabel);
  //
  std::unique_ptr<art::FindManyP<recob::Cluster> > assocClusters = std::unique_ptr<art::FindManyP<recob::Cluster> >(new art::FindManyP<recob::Cluster>(inputPFParticle, e, "pandora"));
  auto const& clHitsAssn = *e.getValidHandle<art::Assns<recob::Cluster, recob::Hit> >("pandora");
  const auto& clusterHitsGroups = util::associated_groups(clHitsAssn);
  //
  const auto& inputHits = e.getValidHandle<vector<Hit> >("gaushit");
  const auto& hittruth = unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>(inputHits,e,"gaushitTruthMatch"));
  //
  auto const& mctruth = *e.getValidHandle<std::vector<simb::MCTruth> >("generator");
  auto const& inputMCP = e.getValidHandle<std::vector<simb::MCParticle> >("largeant");

  bool isNuSim = false;

  Point_t nuvtx(0,0,0);
  if (mctruth[0].NeutrinoSet()) {
    isNuSim = true;
    nuvtx = Point_t(mctruth[0].GetNeutrino().Nu().Position().X(),mctruth[0].GetNeutrino().Nu().Position().Y(),mctruth[0].GetNeutrino().Nu().Position().Z());
    if (0) std::cout << "nu vtx=" << nuvtx << " with daughters=" << mctruth[0].GetNeutrino().Nu().NumberDaughters() << std::endl;
    for (int i=0; i<mctruth[0].NParticles(); ++i) {
      if (mctruth[0].GetParticle(i).StatusCode()!=1) continue;
      if (0) cout << "part pdgid=" << mctruth[0].GetParticle(i).PdgCode() << " pos=(" << mctruth[0].GetParticle(i).Vx() << "," << mctruth[0].GetParticle(i).Vy() << "," << mctruth[0].GetParticle(i).Vz() << ") dir=(" << mctruth[0].GetParticle(i).Px()/mctruth[0].GetParticle(i).P() << "," << mctruth[0].GetParticle(i).Py()/mctruth[0].GetParticle(i).P() << "," << mctruth[0].GetParticle(i).Pz()/mctruth[0].GetParticle(i).P() << ") p=" << mctruth[0].GetParticle(i).P() << " status=" << mctruth[0].GetParticle(i).StatusCode() << " process=" << mctruth[0].GetParticle(i).Process() << endl;
    }
  } else {
    nuvtx = Point_t(mctruth[0].GetParticle(0).Position().X(),mctruth[0].GetParticle(0).Position().Y(),mctruth[0].GetParticle(0).Position().Z());
    if (0) cout << "part pdgid=" << mctruth[0].GetParticle(0).PdgCode() << " pos=(" << mctruth[0].GetParticle(0).Vx() << "," << mctruth[0].GetParticle(0).Vy() << "," << mctruth[0].GetParticle(0).Vz() << ") dir=(" << mctruth[0].GetParticle(0).Px()/mctruth[0].GetParticle(0).P() << "," << mctruth[0].GetParticle(0).Py()/mctruth[0].GetParticle(0).P() << "," << mctruth[0].GetParticle(0).Pz()/mctruth[0].GetParticle(0).P() << ") p=" << mctruth[0].GetParticle(0).P() << " status=" << mctruth[0].GetParticle(0).StatusCode() << " process=" << mctruth[0].GetParticle(0).Process() << " id=" << mctruth[0].GetParticle(0).TrackId() << endl;
  }

  //
  std::vector<art::Ptr<recob::Hit> > gaushits;
  for (size_t ih=0;ih<inputHits->size();ih++) gaushits.push_back({inputHits,ih});

  
  //
  art::Ptr<simb::MCParticle> mcp(inputMCP,0);
  if (mctruth[0].NeutrinoSet()) {
    // for (size_t imcp = 0; imcp<inputMCP->size();imcp++) {
    //   art::Ptr<simb::MCParticle> jmcp(inputMCP,imcp);
    //   if (0) cout << "mcp pdgid=" << jmcp->PdgCode()
    // 		  << " pos=(" << jmcp->Vx() << "," << jmcp->Vy() << "," << jmcp->Vz()
    // 		  << ") dir=(" << jmcp->Px()/jmcp->P() << "," << jmcp->Py()/jmcp->P() << "," << jmcp->Pz()/jmcp->P()
    // 		  << ") p=" << jmcp->P() << " status=" << jmcp->StatusCode() << " process=" << jmcp->Process()
    // 		  << " id=" << jmcp->TrackId() << " key=" << mcp.key() << " mom=" << jmcp->Mother()
    // 	       /*<< " firstd=" << jmcp->FirstDaughter() << " lastd=" << jmcp->LastDaughter()*/
    // 		  << std::endl;
    // }
    //
    bool found_bigel = false;
    for (size_t imcp = 0; imcp<inputMCP->size();imcp++) {
      art::Ptr<simb::MCParticle> tmcp(inputMCP,imcp);
      if (std::abs(tmcp->PdgCode())==11 && tmcp->Process()=="primary") {
	mcp = art::Ptr<simb::MCParticle>(inputMCP,imcp);
	found_bigel = true;
	break;
      }
      if (std::abs(tmcp->PdgCode())==111 && tmcp->Process()=="primary") {
	// found pi0, look at daughter gammas
	size_t fd = tmcp->FirstDaughter();
	size_t ld = tmcp->LastDaughter();
	size_t maxkey = 0;
	float maxP = 0;
	for (size_t imcp = 0; imcp<inputMCP->size();imcp++) {
	  art::Ptr<simb::MCParticle> gmcp(inputMCP,imcp);
	  if (std::abs(gmcp->PdgCode())!=22) continue;
	  if (gmcp->TrackId()<fd || gmcp->TrackId()>ld) continue;
	  if (gmcp->P()>maxP) {
	    maxP = gmcp->P();
	    maxkey = gmcp.key();
	  }
	}
	// found gamma, look at daughter electrons
	art::Ptr<simb::MCParticle> gmcp(inputMCP,maxkey);
	if (gmcp->NumberDaughters()<2) continue;
	fd = gmcp->FirstDaughter();
	ld = gmcp->LastDaughter();
	maxkey = 0;
	maxP = 0;
	for (size_t imcp = 0; imcp<inputMCP->size();imcp++) {
	  art::Ptr<simb::MCParticle> emcp(inputMCP,imcp);
	  if (std::abs(emcp->PdgCode())!=11) continue;
	  if (emcp->TrackId()<fd || emcp->TrackId()>ld) continue;
	  if (emcp->P()>maxP) {
	    maxP = emcp->P();
	    maxkey = emcp.key();
	  }
	}
	//
	if (maxkey>0) {
	  mcp = art::Ptr<simb::MCParticle>(inputMCP,maxkey);
	  found_bigel = true;
	  break;
	}
      }
    }
    if (!found_bigel) return;
  } else {
    // for (size_t imcp = 0; imcp<inputMCP->size();imcp++) {
    // 	art::Ptr<simb::MCParticle> mcp(inputMCP,imcp);
    // 	if (0) cout << "mcp pdgid=" << mcp->PdgCode() << " pos=(" << mcp->Vx() << "," << mcp->Vy() << "," << mcp->Vz() << ") dir=(" << mcp->Px()/mcp->P() << "," << mcp->Py()/mcp->P() << "," << mcp->Pz()/mcp->P() << ") p=" << mcp->P() << " status=" << mcp->StatusCode() << " process=" << mcp->Process() << " id=" << mcp->TrackId() << " key=" << mcp.key() << " mom=" << mcp->Mother() << " firstd=" << mcp->FirstDaughter() << " lastd=" << mcp->LastDaughter() << std::endl;	
    // }
    if (std::abs(mcp->PdgCode())==22) {
      size_t fd = mcp->FirstDaughter();
      size_t ld = mcp->LastDaughter();
      size_t maxkey = 0;
      float maxP = 0;
      for (size_t imcp = 0; imcp<inputMCP->size();imcp++) {
	art::Ptr<simb::MCParticle> mcp(inputMCP,imcp);
	if (mcp->TrackId()<fd || mcp->TrackId()>ld) continue;
	if (mcp->P()>maxP) {
	  maxP = mcp->P();
	  maxkey = mcp.key();
	}
      }
      mcp = art::Ptr<simb::MCParticle>(inputMCP,maxkey);
    }
  }
  if (0) cout << "pick mcp pdgid=" << mcp->PdgCode() << " pos=(" << mcp->Vx() << "," << mcp->Vy() << "," << mcp->Vz() << ") dir=(" << mcp->Px()/mcp->P() << "," << mcp->Py()/mcp->P() << "," << mcp->Pz()/mcp->P() << ") p=" << mcp->P() << " status=" << mcp->StatusCode() << " process=" << mcp->Process() << " id=" << mcp->TrackId() << " key=" << mcp.key() << " mom=" << mcp->Mother() << std::endl;
  int nTotMcpHits = nHitsFromMCParticle(mcp.key(),gaushits,hittruth);
  //
  auto scecorr = SCE->GetPosOffsets(geo::Point_t(mcp->Vx(),mcp->Vy(),mcp->Vz()));
  double g4Ticks = detClocks->TPCG4Time2Tick(mcp->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
  double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)+scecorr.X()+0.6;// factor 0.6 to compensate for inter-plane distance
  double yOffset = scecorr.Y();
  double zOffset = scecorr.Z();
  Point_t mcpos( (mcp->Vx()+xOffset)*(1.114/1.098),mcp->Vy()+yOffset,mcp->Vz()+zOffset);//fixme: factor 1.114/1.098 to correct for drift velocity differences
  Vector_t mcdir(mcp->Px()/mcp->P(),mcp->Py()/mcp->P(),mcp->Pz()/mcp->P());

  //
  // Loop over pfps
  for (unsigned int iPfp = 0; iPfp < inputPFParticle->size(); ++iPfp) {
    //
    if (0) std::cout << "ipfp=" << iPfp << std::endl;
    //
    const art::Ptr<recob::PFParticle> pfp(inputPFParticle, iPfp);
    //
    // Showers associated to PFParticles
    const std::vector<art::Ptr<recob::Shower> >& showers = assocShowers->at(iPfp);
    // if there is more than one shower the logic below to get the hits does not work! this works, at least for uboone
    if (0) std::cout << "showers.size()=" << showers.size() << std::endl;
    if (showers.size()!=1) continue;
    //
    // Get hits for shower (through the chain pfp->clusters->hits)
    std::vector<art::Ptr<recob::Hit> > inHits;
    const std::vector<art::Ptr<recob::Cluster> > clustersRange = assocClusters->at(iPfp);
    int maxH = 0;
    int minH = std::numeric_limits<int>::max();
    for (art::Ptr<recob::Cluster> const& cluster: clustersRange) {
      // for hits we use groupByIndex since it preserves the order (and we can use it since each cluster must have associated hits)
      decltype(auto) hitsRange = util::groupByIndex(clusterHitsGroups, cluster.key());
      for (art::Ptr<recob::Hit> const& hit: hitsRange) inHits.push_back(hit);
      int nh = hitsRange.size();
      if (nh<minH) minH = nh;
      if (nh>maxH) maxH = nh;
    }
    if (minH<5 || maxH<10) continue;
    //
    resetTree();
    //
    run = e.run();
    subrun = e.subRun();
    eventid = e.event();
    //
    mce_ispr = (mcp->Process() == "primary");
    mce_key = mcp.key();
    mce_nh = nTotMcpHits;
    mce_x = mcpos.X();
    mce_y = mcpos.Y();
    mce_z = mcpos.Z();
    mce_ux = mcdir.X();
    mce_uy = mcdir.Y();
    mce_uz = mcdir.Z();
    mce_e = mcp->E();
    mce_isconv = (mcp->Process() == "conv");
    mce_pdg = mcp->PdgCode();
    if (mcp->Process() != "primary") {
      for (size_t imcp = 0; imcp<inputMCP->size();imcp++) {
	art::Ptr<simb::MCParticle> mcm(inputMCP,imcp);
	if (mcp->Mother() == mcm->TrackId()) {
	  mce_mompdg = mcm->PdgCode();
	  mce_momispr = (mcm->Process() == "primary");
	  break;
	}
      }
    }
    //
    // Loop over showers (should be only one)
    for (unsigned int iShower = 0; iShower < showers.size(); ++iShower) {
      //
      if (0) std::cout << "iShower=" << iShower << std::endl;
      // Get the shower and convert/hack it into a trajectory so that the fit is initialized
      art::Ptr<recob::Shower> shower = showers[iShower];
      recob::tracking::Point_t pos(shower->ShowerStart().X(),shower->ShowerStart().Y(),shower->ShowerStart().Z());
      recob::tracking::Vector_t dir(shower->Direction().X(),shower->Direction().Y(),shower->Direction().Z());
      //
      art::Ptr<simb::MCParticle> mcp = getAssocMCParticle(inHits,hittruth);
      if (mcp.isNull()) continue;
      //
      if (0) cout << "shw mcp pdgid=" << mcp->PdgCode() << " pos=(" << mcp->Vx() << "," << mcp->Vy() << "," << mcp->Vz() << ") dir=(" << mcp->Px()/mcp->P() << "," << mcp->Py()/mcp->P() << "," << mcp->Pz()/mcp->P() << ") p=" << mcp->P() << " status=" << mcp->StatusCode() << " process=" << mcp->Process() << " id=" << mcp->TrackId() << " key=" << mcp.key() << std::endl;
      //
      int nFoundMcpHits = nHitsFromMCParticle(mcp.key(),inHits,hittruth);
      int nTotMcpHits = nHitsFromMCParticle(mcp.key(),gaushits,hittruth);
      //if (float(nFoundMcpHits)/float(nTotMcpHits)<0.95) continue;//fixme, look only at showers that are complete
      //
      auto scecorr = SCE->GetPosOffsets(geo::Point_t(mcp->Vx(),mcp->Vy(),mcp->Vz()));
      double g4Ticks = detClocks->TPCG4Time2Tick(mcp->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
      double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)+scecorr.X()+0.6;// factor 0.6 to compensate for inter-plane distance
      double yOffset = scecorr.Y();
      double zOffset = scecorr.Z();
      Point_t mcpos( (mcp->Vx()+xOffset)*(1.114/1.098),mcp->Vy()+yOffset,mcp->Vz()+zOffset);//fixme: factor 1.114/1.098 to correct for drift velocity differences
      Vector_t mcdir(mcp->Px()/mcp->P(),mcp->Py()/mcp->P(),mcp->Pz()/mcp->P());
      //
      Length->Fill(shower->Length());
      //
      dotpdir->Fill(dir.Dot(mcdir));
      //
      dx_assoc->Fill(pos.X()-mcpos.X());
      dy_assoc->Fill(pos.Y()-mcpos.Y());
      dz_assoc->Fill(pos.Z()-mcpos.Z());
      dux_assoc->Fill(dir.X()-mcdir.X());
      duy_assoc->Fill(dir.Y()-mcdir.Y());
      duz_assoc->Fill(dir.Z()-mcdir.Z());
      //
      purity->Fill(float(nFoundMcpHits)/float(inHits.size()));
      completeness->Fill(float(nFoundMcpHits)/float(nTotMcpHits));
      //
      auto dedx2v = dEdx(*shower, clustersRange, clusterHitsGroups, 2.0);
      auto dedx4v = dEdx(*shower, clustersRange, clusterHitsGroups, 4.0);
      if (dedx4v.size()==3 && dedx4v[2]>0) dedxpl2rng4->Fill(dedx4v[2]);
      //
      shw_len = shower->Length();
      shw_nh = inHits.size();
      shw_x = pos.X();
      shw_y = pos.Y();
      shw_z = pos.Z();
      shw_ux = dir.X();
      shw_uy = dir.Y();
      shw_uz = dir.Z();
      if (shower->Energy().size()>0) { 
	shw_e_p0 = shower->Energy()[0]/1000.;
	shw_e_p1 = shower->Energy()[1]/1000.;
	shw_e_p2 = shower->Energy()[2]/1000.;
      }
      if (dedx4v.size()>0) {
	shw_dedx_r4p0 = dedx4v[0];
	shw_dedx_r4p1 = dedx4v[1];
	shw_dedx_r4p2 = dedx4v[2];
      }
      if (dedx2v.size()>0) {
	shw_dedx_r2p0 = dedx2v[0];
	shw_dedx_r2p1 = dedx2v[1];
	shw_dedx_r2p2 = dedx2v[2];
      }
      if (shower->dEdx().size()>0) {
	shw_dedx_defp0 = shower->dEdx()[0];
	shw_dedx_defp1 = shower->dEdx()[1];
	shw_dedx_defp2 = shower->dEdx()[2];
      }
      shw_pur = float(nFoundMcpHits)/float(inHits.size());
      shw_cmpl = float(nFoundMcpHits)/float(nTotMcpHits);
      shw_mce_dotpdir = dir.Dot(mcdir);
      //
      shw_mce_key = mcp.key();
      shw_mce_ispr = (mcp->Process() == "primary");
      //
      shw_mce_isconv = (mcp->Process() == "conv");
      shw_mce_pdg = mcp->PdgCode();
      if (mcp->Process() != "primary") {
	for (size_t imcp = 0; imcp<inputMCP->size();imcp++) {
	  art::Ptr<simb::MCParticle> mcm(inputMCP,imcp);
	  if (mcp->Mother() == mcm->TrackId()) {
	    shw_mce_mompdg = mcm->PdgCode();
	    shw_mce_momispr = (mcm->Process() == "primary");
	    break;
	  }
	}
      }
      //
    }// iShower
    //
    //now figure out if there is a track associated to the same pfp
    for (auto t : inputTracks) {
      if (t->ID()!=iPfp) continue;
      if (0) std::cout << "found track" << std::endl;
      recob::tracking::Point_t tk_pos = t->Start();
      recob::tracking::Vector_t tk_dir = t->StartDirection();
      //
      std::vector<art::Ptr<recob::Hit> > validtrkhits;
      size_t ih=0;
      size_t vh_p[] = {0,0,0};
      size_t vh4_p[] = {0,0,0};
      for (auto h : t.hits()) {
	if (t->HasValidPoint(ih)) {
	  validtrkhits.push_back(h);
	  vh_p[t.hits()[ih]->WireID().Plane]++;
	  if ((t->Length()-t->Length(ih))<4) vh4_p[t.hits()[ih]->WireID().Plane]++;
	}
	ih++;
      }
      //
      art::Ptr<simb::MCParticle> tk_mcp = getAssocMCParticle(validtrkhits,hittruth);
      if (tk_mcp.isNull()) continue;
      //
      auto scecorr = SCE->GetPosOffsets(geo::Point_t(tk_mcp->Vx(),tk_mcp->Vy(),tk_mcp->Vz()));
      double g4Ticks = detClocks->TPCG4Time2Tick(tk_mcp->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
      double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)+scecorr.X()+0.6;// factor 0.6 to compensate for inter-plane distance
      double yOffset = scecorr.Y();
      double zOffset = scecorr.Z();
      Point_t mcpos( (tk_mcp->Vx()+xOffset)*(1.114/1.098),tk_mcp->Vy()+yOffset,tk_mcp->Vz()+zOffset);//fixme: factor 1.114/1.098 to correct for drift velocity differences
      Vector_t mcdir(tk_mcp->Px()/tk_mcp->P(),tk_mcp->Py()/tk_mcp->P(),tk_mcp->Pz()/tk_mcp->P());
      //
      int tk_nFoundMcpHits = nHitsFromMCParticle(tk_mcp.key(),validtrkhits,hittruth);
      int tk_nTotMcpHits = nHitsFromMCParticle(tk_mcp.key(),gaushits,hittruth);
      //
      if (0) cout << "trk tk_mcp pdgid=" << tk_mcp->PdgCode() << " pos=(" << tk_mcp->Vx() << "," << tk_mcp->Vy() << "," << tk_mcp->Vz() << ") dir=(" << tk_mcp->Px()/tk_mcp->P() << "," << tk_mcp->Py()/tk_mcp->P() << "," << tk_mcp->Pz()/tk_mcp->P() << ") p=" << tk_mcp->P() << " status=" << tk_mcp->StatusCode() << " process=" << tk_mcp->Process() << " id=" << tk_mcp->TrackId() << " key=" << tk_mcp.key() << endl;
      //
      // if (float(tk_nFoundMcpHits)/float(validtrkhits.size())<0.1) continue;//fixme
      //
      tk_Length->Fill(t->Length());
      //
      tk_dotpdir->Fill(tk_dir.Dot(mcdir));
      //
      tk_dx_assoc->Fill(tk_pos.X()-mcpos.X());
      tk_dy_assoc->Fill(tk_pos.Y()-mcpos.Y());
      tk_dz_assoc->Fill(tk_pos.Z()-mcpos.Z());
      tk_dux_assoc->Fill(tk_dir.X()-mcdir.X());
      tk_duy_assoc->Fill(tk_dir.Y()-mcdir.Y());
      tk_duz_assoc->Fill(tk_dir.Z()-mcdir.Z());
      //
      //auto c = t->VertexCovarianceGlobal6D();
      auto c = tracking::Plane::Local5DToGlobal6DCovariance(t->VertexCovarianceLocal5D(),false,t->VertexMomentumVector(),t->VertexDirection());
      tk_dx_assoc_pull->Fill( (tk_pos.X()-mcpos.X())/sqrt(c(0,0)) );
      tk_dy_assoc_pull->Fill( (tk_pos.Y()-mcpos.Y())/sqrt(c(1,1)) );
      tk_dz_assoc_pull->Fill( (tk_pos.Z()-mcpos.Z())/sqrt(c(2,2)) );
      tk_dux_assoc_pull->Fill( (tk_dir.X()-mcdir.X())/sqrt(c(3,3)) );
      tk_duy_assoc_pull->Fill( (tk_dir.Y()-mcdir.Y())/sqrt(c(4,4)) );
      tk_duz_assoc_pull->Fill( (tk_dir.Z()-mcdir.Z())/sqrt(c(5,5)) );
      //
      tk_purity->Fill(float(tk_nFoundMcpHits)/float(validtrkhits.size()));
      tk_completeness->Fill(float(tk_nFoundMcpHits)/float(tk_nTotMcpHits));
      //
      auto dedx4v = dEdx(*t, t.hits(), 4.0);
      auto dedx2v = dEdx(*t, t.hits(), 2.0);
      if (dedx4v.size()==3 && dedx4v[2]>0) tk_dedxpl2rng4->Fill(dedx4v[2]);
      //
      trk_len = t->Length();
      trk_nh = t->NumberTrajectoryPoints();
      trk_nvh = validtrkhits.size();
      trk_nvh_p0 = vh_p[0];
      trk_nvh_p1 = vh_p[1];
      trk_nvh_p2 = vh_p[2];
      trk_nvh4 = vh4_p[0]+vh4_p[1]+vh4_p[2];
      trk_nvh4_p0 = vh4_p[0];
      trk_nvh4_p1 = vh4_p[1];
      trk_nvh4_p2 = vh4_p[2];
      trk_x = tk_pos.X();
      trk_y = tk_pos.Y();
      trk_z = tk_pos.Z();
      trk_ux = tk_dir.X();
      trk_uy = tk_dir.Y();
      trk_uz = tk_dir.Z();
      trk_ex = sqrt(c(0,0));
      trk_ey = sqrt(c(1,1));
      trk_ez = sqrt(c(2,2));
      trk_eux = sqrt(c(3,3));
      trk_euy = sqrt(c(4,4));
      trk_euz = sqrt(c(5,5));
      trk_chi2 = t->Chi2();
      trk_ndof = t->Ndof();
      if (dedx4v.size()>0) {
	trk_dedx_r4p0 = dedx4v[0];
	trk_dedx_r4p1 = dedx4v[1];
	trk_dedx_r4p2 = dedx4v[2];
      }
      if (dedx2v.size()>0) {
	trk_dedx_r2p0 = dedx2v[0];
	trk_dedx_r2p1 = dedx2v[1];
	trk_dedx_r2p2 = dedx2v[2];
      }
      trk_pur = float(tk_nFoundMcpHits)/float(validtrkhits.size());
      trk_cmpl = float(tk_nFoundMcpHits)/float(tk_nTotMcpHits);
      trk_mce_dotpdir = tk_dir.Dot(mcdir);
      //
      trk_mce_key = tk_mcp.key();
      trk_mce_ispr = (tk_mcp->Process() == "primary");
      //
      trk_mce_isconv = (tk_mcp->Process() == "conv");
      trk_mce_pdg = tk_mcp->PdgCode();
      if (tk_mcp->Process() != "primary") {
	for (size_t imcp = 0; imcp<inputMCP->size();imcp++) {
	  art::Ptr<simb::MCParticle> mcm(inputMCP,imcp);
	  if (tk_mcp->Mother() == mcm->TrackId()) {
	    trk_mce_mompdg = mcm->PdgCode();
	    trk_mce_momispr = (mcm->Process() == "primary");
	    break;
	  }
	}
      }
      //
    }// inputTracks
    //
    tree->Fill();
    //
  }
  
}


art::Ptr<simb::MCParticle> ShowerPerformance::getAssocMCParticle(const std::vector<art::Ptr<recob::Hit> >& hits,
								 const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& hittruth) const {
  //credit: Wes Ketchum
  std::unordered_map<int,double> trkide;
  double maxe=-1, tote=0;
  art::Ptr<simb::MCParticle> maxp_me; //pointer for the particle match we will calculate
  for (auto h : hits) {
    std::vector<art::Ptr<simb::MCParticle> > particle_vec = hittruth->at(h.key());
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec = hittruth->data(h.key());;
    //loop over particles
    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
      trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
      tote += match_vec[i_p]->energy; //calculate total energy deposited
      if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
	maxe = trkide[ particle_vec[i_p]->TrackId() ];
	maxp_me = particle_vec[i_p];
      }
    }//end loop over particles per hit
  }
  return maxp_me;
}

int ShowerPerformance::nHitsFromMCParticle(size_t mcid,
					   const std::vector<art::Ptr<recob::Hit> >& hits,
					   const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& hittruth) const {

  int count = 0;
  for (auto h : hits) {
    std::vector<art::Ptr<simb::MCParticle> > particle_vec = hittruth->at(h.key());
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec = hittruth->data(h.key());;
    //loop over particles
    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
      if (match_vec[i_p]->isMaxIDE==0) continue;
      if (particle_vec[i_p].key()!=mcid) continue;
      count++;
    }
  }
  return count;
}

template <typename T>
std::vector<float> ShowerPerformance::dEdx(const recob::Shower& shower, const std::vector<art::Ptr<recob::Cluster> >& clusters, const T& clusterHitsGroups, float range) {

  //std::cout << "dEdx begin" << std::endl;
  //std::cout << "clusters.size()=" << clusters.size() << std::endl;

  std::vector<float> result(3,0);
  if (clusters.size()==0) return result;

  // grab shower direction
  auto const& dir3D = shower.Direction();

  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  float wire2cm = geom->WirePitch(0,0,0);
  float time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
  
  // loop through planes
  for (size_t n = 0; n < clusters.size(); n++) {

    auto const& clus = clusters.at(n);
    //std::cout << "cluster n=" << n << " key=" << clus.key() << std::endl;

    // get the hits associated with this cluster
    decltype(auto) hitsRange = util::groupByIndex(clusterHitsGroups, clus.key());

    // get the plane associated with this cluster
    auto const& pl = clus->Plane();      
    const geo::WireGeo& wire = geom->TPC().Plane(pl).MiddleWire();
    TVector3 wireunitperp = wire.Direction();//(wire.GetStart()-wire.GetEnd()).Unit();
    // rotate by 90 degrees around x
    TVector3 wireunit = {wireunitperp[0], -wireunitperp[2], wireunitperp[1]}; 
    float cosPlane = fabs(cos(wireunit.Angle(dir3D)));

    std::vector<float> dedx_v;
    float pitch = 0.3 / cosPlane;

    //std::cout << "nhits=" << hitsRange.size() << std::endl;
    
    // loop over hits
    for (art::Ptr<recob::Hit> const& h : hitsRange) {
      //
      double d2D = sqrt( pow(h->WireID().Wire*wire2cm - clus->StartWire()*wire2cm, 2) + pow(h->PeakTime()*time2cm - clus->StartTick()*time2cm, 2) );
      //std::cout << "w1=" << h->WireID().Wire*wire2cm << " w2=" << clus->StartWire()*wire2cm << " t1=" << h->PeakTime()*time2cm << " t2=" << clus->StartTick()*time2cm << std::endl;
      double d3D = d2D/cosPlane;
      //std::cout << "d2D=" << d2D << " d3D=" << d3D << " range=" << range << std::endl;
      if (d3D>range) continue;
      //
      float dE = h->Integral() * 248.2 * 0.0000236 / 0.6;//fixme// * ChargeCorrection(h.charge, h.w, h.t, resultShower.fDCosStart, resultShower.fXYZStart, pl, energyCalibProvider);
      float dEdx = dE / pitch;
      //std::cout << "integral=" << h->Integral() << " pitch=" << pitch << " dEdx=" << dEdx << std::endl;
      dedx_v.push_back(dEdx);
    }// loop over all hits
    if (dedx_v.size()>2) {
      std::sort( dedx_v.begin(), dedx_v.end() );
      //std::cout << "median dedx=" << dedx_v[dedx_v.size()/2] << std::endl;
      result[pl.Plane] = dedx_v[dedx_v.size()/2];
    }
  }// for all clusters (planes)

  return result;
}

template <typename T>
std::vector<float> ShowerPerformance::dEdx(const recob::Track& track, const T& hits, float range) {

  std::vector<float> result(3,0);
  if (hits.size()==0) return result;

  auto const* geom = ::lar::providerFrom<geo::Geometry>();

  // loop through hits
  std::vector<std::vector<float> > dedx_v(3);
  for (size_t n = 0; n < track.CountValidPoints(); n++) {

    auto const& hit = hits.at(n);

    // get the plane associated with this cluster
    auto const& pl = hit->WireID().planeID();
      
    const geo::WireGeo& wire = geom->TPC().Plane(pl).MiddleWire();
    TVector3 wireunitperp = wire.Direction();//(wire.GetStart()-wire.GetEnd()).Unit();
    // rotate by 90 degrees around x
    TVector3 wireunit = {wireunitperp[0], -wireunitperp[2], wireunitperp[1]}; 

    // grab track direction
    auto const& dir3D = track.DirectionAtPoint<TVector3>(n);
    float cosPlane = fabs(cos(wireunit.Angle(dir3D)));

    float pitch = 0.3 / cosPlane;

    //
    //std::cout << "n=" << n << " pl=" << pl.Plane << " range=" << track.Length()-track.Length(n) << " l0=" << track.Length() << " ln=" << track.Length(n) << std::endl;
    //
    if (track.Length()-track.Length(n)>range) continue;
    //
    float dE = hit->Integral() * 248.2 * 0.0000236 / 0.6;//fixme// * ChargeCorrection(h.charge, h.w, h.t, resultShower.fDCosStart, resultShower.fXYZStart, pl, energyCalibProvider);
    float dEdx = dE / pitch;
    dedx_v[pl.Plane].push_back(dEdx);
  }// loop over all hits

  for (size_t np=0; np<dedx_v.size(); np++) {
    if (dedx_v[np].size()>2) {
      std::sort( dedx_v[np].begin(), dedx_v[np].end() );
      result[np] = dedx_v[np][dedx_v[np].size()/2];
    }
  }

  return result;
}

DEFINE_ART_MODULE(ShowerPerformance)
