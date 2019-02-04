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
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardata/RecoObjects/TrackStatePropagator.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/Utilities/ForEachAssociatedGroup.h"

#include "TTree.h"

// std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;

class TrackingPerformance : public art::EDAnalyzer {
public:
  explicit TrackingPerformance(fhicl::ParameterSet const & p);

  TrackingPerformance(TrackingPerformance const &) = delete;
  TrackingPerformance(TrackingPerformance &&) = delete;
  TrackingPerformance & operator = (TrackingPerformance const &) = delete;
  TrackingPerformance & operator = (TrackingPerformance &&) = delete;

  void analyze(art::Event const & e) override;

  void beginJob() override;

art::Ptr<simb::MCParticle> getAssocMCParticle(const std::vector<art::Ptr<recob::Hit> >&,
                                              const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& hittruth) const;

  int nHitsFromMCParticle(size_t mcid,
                          const std::vector<art::Ptr<recob::Hit> >&,
                          const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& hittruth) const;

  void resetTree();

private:
  //
  std::string inputTracksLabel;
  //
  TTree* tree;
  //
  // event id
  int run, subrun, eventid;
  //
  // generator information
  int gen_nupdg, gen_mode, gen_ccnc, gen_inttype;
  float gen_nue, gen_nux, gen_nuy, gen_nuz, gen_nupx, gen_nupy, gen_nupz,  gen_nux_sccorr, gen_nuy_sccorr, gen_nuz_sccorr;
  //
  // mcparticle information
  int mcp_num;
  std::vector<int> mcp_idx, mcp_pdg, mcp_isprimary;
  std::vector<float> mcp_e, mcp_ke, mcp_p, mcp_start_x, mcp_start_y, mcp_start_z, mcp_start_px, mcp_start_py, mcp_start_pz;
  std::vector<float> mcp_start_x_sccorr, mcp_start_y_sccorr, mcp_start_z_sccorr;
  std::vector<float> mcp_end_x, mcp_end_y, mcp_end_z;
  //
  // pandora vertex information
  int vtx_nupdg, vtx_nprimaries, vtx_ntrks;
  float vtx_x, vtx_y, vtx_z;
  std::vector<int> vtx_pdgprimaries;
  //
  // pandora track information
  std::vector<int> trk_idx, trk_nhits, trk_nhits_p0, trk_nhits_p1, trk_nhits_p2;
  std::vector<float> trk_start_x, trk_start_y, trk_start_z, trk_start_ux, trk_start_uy, trk_start_uz;
  std::vector<float> trk_end_x, trk_end_y, trk_end_z, trk_end_ux, trk_end_uy, trk_end_uz;
  std::vector<float> trk_length;
  std::vector<int> trk_mcp_idx, trk_mcp_ntoth, trk_mcp_ntrkh;
  std::vector<float> trk_mcp_length;
  //
  //
};


TrackingPerformance::TrackingPerformance(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
  , inputTracksLabel(p.get<std::string>("inputTracksLabel"))
{}

void TrackingPerformance::resetTree() {
  //
  // event id
  run = -999;
  subrun = -999;
  eventid = -999;
  //
  // generator information
  gen_nupdg = -999;
  gen_mode = -999;
  gen_ccnc = -999;
  gen_inttype = -999;
  gen_nue = -999;
  gen_nux = -999;
  gen_nuy = -999;
  gen_nuz = -999;
  gen_nupx = -999;
  gen_nupy = -999;
  gen_nupz = -999;
  gen_nux_sccorr = -999;
  gen_nuy_sccorr = -999;
  gen_nuz_sccorr = -999;
  //
  // mcparticle information
  mcp_num = -999;
  mcp_idx.clear();
  mcp_pdg.clear();
  mcp_isprimary.clear();
  mcp_e.clear();
  mcp_ke.clear();
  mcp_p.clear();
  mcp_start_x.clear();
  mcp_start_y.clear();
  mcp_start_z.clear();
  mcp_start_px.clear();
  mcp_start_py.clear();
  mcp_start_pz.clear();
  mcp_start_x_sccorr.clear();
  mcp_start_y_sccorr.clear();
  mcp_start_z_sccorr.clear();
  mcp_end_x.clear();
  mcp_end_y.clear();
  mcp_end_z.clear();
  //
  // pandora vertex information
  vtx_nupdg = -999;
  vtx_nprimaries = -999;
  vtx_ntrks = -999;
  vtx_x = -999;
  vtx_y = -999;
  vtx_z = -999;
  vtx_pdgprimaries.clear();
  //
  // pandora track information
  trk_idx.clear();
  trk_nhits.clear();
  trk_nhits_p0.clear();
  trk_nhits_p1.clear();
  trk_nhits_p2.clear();
  trk_start_x.clear();
  trk_start_y.clear();
  trk_start_z.clear();
  trk_start_ux.clear();
  trk_start_uy.clear();
  trk_start_uz.clear();
  trk_end_x.clear();
  trk_end_y.clear();
  trk_end_z.clear();
  trk_end_ux.clear();
  trk_end_uy.clear();
  trk_end_uz.clear();
  trk_length.clear();
  trk_mcp_idx.clear();
  trk_mcp_ntoth.clear();
  trk_mcp_ntrkh.clear();
  trk_mcp_length.clear();
}

void TrackingPerformance::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  //
  tree = tfs->make<TTree>("tree", "tree");
  //
  // event id
  tree->Branch("run",&run,"run/I");
  tree->Branch("subrun",&subrun,"subrun/I");
  tree->Branch("eventid",&eventid,"eventid/I");
  //
  // generator information
  tree->Branch("gen_nupdg",&gen_nupdg,"gen_nupdg/I");
  tree->Branch("gen_mode",&gen_mode,"gen_mode/I");
  tree->Branch("gen_ccnc",&gen_ccnc,"gen_ccnc/I");
  tree->Branch("gen_inttype",&gen_inttype,"gen_inttype/I");
  tree->Branch("gen_nue",&gen_nue,"gen_nue/F");
  tree->Branch("gen_nux",&gen_nux,"gen_nux/F");
  tree->Branch("gen_nuy",&gen_nuy,"gen_nuy/F");
  tree->Branch("gen_nuz",&gen_nuz,"gen_nuz/F");
  tree->Branch("gen_nupx",&gen_nupx,"gen_nupx/F");
  tree->Branch("gen_nupy",&gen_nupy,"gen_nupy/F");
  tree->Branch("gen_nupz",&gen_nupz,"gen_nupz/F");
  tree->Branch("gen_nux_sccorr",&gen_nux_sccorr,"gen_nux_sccorr/F");
  tree->Branch("gen_nuy_sccorr",&gen_nuy_sccorr,"gen_nuy_sccorr/F");
  tree->Branch("gen_nuz_sccorr",&gen_nuz_sccorr,"gen_nuz_sccorr/F");
  //
  // mcparticle information
  tree->Branch("mcp_num",&mcp_num,"mcp_num/I");
  tree->Branch("mcp_idx",&mcp_idx);
  tree->Branch("mcp_pdg",&mcp_pdg);
  tree->Branch("mcp_isprimary",&mcp_isprimary);
  tree->Branch("mcp_e",&mcp_e);
  tree->Branch("mcp_ke",&mcp_ke);
  tree->Branch("mcp_p",&mcp_p);
  tree->Branch("mcp_start_x",&mcp_start_x);
  tree->Branch("mcp_start_y",&mcp_start_y);
  tree->Branch("mcp_start_z",&mcp_start_z);
  tree->Branch("mcp_start_px",&mcp_start_px);
  tree->Branch("mcp_start_py",&mcp_start_py);
  tree->Branch("mcp_start_pz",&mcp_start_pz);
  tree->Branch("mcp_start_x_sccorr",&mcp_start_x_sccorr);
  tree->Branch("mcp_start_y_sccorr",&mcp_start_y_sccorr);
  tree->Branch("mcp_start_z_sccorr",&mcp_start_z_sccorr);
  tree->Branch("mcp_end_x",&mcp_end_x);
  tree->Branch("mcp_end_y",&mcp_end_y);
  tree->Branch("mcp_end_z",&mcp_end_z);
  //
  // pandora vertex information
  tree->Branch("vtx_nupdg",&vtx_nupdg,"vtx_nupdg/I");
  tree->Branch("vtx_nprimaries",&vtx_nprimaries,"vtx_nprimaries/I");
  tree->Branch("vtx_ntrks",&vtx_ntrks,"vtx_ntrks/I");
  tree->Branch("vtx_x",&vtx_x,"vtx_x/F");
  tree->Branch("vtx_y",&vtx_y,"vtx_y/F");
  tree->Branch("vtx_z",&vtx_z,"vtx_z/F");
  tree->Branch("vtx_pdgprimaries",&vtx_pdgprimaries);
  //
  // pandora track information
  tree->Branch("trk_idx",&trk_idx);
  tree->Branch("trk_nhits",&trk_nhits);
  tree->Branch("trk_nhits_p0",&trk_nhits_p0);
  tree->Branch("trk_nhits_p1",&trk_nhits_p1);
  tree->Branch("trk_nhits_p2",&trk_nhits_p2);
  tree->Branch("trk_start_x",&trk_start_x);
  tree->Branch("trk_start_y",&trk_start_y);
  tree->Branch("trk_start_z",&trk_start_z);
  tree->Branch("trk_start_ux",&trk_start_ux);
  tree->Branch("trk_start_uy",&trk_start_uy);
  tree->Branch("trk_start_uz",&trk_start_uz);
  tree->Branch("trk_end_x",&trk_end_x);
  tree->Branch("trk_end_y",&trk_end_y);
  tree->Branch("trk_end_z",&trk_end_z);
  tree->Branch("trk_end_ux",&trk_end_ux);
  tree->Branch("trk_end_uy",&trk_end_uy);
  tree->Branch("trk_end_uz",&trk_end_uz);
  tree->Branch("trk_length",&trk_length);
  tree->Branch("trk_mcp_idx",&trk_mcp_idx);
  tree->Branch("trk_mcp_ntoth",&trk_mcp_ntoth);
  tree->Branch("trk_mcp_ntrkh",&trk_mcp_ntrkh);
  tree->Branch("trk_mcp_length",&trk_mcp_length);
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
  const auto& inputPFParticle = e.getValidHandle<vector<PFParticle> >("pandoraNu");
  const auto& inputTracks = e.getValidHandle<vector<Track> >(inputTracksLabel);
  //
  const auto& assocVertex = unique_ptr<art::FindManyP<Vertex> >(new art::FindManyP<Vertex>(inputPFParticle, e, inputTracksLabel));
  //
  const auto& assocTracks = unique_ptr<art::FindManyP<Track> >(new art::FindManyP<Track>(inputPFParticle, e, inputTracksLabel));
  const auto& tkHitsAssn = unique_ptr<art::FindManyP<Hit> >(new art::FindManyP<Hit>(inputTracks, e, inputTracksLabel));
  //
  const auto& inputHits = e.getValidHandle<vector<Hit> >("gaushit");
  const auto& hittruth = unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>(inputHits,e,"gaushitTruthMatch"));
  //
  auto const& mctruth = *e.getValidHandle<std::vector<simb::MCTruth> >("generator");
  auto const& mcparts = *e.getValidHandle<std::vector<simb::MCParticle> >("largeant");
  //

  Point_t nuvtx(mctruth[0].GetNeutrino().Nu().Position().X(),mctruth[0].GetNeutrino().Nu().Position().Y(),mctruth[0].GetNeutrino().Nu().Position().Z());
  if (1) std::cout << "nu vtx=" << nuvtx
                   << " mode=" << mctruth[0].GetNeutrino().Mode()
                   << " CCNC=" << mctruth[0].GetNeutrino().CCNC()
                   << " int_type=" << mctruth[0].GetNeutrino().InteractionType()
                   << " Qsqr=" << mctruth[0].GetNeutrino().QSqr()
                   << " with daughters=" << mctruth[0].GetNeutrino().Nu().NumberDaughters() << std::endl;
  // for (int i=0; i<mctruth[0].NParticles(); ++i) {
  //   if (mctruth[0].GetParticle(i).StatusCode()!=1) continue;
  //   if (1) cout << "part pdgid=" << mctruth[0].GetParticle(i).PdgCode() << " id=" << mctruth[0].GetParticle(i).TrackId() << " pos=(" << mctruth[0].GetParticle(i).Vx() << "," << mctruth[0].GetParticle(i).Vy() << "," << mctruth[0].GetParticle(i).Vz() << ") dir=(" << mctruth[0].GetParticle(i).Px()/mctruth[0].GetParticle(i).P() << "," << mctruth[0].GetParticle(i).Py()/mctruth[0].GetParticle(i).P() << "," << mctruth[0].GetParticle(i).Pz()/mctruth[0].GetParticle(i).P() << ") p=" << mctruth[0].GetParticle(i).P() << " status=" << mctruth[0].GetParticle(i).StatusCode() << " process=" << mctruth[0].GetParticle(i).Process() << endl;
  // }
  for (auto m : mcparts) {
    if (m.P()<0.100) continue;
    if (1) cout << "mcpart pdgid=" << m.PdgCode() << " id=" << m.TrackId() << " pos=(" << m.Vx() << "," << m.Vy() << "," << m.Vz() << ") dir=(" << m.Px()/m.P() << "," << m.Py()/m.P() << "," << m.Pz()/m.P() << ") p=" << m.P() << " ke=" << m.E()-m.Mass() << " status=" << m.StatusCode() << " process=" << m.Process() << endl;
  }

  for (size_t iPF = 0; iPF < inputPFParticle->size(); ++iPF) {

    art::Ptr<PFParticle> pfp(inputPFParticle, iPF);
    if (pfp->IsPrimary()==false) continue;
    //
    resetTree();
    //
    run = e.run();
    subrun = e.subRun();
    eventid = e.event();
    //
    //
    auto scecorr_nu = SCE->GetPosOffsets(mctruth[0].GetNeutrino().Nu().Vx(),mctruth[0].GetNeutrino().Nu().Vy(),mctruth[0].GetNeutrino().Nu().Vz());
    double g4Ticks_nu = detClocks->TPCG4Time2Tick(mctruth[0].GetNeutrino().Nu().T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
    double xOffset_nu = theDetector->ConvertTicksToX(g4Ticks_nu, 0, 0, 0)-scecorr_nu[0];
    double yOffset_nu = scecorr_nu[1];
    double zOffset_nu = scecorr_nu[2];
    //
    gen_nupdg = mctruth[0].GetNeutrino().Nu().PdgCode();
    gen_mode = mctruth[0].GetNeutrino().Mode();
    gen_ccnc = mctruth[0].GetNeutrino().CCNC();
    gen_inttype = mctruth[0].GetNeutrino().InteractionType();
    gen_nue = mctruth[0].GetNeutrino().Nu().E();
    gen_nux = mctruth[0].GetNeutrino().Nu().Vx();
    gen_nuy = mctruth[0].GetNeutrino().Nu().Vy();
    gen_nuz = mctruth[0].GetNeutrino().Nu().Vz();
    gen_nupx = mctruth[0].GetNeutrino().Nu().Px();
    gen_nupy = mctruth[0].GetNeutrino().Nu().Py();
    gen_nupz = mctruth[0].GetNeutrino().Nu().Pz();
    gen_nux_sccorr = mctruth[0].GetNeutrino().Nu().Vx()+xOffset_nu;
    gen_nuy_sccorr = mctruth[0].GetNeutrino().Nu().Vy()+yOffset_nu;
    gen_nuz_sccorr = mctruth[0].GetNeutrino().Nu().Vz()+zOffset_nu;
    //
    // mcparticle information
    mcp_num = 0;
    for (auto m : mcparts) {
      if (abs(m.PdgCode())==13 && (m.E()-m.Mass())<0.035) continue;
      if (abs(m.PdgCode())==2212 && (m.E()-m.Mass())<0.060) continue;
      if (abs(m.PdgCode())!=13 && abs(m.PdgCode())!=2212 && m.P()<0.100) continue;
      auto scecorr_mcp = SCE->GetPosOffsets(m.Vx(),m.Vy(),m.Vz());
      double g4Ticks_mcp = detClocks->TPCG4Time2Tick(m.T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
      double xOffset_mcp = theDetector->ConvertTicksToX(g4Ticks_mcp, 0, 0, 0)-scecorr_mcp[0];
      double yOffset_mcp = scecorr_mcp[1];
      double zOffset_mcp = scecorr_mcp[2];
      mcp_num++;
      mcp_idx.push_back(m.TrackId());
      mcp_pdg.push_back(m.PdgCode());
      mcp_isprimary.push_back(m.Process()=="primary");
      mcp_e.push_back(m.E());
      mcp_ke.push_back(m.E()-m.Mass());
      mcp_p.push_back(m.P());
      mcp_start_x.push_back(m.Vx());
      mcp_start_y.push_back(m.Vy());
      mcp_start_z.push_back(m.Vz());
      mcp_start_px.push_back(m.Px());
      mcp_start_py.push_back(m.Py());
      mcp_start_pz.push_back(m.Pz());
      mcp_start_x_sccorr.push_back(m.Vx()+xOffset_mcp);
      mcp_start_y_sccorr.push_back(m.Vy()+yOffset_mcp);
      mcp_start_z_sccorr.push_back(m.Vz()+zOffset_mcp);
      mcp_end_x.push_back(m.EndX());
      mcp_end_y.push_back(m.EndY());
      mcp_end_z.push_back(m.EndZ());
    }
    //
    auto vtx = assocVertex->at(iPF);
    double xyz[3];
    vtx[0]->XYZ(xyz);
    if (1) std::cout << "pfp#" << iPF << " PdgCode=" << pfp->PdgCode()
              << " IsPrimary=" << pfp->IsPrimary()
              << " NumDaughters=" << pfp->NumDaughters()
              << " vtx=" << Point_t(xyz[0],xyz[1],xyz[2])
              << std::endl;
    auto& pfd = pfp->Daughters();

    vtx_nupdg = pfp->PdgCode();
    vtx_nprimaries = pfp->NumDaughters();
    vtx_x = xyz[0];
    vtx_y = xyz[1];
    vtx_z = xyz[2];
    vtx_ntrks = 0;
    for (auto ipfd : pfd) {
      vector< art::Ptr<Track> > pftracks = assocTracks->at(ipfd);
      art::Ptr<PFParticle> pfpd(inputPFParticle, ipfd);
      if (1) cout << "pfp id=" << ipfd << " pdg=" << pfpd->PdgCode() << endl;
      vtx_pdgprimaries.push_back(pfpd->PdgCode());
      vtx_ntrks+=pftracks.size();
    }

    for (auto ipfd : pfd) {
      vector< art::Ptr<Track> > pftracks = assocTracks->at(ipfd);
      art::Ptr<PFParticle> pfpd(inputPFParticle, ipfd);
      if (1) cout << "pfp id=" << ipfd << " pdg=" << pfpd->PdgCode() << endl;

      for (auto t : pftracks) {
        //
        if (1) std::cout << "track id=" << t->ID() << " nhits=" << t->NumberTrajectoryPoints()
                         << " len=" << t->Length()
                         << " start=" << t->Trajectory().Vertex() << " end=" << t->Trajectory().End()
                         << " startdir=" << t->Trajectory().VertexDirection() << " enddir=" << t->Trajectory().EndDirection()
                         << std::endl;
        //
        int nhp0 = 0, nhp1 = 0, nhp2 = 0;
        std::vector<art::Ptr<recob::Hit> > alltrkhits;
        auto hitsRange = tkHitsAssn->at(t.key());
        for (art::Ptr<recob::Hit> const& hit: hitsRange) {
          if (hit->WireID().Plane==0) nhp0++;
          if (hit->WireID().Plane==1) nhp1++;
          if (hit->WireID().Plane==2) nhp2++;
          //find match in gaushit (needed because the mc truth is broken for the cosmic-removed hit collection)
          for (unsigned int ig=0; ig<inputHits->size(); ig++) {
            const auto& gahit = inputHits->at(ig);
            if (hit->Channel()==gahit.Channel() && std::abs(hit->PeakTime()-gahit.PeakTime())<0.000001) {
              alltrkhits.push_back(art::Ptr<recob::Hit>(inputHits,ig));
              break;
            }
          }
        }

        //
        trk_idx.push_back(t->ID());
        trk_nhits.push_back(alltrkhits.size());
        trk_nhits_p0.push_back(nhp0);
        trk_nhits_p1.push_back(nhp1);
        trk_nhits_p2.push_back(nhp2);
      	trk_start_x.push_back(t->Start().X());
	trk_start_y.push_back(t->Start().Y());
	trk_start_z.push_back(t->Start().Z());
	trk_start_ux.push_back(t->StartDirection().X());
	trk_start_uy.push_back(t->StartDirection().Y());
	trk_start_uz.push_back(t->StartDirection().Z());
      	trk_end_x.push_back(t->End().X());
	trk_end_y.push_back(t->End().Y());
	trk_end_z.push_back(t->End().Z());
	trk_end_ux.push_back(t->EndDirection().X());
	trk_end_uy.push_back(t->EndDirection().Y());
	trk_end_uz.push_back(t->EndDirection().Z());
	trk_length.push_back(t->Length());
	//

        art::Ptr<simb::MCParticle> mcp = getAssocMCParticle(alltrkhits,hittruth);
        if (mcp.isNull()) continue;

        int nFoundMcpHits = nHitsFromMCParticle(mcp.key(),alltrkhits,hittruth);
        std::vector<art::Ptr<recob::Hit> > gaushits;
        for (size_t ih=0;ih<inputHits->size();ih++) gaushits.push_back({inputHits,ih});
        int nTotMcpHits = nHitsFromMCParticle(mcp.key(),gaushits,hittruth);

        auto mcl = mcp->Trajectory().TotalLength();
        //
        auto scecorr = SCE->GetPosOffsets(mcp->Vx(),mcp->Vy(),mcp->Vz());
        double g4Ticks = detClocks->TPCG4Time2Tick(mcp->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
        double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr[0];
        double yOffset = scecorr[1];
        double zOffset = scecorr[2];
        Point_t mcpos(mcp->Vx()+xOffset,mcp->Vy()+yOffset,mcp->Vz()+zOffset);
        Vector_t mcdir(mcp->Px()/mcp->P(),mcp->Py()/mcp->P(),mcp->Pz()/mcp->P());
        //
        if (1) cout << "matched particle with id=" << mcp->TrackId() << " nFoundMcpHits=" << nFoundMcpHits << " nTotMcpHits=" << nTotMcpHits
		    << " status=" << mcp->StatusCode()
		    << " length=" << mcl
		    << " pos=(" << mcp->Vx() << "," << mcp->Vy() << "," << mcp->Vz() 
		    << ") dir=(" << mcp->Px()/mcp->P() << "," << mcp->Py()/mcp->P() << "," << mcp->Pz()/mcp->P() << ") p=" << mcp->P() 
		    << " sce_corr_pos=" << mcpos
		    << " pdg=" << mcp->PdgCode() << " mom=" << mcp->Mother() << " proc=" << mcp->Process() 
		    << endl;
        if (1) cout << "scecorr=" << scecorr[0] << ", " << scecorr[1] << ", " << scecorr[2] << endl;
	//
	trk_mcp_idx.push_back(mcp->TrackId());
	trk_mcp_ntoth.push_back(nTotMcpHits);
	trk_mcp_ntrkh.push_back(nFoundMcpHits);
	trk_mcp_length.push_back(mcl);
      }
    }
    tree->Fill();
  }


}

art::Ptr<simb::MCParticle> TrackingPerformance::getAssocMCParticle(const std::vector<art::Ptr<recob::Hit> >& hits,
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

int TrackingPerformance::nHitsFromMCParticle(size_t mcid,
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

DEFINE_ART_MODULE(TrackingPerformance)
