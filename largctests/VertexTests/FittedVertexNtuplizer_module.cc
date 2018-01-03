////////////////////////////////////////////////////////////////////////
// Class:       FittedVertexNtuplizer
// Plugin Type: analyzer (art v2_05_00)
// File:        FittedVertexNtuplizer_module.cc
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
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingPlane.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "lardata/RecoObjects/TrackingPlaneHelper.h"

#include "larreco/RecoAlg/Geometric3DVertexFitter.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

class FittedVertexNtuplizer;

class FittedVertexNtuplizer : public art::EDAnalyzer {
public:
  struct Inputs {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<art::InputTag> inputPFLabel {
      Name("inputPFLabel"),
	Comment("Label of recob::PFParticle Collection to be fit")
	};
    fhicl::Atom<art::InputTag> fittedTracksLabel {
      Name("fittedTracksLabel"),
	Comment("Label of fitted recob::Track Collection associated to PFParticles")
	};
  };
  
  struct Config {
    using Name = fhicl::Name;
    fhicl::Table<FittedVertexNtuplizer::Inputs> inputs {
      Name("inputs"),
	};
    fhicl::Table<trkf::Geometric3DVertexFitter::Config> options {
      Name("options")
	};
    fhicl::Table<trkf::TrackStatePropagator::Config> propagator {
      Name("propagator")
	};
  };
  using Parameters = art::EDAnalyzer::Table<Config>;
  
  explicit FittedVertexNtuplizer(Parameters const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FittedVertexNtuplizer(FittedVertexNtuplizer const &) = delete;
  FittedVertexNtuplizer(FittedVertexNtuplizer &&) = delete;
  FittedVertexNtuplizer & operator = (FittedVertexNtuplizer const &) = delete;
  FittedVertexNtuplizer & operator = (FittedVertexNtuplizer &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  void beginJob() override;
  void resetTree();
private:
  //
  art::InputTag PFInputTag;
  art::InputTag fittedTrackInputTag;
  trkf::Geometric3DVertexFitter fitter;
  //
  TTree* tree;
  //
  // Event
  //
  int run, subrun, eventid;
  int passSelII;
  //
  // MC Vertex
  //
  float vtx_gen_x;
  float vtx_gen_y;
  float vtx_gen_z;
  int vtx_gen_ntks;
  int gen_nuflav;
  float gen_numom;
  float gen_lepmom;
  //
  // Pandora Vertex
  //
  float vtx_pan_x;
  float vtx_pan_y;
  float vtx_pan_z;
  int vtx_pan_ntks;
  std::vector<int> vtx_pan_tknh;
  std::vector<float> vtx_pan_tklen;
  std::vector<float> vtx_pan_tkvx;
  std::vector<float> vtx_pan_tkvy;
  std::vector<float> vtx_pan_tkvz;
  std::vector<float> vtx_pan_tkvux;
  std::vector<float> vtx_pan_tkvuy;
  std::vector<float> vtx_pan_tkvuz;
  std::vector<float> vtx_pan_tkex;
  std::vector<float> vtx_pan_tkey;
  std::vector<float> vtx_pan_tkez;
  std::vector<float> vtx_pan_tkeux;
  std::vector<float> vtx_pan_tkeuy;
  std::vector<float> vtx_pan_tkeuz;
  int vtx_pan_tknh_max;
  int vtx_pan_tknh_min;
  float vtx_pan_tklen_max;
  float vtx_pan_tklen_min;
  //
  // Fitted Vertex
  //
  float vtx_fit_x;
  float vtx_fit_y;
  float vtx_fit_z;
  int vtx_fit_ntks;
  float vtx_fit_cxx;
  float vtx_fit_cyy;
  float vtx_fit_czz;
  float vtx_fit_cxy;
  float vtx_fit_cxz;
  float vtx_fit_cyz;
  float vtx_fit_chi2;
  int vtx_fit_ndof;
  std::vector<float> vtx_fit_tkpids;
  std::vector<float> vtx_fit_tkdists;
  std::vector<int> vtx_fit_tknh;
  std::vector<float> vtx_fit_tklen;
  std::vector<float> vtx_fit_tkvx;
  std::vector<float> vtx_fit_tkvy;
  std::vector<float> vtx_fit_tkvz;
  std::vector<float> vtx_fit_tkvux;
  std::vector<float> vtx_fit_tkvuy;
  std::vector<float> vtx_fit_tkvuz;
  std::vector<float> vtx_fit_tkex;
  std::vector<float> vtx_fit_tkey;
  std::vector<float> vtx_fit_tkez;
  std::vector<float> vtx_fit_tkeux;
  std::vector<float> vtx_fit_tkeuy;
  std::vector<float> vtx_fit_tkeuz;
  std::vector<float> vtx_fit_tkchi2;
  std::vector<float> vtx_fit_tkndof;
  std::vector<float> vtx_fit_tkubchi2;
  std::vector<float> vtx_fit_tkubip;
  std::vector<float> vtx_fit_tkubsip;
  int vtx_fit_tknh_max;
  int vtx_fit_tknh_min;
  float vtx_fit_tklen_max;
  float vtx_fit_tklen_min;
  float vtx_fit_tkdists_max;
  float vtx_fit_tkdists_min;
  float vtx_fit_tkubchi2_max;
  float vtx_fit_tkubchi2_min;
  float vtx_fit_tkubip_max;
  float vtx_fit_tkubip_min;
  float vtx_fit_tkubsip_max;
  float vtx_fit_tkubsip_min;
};

void FittedVertexNtuplizer::beginJob()
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
  // MC Vertex
  //
  tree->Branch("vtx_gen_x", &vtx_gen_x,"vtx_gen_x/F");
  tree->Branch("vtx_gen_y", &vtx_gen_y,"vtx_gen_y/F");
  tree->Branch("vtx_gen_z", &vtx_gen_z,"vtx_gen_z/F");
  tree->Branch("vtx_gen_ntks", &vtx_gen_ntks,"vtx_gen_ntks/I");
  tree->Branch("gen_nuflav", &gen_nuflav,"gen_nuflav/I");
  tree->Branch("gen_numom", &gen_numom,"gen_numom/F");
  tree->Branch("gen_lepmom", &gen_lepmom,"gen_lepmom/F");
  //
  // Pandora Vertex
  //
  tree->Branch("vtx_pan_x", &vtx_pan_x,"vtx_pan_x/F");
  tree->Branch("vtx_pan_y", &vtx_pan_y,"vtx_pan_y/F");
  tree->Branch("vtx_pan_z", &vtx_pan_z,"vtx_pan_z/F");
  tree->Branch("vtx_pan_ntks", &vtx_pan_ntks,"vtx_pan_ntks/I");
  tree->Branch("vtx_pan_tknh"  , &vtx_pan_tknh );
  tree->Branch("vtx_pan_tklen" , &vtx_pan_tklen);
  tree->Branch("vtx_pan_tkvx"  , &vtx_pan_tkvx );
  tree->Branch("vtx_pan_tkvy"  , &vtx_pan_tkvy );
  tree->Branch("vtx_pan_tkvz"  , &vtx_pan_tkvz );
  tree->Branch("vtx_pan_tkvux" , &vtx_pan_tkvux);
  tree->Branch("vtx_pan_tkvuy" , &vtx_pan_tkvuy);
  tree->Branch("vtx_pan_tkvuz" , &vtx_pan_tkvuz);
  tree->Branch("vtx_pan_tkex"  , &vtx_pan_tkex );
  tree->Branch("vtx_pan_tkey"  , &vtx_pan_tkey );
  tree->Branch("vtx_pan_tkez"  , &vtx_pan_tkez );
  tree->Branch("vtx_pan_tkeux" , &vtx_pan_tkeux);
  tree->Branch("vtx_pan_tkeuy" , &vtx_pan_tkeuy);
  tree->Branch("vtx_pan_tkeuz" , &vtx_pan_tkeuz);
  tree->Branch("vtx_pan_tknh_max", &vtx_pan_tknh_max,"vtx_pan_tknh_max/I");
  tree->Branch("vtx_pan_tknh_min", &vtx_pan_tknh_min,"vtx_pan_tknh_min/I");
  tree->Branch("vtx_pan_tklen_max", &vtx_pan_tklen_max,"vtx_pan_tklen_max/F");
  tree->Branch("vtx_pan_tklen_min", &vtx_pan_tklen_min,"vtx_pan_tklen_min/F");
  //
  // Fitted Vertex
  //
  tree->Branch("vtx_fit_x", &vtx_fit_x,"vtx_fit_x/F");
  tree->Branch("vtx_fit_y", &vtx_fit_y,"vtx_fit_y/F");
  tree->Branch("vtx_fit_z", &vtx_fit_z,"vtx_fit_z/F");
  tree->Branch("vtx_fit_ntks", &vtx_fit_ntks,"vtx_fit_ntks/I");
  tree->Branch("vtx_fit_tknh"  , &vtx_fit_tknh );
  tree->Branch("vtx_fit_tklen" , &vtx_fit_tklen);
  tree->Branch("vtx_fit_tkvx"  , &vtx_fit_tkvx );
  tree->Branch("vtx_fit_tkvy"  , &vtx_fit_tkvy );
  tree->Branch("vtx_fit_tkvz"  , &vtx_fit_tkvz );
  tree->Branch("vtx_fit_tkvux" , &vtx_fit_tkvux);
  tree->Branch("vtx_fit_tkvuy" , &vtx_fit_tkvuy);
  tree->Branch("vtx_fit_tkvuz" , &vtx_fit_tkvuz);
  tree->Branch("vtx_fit_tkex"  , &vtx_fit_tkex );
  tree->Branch("vtx_fit_tkey"  , &vtx_fit_tkey );
  tree->Branch("vtx_fit_tkez"  , &vtx_fit_tkez );
  tree->Branch("vtx_fit_tkeux" , &vtx_fit_tkeux);
  tree->Branch("vtx_fit_tkeuy" , &vtx_fit_tkeuy);
  tree->Branch("vtx_fit_tkeuz" , &vtx_fit_tkeuz);
  tree->Branch("vtx_fit_cxx", &vtx_fit_cxx,"vtx_fit_cxx/F");
  tree->Branch("vtx_fit_cyy", &vtx_fit_cyy,"vtx_fit_cyy/F");
  tree->Branch("vtx_fit_czz", &vtx_fit_czz,"vtx_fit_czz/F");
  tree->Branch("vtx_fit_cxy", &vtx_fit_cxy,"vtx_fit_cxy/F");
  tree->Branch("vtx_fit_cxz", &vtx_fit_cxz,"vtx_fit_cxz/F");
  tree->Branch("vtx_fit_cyz", &vtx_fit_cyz,"vtx_fit_cyz/F");
  tree->Branch("vtx_fit_chi2", &vtx_fit_chi2,"vtx_fit_chi2/F");
  tree->Branch("vtx_fit_ndof", &vtx_fit_ndof,"vtx_fit_ndof/I");
  tree->Branch("vtx_fit_tkpids" , &vtx_fit_tkpids);
  tree->Branch("vtx_fit_tkdists", &vtx_fit_tkdists);
  tree->Branch("vtx_fit_tkchi2", &vtx_fit_tkchi2);
  tree->Branch("vtx_fit_tkndof", &vtx_fit_tkndof);
  tree->Branch("vtx_fit_tkubchi2", &vtx_fit_tkubchi2);
  tree->Branch("vtx_fit_tkubip", &vtx_fit_tkubip);
  tree->Branch("vtx_fit_tkubsip", &vtx_fit_tkubsip);
  tree->Branch("vtx_fit_tknh_max", &vtx_fit_tknh_max,"vtx_fit_tknh_max/I");
  tree->Branch("vtx_fit_tknh_min", &vtx_fit_tknh_min,"vtx_fit_tknh_min/I");
  tree->Branch("vtx_fit_tklen_max", &vtx_fit_tklen_max,"vtx_fit_tklen_max/F");
  tree->Branch("vtx_fit_tklen_min", &vtx_fit_tklen_min,"vtx_fit_tklen_min/F");
  tree->Branch("vtx_fit_tkdists_max", &vtx_fit_tkdists_max,"vtx_fit_tkdists_max/F");
  tree->Branch("vtx_fit_tkdists_min", &vtx_fit_tkdists_min,"vtx_fit_tkdists_min/F");
  tree->Branch("vtx_fit_tkubchi2_max", &vtx_fit_tkubchi2_max,"vtx_fit_tkubchi2_max/F");
  tree->Branch("vtx_fit_tkubchi2_min", &vtx_fit_tkubchi2_min,"vtx_fit_tkubchi2_min/F");
  tree->Branch("vtx_fit_tkubip_max", &vtx_fit_tkubip_max,"vtx_fit_tkubip_max/F");
  tree->Branch("vtx_fit_tkubip_min", &vtx_fit_tkubip_min,"vtx_fit_tkubip_min/F");
  tree->Branch("vtx_fit_tkubsip_max", &vtx_fit_tkubsip_max,"vtx_fit_tkubsip_max/F");
  tree->Branch("vtx_fit_tkubsip_min", &vtx_fit_tkubsip_min,"vtx_fit_tkubsip_min/F");
}

FittedVertexNtuplizer::FittedVertexNtuplizer(Parameters const & p)
  : EDAnalyzer(p)
  , PFInputTag(p().inputs().inputPFLabel())
  , fittedTrackInputTag(p().inputs().fittedTracksLabel())
  , fitter(p().options,p().propagator)
{}

void FittedVertexNtuplizer::resetTree() {
  //
  // Event
  //
  run = -999;
  subrun = -999;
  eventid = -999;
  passSelII = -999;
  //
  // MC Vertex
  //
  vtx_gen_x = -999;
  vtx_gen_y = -999;
  vtx_gen_z = -999;
  vtx_gen_ntks = -999;
  gen_nuflav = -999;
  gen_numom = -999;
  gen_lepmom = -999;
  //
  // Pandora Vertex
  //
  vtx_pan_x = -999;
  vtx_pan_y = -999;
  vtx_pan_z = -999;
  vtx_pan_ntks = -999;
  vtx_pan_tknh .clear();
  vtx_pan_tklen.clear();
  vtx_pan_tkvx .clear();
  vtx_pan_tkvy .clear();
  vtx_pan_tkvz .clear();
  vtx_pan_tkvux.clear();
  vtx_pan_tkvuy.clear();
  vtx_pan_tkvuz.clear();
  vtx_pan_tkex .clear();
  vtx_pan_tkey .clear();
  vtx_pan_tkez .clear();
  vtx_pan_tkeux.clear();
  vtx_pan_tkeuy.clear();
  vtx_pan_tkeuz.clear();
  vtx_pan_tknh_max = -999;
  vtx_pan_tknh_min = -999;
  vtx_pan_tklen_max = -999;
  vtx_pan_tklen_min = -999;
  //
  // Fitted Vertex
  //
  vtx_fit_x = -999;
  vtx_fit_y = -999;
  vtx_fit_z = -999;
  vtx_fit_ntks = -999;
  vtx_fit_tknh .clear();
  vtx_fit_tklen.clear();
  vtx_fit_tkvx .clear();
  vtx_fit_tkvy .clear();
  vtx_fit_tkvz .clear();
  vtx_fit_tkvux.clear();
  vtx_fit_tkvuy.clear();
  vtx_fit_tkvuz.clear();
  vtx_fit_tkex .clear();
  vtx_fit_tkey .clear();
  vtx_fit_tkez .clear();
  vtx_fit_tkeux.clear();
  vtx_fit_tkeuy.clear();
  vtx_fit_tkeuz.clear();
  vtx_fit_cxx = -999;
  vtx_fit_cyy = -999;
  vtx_fit_czz = -999;
  vtx_fit_cxy = -999;
  vtx_fit_cxz = -999;
  vtx_fit_cyz = -999;
  vtx_fit_chi2 = -999;
  vtx_fit_ndof = -999;
  vtx_fit_tkpids.clear();
  vtx_fit_tkdists.clear();
  vtx_fit_tkchi2.clear();
  vtx_fit_tkndof.clear();
  vtx_fit_tkubchi2.clear();
  vtx_fit_tkubip.clear();
  vtx_fit_tkubsip.clear();
  vtx_fit_tknh_max = -999;
  vtx_fit_tknh_min = -999;
  vtx_fit_tklen_max = -999;
  vtx_fit_tklen_min = -999;
  vtx_fit_tkdists_max = -999;
  vtx_fit_tkdists_min = -999;
  vtx_fit_tkubchi2_max = -999;
  vtx_fit_tkubchi2_min = -999;
  vtx_fit_tkubip_max = -999;
  vtx_fit_tkubip_min = -999;
  vtx_fit_tkubsip_max = -999;
  vtx_fit_tkubsip_min = -999;
}

void FittedVertexNtuplizer::analyze(art::Event const & e)
{
  //
  using namespace std;
  using namespace trkf;
  using namespace recob::tracking;
  //
  art::ValidHandle<art::TriggerResults> filter = e.getValidHandle<art::TriggerResults>("TriggerResults");
  size_t ntp =  art::ServiceHandle<art::TriggerNamesService>()->size();
  size_t ftp = ntp;
  for (size_t itp=0;itp<ntp;itp++) {
    //std::cout << art::ServiceHandle<art::TriggerNamesService>()->getTrigPath(itp) << " " << filter->at(itp).accept()  << std::endl;
    if (art::ServiceHandle<art::TriggerNamesService>()->getTrigPath(itp)=="filt") ftp = itp; 
  }
  assert(ftp<ntp);
  //
  detinfo::DetectorProperties const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
  detinfo::DetectorClocks const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>(); 
  TrackStatePropagator prop(1.0,0.1,10,10.,0.01,false);
  //
  art::ValidHandle<std::vector<recob::PFParticle> > inputPFParticle = e.getValidHandle<std::vector<recob::PFParticle> >(PFInputTag);
  auto assocTracks = std::unique_ptr<art::FindManyP<recob::Track> >(new art::FindManyP<recob::Track>(inputPFParticle, e, PFInputTag));
  auto assocVertices = std::unique_ptr<art::FindManyP<recob::Vertex> >(new art::FindManyP<recob::Vertex>(inputPFParticle, e, PFInputTag));
  //
  auto assocFittedTracks = std::unique_ptr<art::FindManyP<recob::Track> >(new art::FindManyP<recob::Track>(inputPFParticle, e, fittedTrackInputTag));
  //
  art::InputTag TruthInputTag("generator");
  const std::vector<simb::MCTruth>* mcTruth = 0;
  if (e.isRealData()==0) {
    mcTruth = e.getValidHandle<std::vector<simb::MCTruth> >(TruthInputTag).product();
  }
  //
  for (unsigned int iPF = 0; iPF < inputPFParticle->size(); ++iPF) {
    //
    art::Ptr<recob::PFParticle> pfp(inputPFParticle, iPF);
    if (pfp->IsPrimary()==false || pfp->NumDaughters()<2) continue;
    //
    vector< art::Ptr<recob::Track> > pandora_tracks;
    vector< art::Ptr<recob::Track> > fitted_tracks;
    auto& pfd = pfp->Daughters();
    for (auto ipfd : pfd) {
      art::Ptr<recob::PFParticle> pfpd(inputPFParticle, ipfd);
      //
      // pandora tracks
      vector< art::Ptr<recob::Track> > patracks = assocTracks->at(ipfd);
      for (auto t : patracks) {
	pandora_tracks.push_back(t);
      }
      // kalmantrack tracks
      vector< art::Ptr<recob::Track> > kftracks = assocFittedTracks->at(ipfd);
      for (auto t : kftracks) {
	fitted_tracks.push_back(t);
      }
    }
    if (pandora_tracks.size()<2) continue;
    //
    resetTree();
    //
    run = e.run();
    subrun = e.subRun();
    eventid = e.event();
    passSelII = filter->at(ftp).accept();
    //
    const std::vector<art::Ptr<recob::Vertex> >& pandoraVertices = assocVertices->at(iPF);
    for (unsigned int iVertex = 0; iVertex < pandoraVertices.size(); ++iVertex) {
      art::Ptr<recob::Vertex> pvertex = pandoraVertices[iVertex];
      double xyz[3];
      pvertex->XYZ(xyz);
      vtx_pan_x = xyz[0];
      vtx_pan_y = xyz[1];
      vtx_pan_z = xyz[2];
      break;
    }
    if (0) cout << "pandora vertex pos=" << Point_t(vtx_pan_x,vtx_pan_y,vtx_pan_z) << " ntracks=" <<  pandora_tracks.size() << endl;
    vtx_pan_ntks = pandora_tracks.size();
    for (size_t itk=0;itk<pandora_tracks.size();++itk) {
      const auto* tk = pandora_tracks[itk].get();
      vtx_pan_tknh .push_back(tk->CountValidPoints());
      vtx_pan_tklen.push_back(tk->Length());
      vtx_pan_tkvx .push_back(tk->Start().X());
      vtx_pan_tkvy .push_back(tk->Start().Y());
      vtx_pan_tkvz .push_back(tk->Start().Z());
      vtx_pan_tkvux.push_back(tk->StartDirection().X());
      vtx_pan_tkvuy.push_back(tk->StartDirection().Y());
      vtx_pan_tkvuz.push_back(tk->StartDirection().Z());
      vtx_pan_tkex .push_back(tk->End().X());
      vtx_pan_tkey .push_back(tk->End().Y());
      vtx_pan_tkez .push_back(tk->End().Z());
      vtx_pan_tkeux.push_back(tk->EndDirection().X());
      vtx_pan_tkeuy.push_back(tk->EndDirection().Y());
      vtx_pan_tkeuz.push_back(tk->EndDirection().Z());
    }
    auto tmp_vtx_pan_tknh = vtx_pan_tknh;
    sort(tmp_vtx_pan_tknh.begin(),tmp_vtx_pan_tknh.end());
    vtx_pan_tknh_max = tmp_vtx_pan_tknh.back();
    vtx_pan_tknh_min = tmp_vtx_pan_tknh.front();
    auto tmp_vtx_pan_tklen = vtx_pan_tklen;
    sort(tmp_vtx_pan_tklen.begin(),tmp_vtx_pan_tklen.end());
    vtx_pan_tklen_max = tmp_vtx_pan_tklen.back();
    vtx_pan_tklen_min = tmp_vtx_pan_tklen.front();
    //
    trkf::VertexWrapper fitvtx = (fitted_tracks.size()>1 ? fitter.fitTracks(fitted_tracks) : trkf::VertexWrapper());
    if (fitvtx.isValid()) {
      vtx_fit_x = fitvtx.position().X();
      vtx_fit_y = fitvtx.position().Y();
      vtx_fit_z = fitvtx.position().Z();
      if (0) cout << "fitted vertex pos=" << recob::tracking::Point_t(vtx_fit_x,vtx_fit_y,vtx_fit_z) << " ntracks=" <<  fitvtx.tracks().size() << endl;
      vtx_fit_ntks = fitvtx.tracks().size();
      vtx_fit_cxx = fitvtx.covariance()(0,0);
      vtx_fit_cyy = fitvtx.covariance()(1,1);
      vtx_fit_czz = fitvtx.covariance()(2,2);
      vtx_fit_cxy = fitvtx.covariance()(0,1);
      vtx_fit_cxz = fitvtx.covariance()(0,2);
      vtx_fit_cyz = fitvtx.covariance()(1,2);
      vtx_fit_chi2 = fitvtx.vertex().chi2();
      vtx_fit_ndof = fitvtx.vertex().ndof();
      for (size_t itk=0;itk<fitvtx.tracks().size();++itk) {
	const auto tk = fitvtx.tracks()[itk];
	vtx_fit_tkpids.push_back(tk.get().ParticleId());
	vtx_fit_tknh .push_back(tk.get().CountValidPoints());
	vtx_fit_tklen.push_back(tk.get().Length());
	vtx_fit_tkvx .push_back(tk.get().Start().X());
	vtx_fit_tkvy .push_back(tk.get().Start().Y());
	vtx_fit_tkvz .push_back(tk.get().Start().Z());
	vtx_fit_tkvux.push_back(tk.get().StartDirection().X());
	vtx_fit_tkvuy.push_back(tk.get().StartDirection().Y());
	vtx_fit_tkvuz.push_back(tk.get().StartDirection().Z());
	vtx_fit_tkex .push_back(tk.get().End().X());
	vtx_fit_tkey .push_back(tk.get().End().Y());
	vtx_fit_tkez .push_back(tk.get().End().Z());
	vtx_fit_tkeux.push_back(tk.get().EndDirection().X());
	vtx_fit_tkeuy.push_back(tk.get().EndDirection().Y());
	vtx_fit_tkeuz.push_back(tk.get().EndDirection().Z());
	vtx_fit_tkchi2.push_back(tk.get().Chi2());
	vtx_fit_tkndof.push_back(tk.get().Ndof());
	trkf::VertexWrapper ubvtx = fitter.unbiasedVertex(fitvtx, tk);
	vtx_fit_tkdists.push_back( fitter.pDistUnbiased(ubvtx, tk) );
	vtx_fit_tkubchi2.push_back( fitter.chi2Unbiased(ubvtx, tk) );
	vtx_fit_tkubip.push_back( fitter.ipUnbiased(ubvtx, tk) );
	vtx_fit_tkubsip.push_back( fitter.sipUnbiased(ubvtx, tk) );
      }
      auto tmp_vtx_fit_tknh = vtx_fit_tknh;
      sort(tmp_vtx_fit_tknh.begin(),tmp_vtx_fit_tknh.end());
      vtx_fit_tknh_max = tmp_vtx_fit_tknh.back();
      vtx_fit_tknh_min = tmp_vtx_fit_tknh.front();
      auto tmp_vtx_fit_tklen = vtx_fit_tklen;
      sort(tmp_vtx_fit_tklen.begin(),tmp_vtx_fit_tklen.end());
      vtx_fit_tklen_max = tmp_vtx_fit_tklen.back();
      vtx_fit_tklen_min = tmp_vtx_fit_tklen.front();
      auto tmp_vtx_fit_tkdists = vtx_fit_tkdists;
      sort(tmp_vtx_fit_tkdists.begin(),tmp_vtx_fit_tkdists.end());
      vtx_fit_tkdists_max = tmp_vtx_fit_tkdists.back();
      vtx_fit_tkdists_min = tmp_vtx_fit_tkdists.front();
      auto tmp_vtx_fit_tkubchi2 = vtx_fit_tkubchi2;
      sort(tmp_vtx_fit_tkubchi2.begin(),tmp_vtx_fit_tkubchi2.end());
      vtx_fit_tkubchi2_max = tmp_vtx_fit_tkubchi2.back();
      vtx_fit_tkubchi2_min = tmp_vtx_fit_tkubchi2.front();
      auto tmp_vtx_fit_tkubip = vtx_fit_tkubip;
      sort(tmp_vtx_fit_tkubip.begin(),tmp_vtx_fit_tkubip.end());
      vtx_fit_tkubip_max = tmp_vtx_fit_tkubip.back();
      vtx_fit_tkubip_min = tmp_vtx_fit_tkubip.front();
      auto tmp_vtx_fit_tkubsip = vtx_fit_tkubsip;
      sort(tmp_vtx_fit_tkubsip.begin(),tmp_vtx_fit_tkubsip.end());
      vtx_fit_tkubsip_max = tmp_vtx_fit_tkubsip.back();
      vtx_fit_tkubsip_min = tmp_vtx_fit_tkubsip.front();
    }
    //
    if (e.isRealData()==0) {
      const auto& mct = mcTruth->at(0).GetNeutrino();
      auto scecorr = SCE->GetPosOffsets(geo::Point_t(mct.Lepton().Vx(),mct.Lepton().Vy(),mct.Lepton().Vz()));
      double ticks = detClocks->TPCG4Time2Tick(mct.Lepton().T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
      double xOffset = theDetector->ConvertTicksToX(ticks, 0, 0, 0)-scecorr.X();
      double yOffset = scecorr.Y();
      double zOffset = scecorr.Z();
      vtx_gen_x = mct.Lepton().Vx()+xOffset;
      vtx_gen_y = mct.Lepton().Vy()+yOffset;
      vtx_gen_z = mct.Lepton().Vz()+zOffset;
      vtx_gen_ntks = mct.Nu().NumberDaughters();
      gen_nuflav = mct.Nu().PdgCode();
      gen_numom = mct.Nu().P();
      gen_lepmom = mct.Lepton().P();
    }
    //
    tree->Fill();
  }
}

DEFINE_ART_MODULE(FittedVertexNtuplizer)
