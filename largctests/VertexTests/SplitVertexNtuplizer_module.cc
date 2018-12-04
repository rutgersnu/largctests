////////////////////////////////////////////////////////////////////////
// Class:       SplitVertexNtuplizer
// Plugin Type: analyzer (art v2_05_00)
// File:        SplitVertexNtuplizer_module.cc
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
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/VertexAssnMeta.h"
#include "lardataobj/RecoBase/Track.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

class SplitVertexNtuplizer;

class SplitVertexNtuplizer : public art::EDAnalyzer {
public:
  explicit SplitVertexNtuplizer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SplitVertexNtuplizer(SplitVertexNtuplizer const &) = delete;
  SplitVertexNtuplizer(SplitVertexNtuplizer &&) = delete;
  SplitVertexNtuplizer & operator = (SplitVertexNtuplizer const &) = delete;
  SplitVertexNtuplizer & operator = (SplitVertexNtuplizer &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  void beginJob() override;
  void resetTree();
private:
  //
  art::InputTag inputPFLabel;
  art::InputTag inputTracksLabel;
  art::InputTag inputVertexLabel;
  art::InputTag inputVertexLabel1;
  art::InputTag inputVertexLabel2;
  //
  TTree* tree;
  //
  // Event
  //
  int  run, subrun, eventid;
  int  passSelII;
  //
  // Original vertex
  //
  float vchi2;
  int   vstatus;
  int   vndof;
  std::vector<float> vpos;
  std::vector<float> vcov;
  std::vector<int>   vtkey;
  std::vector<float> vtpd;
  std::vector<float> vtip;
  std::vector<float> vtipe;
  std::vector<float> vtsip;
  std::vector<float> vtchi2;
  std::vector<int>   vtstatus;
  std::vector<float> tk_vtx_x;
  std::vector<float> tk_vtx_y;
  std::vector<float> tk_vtx_z;
  std::vector<float> tk_vtx_ux;
  std::vector<float> tk_vtx_uy;
  std::vector<float> tk_vtx_uz;
  std::vector<float> tk_vtx_xe;
  std::vector<float> tk_vtx_ye;
  std::vector<float> tk_vtx_ze;
  std::vector<float> tk_vtx_uxe;
  std::vector<float> tk_vtx_uye;
  std::vector<float> tk_vtx_uze;
  std::vector<float> tk_end_x;
  std::vector<float> tk_end_y;
  std::vector<float> tk_end_z;
  std::vector<float> tk_end_ux;
  std::vector<float> tk_end_uy;
  std::vector<float> tk_end_uz;
  std::vector<float> tk_end_xe;
  std::vector<float> tk_end_ye;
  std::vector<float> tk_end_ze;
  std::vector<float> tk_end_uxe;
  std::vector<float> tk_end_uye;
  std::vector<float> tk_end_uze;
  std::vector<int>   tnpoints;
  std::vector<int>   tnvalid;
  std::vector<float> tlength;
  std::vector<float> tchi2;
  std::vector<int>   tndof;
  std::vector<int>   tpdg;
  std::vector<int>   tcontain;
  std::vector<float> trvtx;
  std::vector<float> tmom;
  //
  // Split vertex 1
  //
  float v1chi2;
  int   v1status;
  int   v1ndof;
  std::vector<float> v1pos;
  std::vector<float> v1cov;
  std::vector<int>   v1tkey;
  std::vector<int>   v1tidx;
  std::vector<float> v1tpd;
  std::vector<float> v1tip;
  std::vector<float> v1tipe;
  std::vector<float> v1tsip;
  std::vector<float> v1tchi2;
  std::vector<int>   v1tstatus;
  //
  // Split vertex 2
  //
  float v2chi2;
  int   v2status;
  int   v2ndof;
  std::vector<float> v2pos;
  std::vector<float> v2cov;
  std::vector<int>   v2tkey;
  std::vector<int>   v2tidx;
  std::vector<float> v2tpd;
  std::vector<float> v2tip;
  std::vector<float> v2tipe;
  std::vector<float> v2tsip;
  std::vector<float> v2tchi2;
  std::vector<int>   v2tstatus;
  //
  // MC 
  //
  int mcntks;
  std::vector<float> mcpos;
  std::vector<float> mcpos_unc;
  std::vector<float> mc_vtx_x;
  std::vector<float> mc_vtx_y;
  std::vector<float> mc_vtx_z;
  std::vector<float> mc_vtx_ux;
  std::vector<float> mc_vtx_uy;
  std::vector<float> mc_vtx_uz;
  std::vector<int>   mcpdg;
  std::vector<float> mcmom;
};

void SplitVertexNtuplizer::beginJob()
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
  // Original vertex
  //
  tree->Branch("vchi2"  , &vchi2,   "vchi2/F");
  tree->Branch("vstatus", &vstatus, "vstatus/I");
  tree->Branch("vndof"  , &vndof,   "vndof/I");
  tree->Branch("vpos"     , &vpos     );
  tree->Branch("vcov"     , &vcov     );
  tree->Branch("vtkey"    , &vtkey    );
  tree->Branch("vtpd"     , &vtpd     );
  tree->Branch("vtip"     , &vtip     );
  tree->Branch("vtipe"    , &vtipe    );
  tree->Branch("vtsip"    , &vtsip    );
  tree->Branch("vtchi2"   , &vtchi2   );
  tree->Branch("vtstatus" , &vtstatus );
  //
  tree->Branch("tk_vtx_x"   , &tk_vtx_x   );
  tree->Branch("tk_vtx_y"   , &tk_vtx_y   );
  tree->Branch("tk_vtx_z"   , &tk_vtx_z   );
  tree->Branch("tk_vtx_ux"  , &tk_vtx_ux  );
  tree->Branch("tk_vtx_uy"  , &tk_vtx_uy  );
  tree->Branch("tk_vtx_uz"  , &tk_vtx_uz  );
  tree->Branch("tk_vtx_xe"  , &tk_vtx_xe  );
  tree->Branch("tk_vtx_ye"  , &tk_vtx_ye  );
  tree->Branch("tk_vtx_ze"  , &tk_vtx_ze  );
  tree->Branch("tk_vtx_uxe" , &tk_vtx_uxe );
  tree->Branch("tk_vtx_uye" , &tk_vtx_uye );
  tree->Branch("tk_vtx_uze" , &tk_vtx_uze );
  tree->Branch("tk_end_x"   , &tk_end_x   );
  tree->Branch("tk_end_y"   , &tk_end_y   );
  tree->Branch("tk_end_z"   , &tk_end_z   );
  tree->Branch("tk_end_ux"  , &tk_end_ux  );
  tree->Branch("tk_end_uy"  , &tk_end_uy  );
  tree->Branch("tk_end_uz"  , &tk_end_uz  );
  tree->Branch("tk_end_xe"  , &tk_end_xe  );
  tree->Branch("tk_end_ye"  , &tk_end_ye  );
  tree->Branch("tk_end_ze"  , &tk_end_ze  );
  tree->Branch("tk_end_uxe" , &tk_end_uxe );
  tree->Branch("tk_end_uye" , &tk_end_uye );
  tree->Branch("tk_end_uze" , &tk_end_uze );
  tree->Branch("tnpoints"   , &tnpoints   );
  tree->Branch("tnvalid"    , &tnvalid    );
  tree->Branch("tlength"    , &tlength    );
  tree->Branch("tchi2"	    , &tchi2      );
  tree->Branch("tndof"	    , &tndof      );
  tree->Branch("tpdg"	    , &tpdg       );
  tree->Branch("tcontain"   , &tcontain   );
  tree->Branch("trvtx"	    , &trvtx      );
  tree->Branch("tmom"       , &tmom       );
  //
  // Split vertex 1
  //
  tree->Branch("v1chi2"  , &v1chi2,   "v1chi2/F");
  tree->Branch("v1status", &v1status, "v1status/I");
  tree->Branch("v1ndof"  , &v1ndof,   "v1ndof/I");
  tree->Branch("v1pos"     , &v1pos     );
  tree->Branch("v1cov"     , &v1cov     );
  tree->Branch("v1tkey"    , &v1tkey    );
  tree->Branch("v1tidx"    , &v1tidx    );
  tree->Branch("v1tpd"     , &v1tpd     );
  tree->Branch("v1tip"     , &v1tip     );
  tree->Branch("v1tipe"    , &v1tipe    );
  tree->Branch("v1tsip"    , &v1tsip    );
  tree->Branch("v1tchi2"   , &v1tchi2   );
  tree->Branch("v1tstatus" , &v1tstatus );
  //
  // Split vertex 2
  //
  tree->Branch("v2chi2"  , &v2chi2,   "v2chi2/F");
  tree->Branch("v2status", &v2status, "v2status/I");
  tree->Branch("v2ndof"  , &v2ndof,   "v2ndof/I");
  tree->Branch("v2pos"     , &v2pos     );
  tree->Branch("v2cov"     , &v2cov     );
  tree->Branch("v2tkey"    , &v2tkey    );
  tree->Branch("v2tidx"    , &v2tidx    );
  tree->Branch("v2tpd"     , &v2tpd     );
  tree->Branch("v2tip"     , &v2tip     );
  tree->Branch("v2tipe"    , &v2tipe    );
  tree->Branch("v2tsip"    , &v2tsip    );
  tree->Branch("v2tchi2"   , &v2tchi2   );
  tree->Branch("v2tstatus" , &v2tstatus );
  //
  // MC
  //
  tree->Branch("mcntks", &mcntks,"mcntks/I");
  tree->Branch("mcpos" , &mcpos );
  tree->Branch("mcpos_unc" , &mcpos_unc );
  tree->Branch("mc_vtx_x" , &mc_vtx_x );
  tree->Branch("mc_vtx_y" , &mc_vtx_y );
  tree->Branch("mc_vtx_z" , &mc_vtx_z );
  tree->Branch("mc_vtx_ux" , &mc_vtx_ux );
  tree->Branch("mc_vtx_uy" , &mc_vtx_uy );
  tree->Branch("mc_vtx_uz" , &mc_vtx_uz );
  tree->Branch("mcpdg" , &mcpdg );
  tree->Branch("mcmom" , &mcmom );
}

SplitVertexNtuplizer::SplitVertexNtuplizer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  inputPFLabel(p.get<art::InputTag>("inputPFLabel","pandoraNu")),
  inputTracksLabel(p.get<art::InputTag>("inputTracksLabel")),
  inputVertexLabel(p.get<art::InputTag>("inputVertexLabel")),
  inputVertexLabel1(p.get<art::InputTag>("inputVertexLabel1st")),
  inputVertexLabel2(p.get<art::InputTag>("inputVertexLabel2nd"))
{}

void SplitVertexNtuplizer::resetTree() {
  //
  // Event
  //
  run = -999;
  subrun = -999;
  eventid = -999;
  passSelII = -999;
  //
  // Original vertex
  //  
  vchi2   = -999;
  vstatus = -999;
  vndof   = -999;
  vpos.clear();
  vcov.clear();
  vtkey.clear();
  vtpd.clear();
  vtip.clear();
  vtipe.clear();
  vtsip.clear();
  vtchi2.clear();
  vtstatus.clear();
  tk_vtx_x.clear();
  tk_vtx_y.clear();
  tk_vtx_z.clear();
  tk_vtx_ux.clear();
  tk_vtx_uy.clear();
  tk_vtx_uz.clear();
  tk_vtx_xe.clear();
  tk_vtx_ye.clear();
  tk_vtx_ze.clear();
  tk_vtx_uxe.clear();
  tk_vtx_uye.clear();
  tk_vtx_uze.clear();
  tk_end_x.clear();
  tk_end_y.clear();
  tk_end_z.clear();
  tk_end_ux.clear();
  tk_end_uy.clear();
  tk_end_uz.clear();
  tk_end_xe.clear();
  tk_end_ye.clear();
  tk_end_ze.clear();
  tk_end_uxe.clear();
  tk_end_uye.clear();
  tk_end_uze.clear();
  tnpoints.clear();
  tnvalid.clear();
  tlength.clear();
  tchi2.clear();
  tndof.clear();
  tpdg.clear();
  tcontain.clear();
  trvtx.clear();
  tmom.clear();
  //
  // Split vertex 1
  //  
  v1chi2   = -999;
  v1status = -999;
  v1ndof   = -999;
  v1pos.clear();
  v1cov.clear();
  v1tkey.clear();
  v1tidx.clear();
  v1tpd.clear();
  v1tip.clear();
  v1tipe.clear();
  v1tsip.clear();
  v1tchi2.clear();
  v1tstatus.clear();
  //
  // Split vertex 2
  //  
  v2chi2   = -999;
  v2status = -999;
  v2ndof   = -999;
  v2pos.clear();
  v2cov.clear();
  v2tkey.clear();
  v2tidx.clear();
  v2tpd.clear();
  v2tip.clear();
  v2tipe.clear();
  v2tsip.clear();
  v2tchi2.clear();
  v2tstatus.clear();
  //
  // MC
  //
  mcntks = -999;
  mcpos.clear();
  mcpos_unc.clear();
  mc_vtx_x.clear();
  mc_vtx_y.clear();
  mc_vtx_z.clear();
  mc_vtx_ux.clear();
  mc_vtx_uy.clear();
  mc_vtx_uz.clear();
  mcpdg.clear();
  mcmom.clear();
}

void SplitVertexNtuplizer::analyze(art::Event const & e)
{
  //
  using namespace std;
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
  //
  art::ValidHandle<std::vector<recob::PFParticle> > inputPFParticle = e.getValidHandle<std::vector<recob::PFParticle> >(inputPFLabel);
  //
  auto assocVertex    = std::unique_ptr<art::FindManyP<recob::Vertex> >(new art::FindManyP<recob::Vertex>(inputPFParticle, e, inputVertexLabel));
  auto assocVertex1 = std::unique_ptr<art::FindManyP<recob::Vertex> >(new art::FindManyP<recob::Vertex>(inputPFParticle, e, inputVertexLabel1));
  auto assocVertex2 = std::unique_ptr<art::FindManyP<recob::Vertex> >(new art::FindManyP<recob::Vertex>(inputPFParticle, e, inputVertexLabel2));
  //
  art::ValidHandle<std::vector<recob::Vertex> > inputVertex    = e.getValidHandle<std::vector<recob::Vertex> >(inputVertexLabel);
  art::ValidHandle<std::vector<recob::Vertex> > inputVertex1 = e.getValidHandle<std::vector<recob::Vertex> >(inputVertexLabel1);
  art::ValidHandle<std::vector<recob::Vertex> > inputVertex2 = e.getValidHandle<std::vector<recob::Vertex> >(inputVertexLabel2);
  auto assocTracks    = std::unique_ptr<art::FindManyP<recob::Track, recob::VertexAssnMeta> >(new art::FindManyP<recob::Track, recob::VertexAssnMeta>(inputVertex, e, inputVertexLabel));
  auto assocTracks1 = std::unique_ptr<art::FindManyP<recob::Track, recob::VertexAssnMeta> >(new art::FindManyP<recob::Track, recob::VertexAssnMeta>(inputVertex1, e, inputVertexLabel1));
  auto assocTracks2 = std::unique_ptr<art::FindManyP<recob::Track, recob::VertexAssnMeta> >(new art::FindManyP<recob::Track, recob::VertexAssnMeta>(inputVertex2, e, inputVertexLabel2));
  //
  art::InputTag McTruthInputTag("generator");
  const std::vector<simb::MCTruth>* mcTruth = 0;
  if (e.isRealData()==0) {
    mcTruth = e.getValidHandle<std::vector<simb::MCTruth> >(McTruthInputTag).product();
  }
  //
  for (unsigned int iPF = 0; iPF < inputPFParticle->size(); ++iPF) {
    //
    const std::vector<art::Ptr<recob::Vertex> >& VertexVec  = assocVertex->at(iPF);
    const std::vector<art::Ptr<recob::Vertex> >& VertexVec1 = assocVertex1->at(iPF);
    const std::vector<art::Ptr<recob::Vertex> >& VertexVec2 = assocVertex2->at(iPF);
    //
    if (VertexVec.size()!=1 || VertexVec1.size()!=1 || VertexVec2.size()!=1) continue;
    //
    resetTree();
    //
    run = e.run();
    subrun = e.subRun();
    eventid = e.event();
    passSelII = filter->at(ftp).accept();

    auto vtx  = VertexVec[0];
    if (0) std::cout << "original pos=" << vtx->position() << std::endl;
    auto tracks = assocTracks->at(vtx.key());
    auto metas  = assocTracks->data(vtx.key());
    vchi2   = vtx->chi2();
    vstatus = vtx->status();
    vndof   = vtx->ndof();
    vpos = { float(vtx->position().X()), float(vtx->position().Y()), float(vtx->position().Z()) };
    vcov.assign( vtx->covariance().begin(), vtx->covariance().end() );
    for (size_t it=0; it<tracks.size(); ++it) {
      auto track = tracks[it];
      auto meta = metas[it];
      if (0) std::cout << "track with npoints=" << track->NPoints() << " pos=" << track->Start() << " dir=" << track->StartDirection() << " sip=" << meta->impactParamSig() << " dist=" << meta->propDist() << " status=" << meta->status() << std::endl;
      vtkey.push_back(track.key());
      vtpd.push_back(meta->propDist());
      vtip.push_back(meta->impactParam());
      vtipe.push_back(meta->impactParamErr());
      vtsip.push_back(meta->impactParamSig());
      vtchi2.push_back(meta->chi2());
      vtstatus.push_back(meta->status());
      //
      tk_vtx_x.push_back(track->Start().X());
      tk_vtx_y.push_back(track->Start().Y());
      tk_vtx_z.push_back(track->Start().Z());
      tk_vtx_ux.push_back(track->StartDirection().X());
      tk_vtx_uy.push_back(track->StartDirection().Y());
      tk_vtx_uz.push_back(track->StartDirection().Z());
      auto c6vtx = track->VertexCovarianceGlobal6D();
      tk_vtx_xe.push_back(sqrt(c6vtx(0,0)));
      tk_vtx_ye.push_back(sqrt(c6vtx(1,1)));
      tk_vtx_ze.push_back(sqrt(c6vtx(2,2)));
      tk_vtx_uxe.push_back(sqrt(c6vtx(3,3)));
      tk_vtx_uye.push_back(sqrt(c6vtx(4,4)));
      tk_vtx_uze.push_back(sqrt(c6vtx(5,5)));
      tk_end_x.push_back(track->End().X());
      tk_end_y.push_back(track->End().Y());
      tk_end_z.push_back(track->End().Z());
      tk_end_ux.push_back(track->EndDirection().X());
      tk_end_uy.push_back(track->EndDirection().Y());
      tk_end_uz.push_back(track->EndDirection().Z());
      auto c6end = track->EndCovarianceGlobal6D();
      tk_end_xe.push_back(sqrt(c6end(0,0)));
      tk_end_ye.push_back(sqrt(c6end(1,1)));
      tk_end_ze.push_back(sqrt(c6end(2,2)));
      tk_end_uxe.push_back(sqrt(c6end(3,3)));
      tk_end_uye.push_back(sqrt(c6end(4,4)));
      tk_end_uze.push_back(sqrt(c6end(5,5)));
      tnpoints.push_back(track->NPoints());
      tnvalid.push_back(track->CountValidPoints());
      tlength.push_back(track->Length());
      tchi2.push_back(track->Chi2());
      tndof.push_back(track->Ndof());
      tpdg.push_back(track->ParticleId());
      auto starttk = track->Start();
      auto endtk = track->End();
      bool tk_contain = (starttk.X()>30.  && starttk.X()<230.  && endtk.X()>30.  && endtk.X()<230. &&
			 starttk.Y()>-85. && starttk.Y()<85.   && endtk.Y()>-85. && endtk.Y()<85.  &&
			 starttk.Z()>30.  && starttk.Z()<1010. && endtk.Z()>30.  && endtk.Z()<1010.);
      tcontain.push_back(tk_contain);
      trvtx.push_back( (track->Start()-vtx->position()).R() );
      tmom.push_back(track->StartMomentum());
    }
    //
    auto vtx1 = VertexVec1[0];
    if (0) std::cout << "split1 pos=" << vtx1->position() << std::endl;
    auto tracks1 = assocTracks1->at(vtx1.key());
    auto metas1  = assocTracks1->data(vtx1.key());
    v1chi2   = vtx1->chi2();
    v1status = vtx1->status();
    v1ndof   = vtx1->ndof();
    v1pos = { float(vtx1->position().X()), float(vtx1->position().Y()), float(vtx1->position().Z()) };
    v1cov.assign( vtx1->covariance().begin(), vtx1->covariance().end() );
    for (size_t it1=0; it1<tracks1.size(); ++it1) {
      auto track1 = tracks1[it1];
      auto meta1 = metas1[it1];
      if (0) std::cout << "track1 with npoints=" << track1->NPoints() << " pos=" << track1->Start() << " sip=" << meta1->impactParamSig() << " dist=" << meta1->propDist()  << " status=" << meta1->status() << std::endl;
      v1tkey.push_back(track1.key());
      v1tidx.push_back( std::find(vtkey.begin(),vtkey.end(),track1.key())-vtkey.begin() );
      v1tpd.push_back(meta1->propDist());
      v1tip.push_back(meta1->impactParam());
      v1tipe.push_back(meta1->impactParamErr());
      v1tsip.push_back(meta1->impactParamSig());
      v1tchi2.push_back(meta1->chi2());
      v1tstatus.push_back(meta1->status());
    }
    //
    auto vtx2 = VertexVec2[0];
    if (0) std::cout << "split2 pos=" << vtx2->position() << std::endl;
    auto tracks2 = assocTracks2->at(vtx2.key());
    auto metas2  = assocTracks2->data(vtx2.key());
    v2chi2   = vtx2->chi2();
    v2status = vtx2->status();
    v2ndof   = vtx2->ndof();
    v2pos = { float(vtx2->position().X()), float(vtx2->position().Y()), float(vtx2->position().Z()) };
    v2cov.assign( vtx2->covariance().begin(), vtx2->covariance().end() );
    for (size_t it2=0; it2<tracks2.size(); ++it2) {
      auto track2 = tracks2[it2];
      auto meta2 = metas2[it2];
      if (0) std::cout << "track2 with npoints=" << track2->NPoints() << " pos=" << track2->Start() << " sip=" << meta2->impactParamSig() << " dist=" << meta2->propDist()  << " status=" << meta2->status() << std::endl;
      v2tkey.push_back(track2.key());
      v2tidx.push_back( std::find(vtkey.begin(),vtkey.end(),track2.key())-vtkey.begin() );
      v2tpd.push_back(meta2->propDist());
      v2tip.push_back(meta2->impactParam());
      v2tipe.push_back(meta2->impactParamErr());
      v2tsip.push_back(meta2->impactParamSig());
      v2tchi2.push_back(meta2->chi2());
      v2tstatus.push_back(meta2->status());
    }
    //
    // MC
    //
    if (e.isRealData()==0) {
      Point_t nuvtxbare(mcTruth->at(0).GetNeutrino().Nu().Position().X(),mcTruth->at(0).GetNeutrino().Nu().Position().Y(),mcTruth->at(0).GetNeutrino().Nu().Position().Z());
      auto scecorr = SCE->GetPosOffsets(nuvtxbare);
      double g4Ticks = detClocks->TPCG4Time2Tick(mcTruth->at(0).GetNeutrino().Nu().T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
      double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
      double yOffset = scecorr.Y();
      double zOffset = scecorr.Z();
      Point_t nuvtx(nuvtxbare.X()+xOffset,nuvtxbare.Y()+yOffset,nuvtxbare.Z()+zOffset);
      if (0) std::cout << "nu vtx=" << nuvtxbare << " corr=" << nuvtx << std::endl;
      unsigned int ngp = 0;
      std::vector<int> pdgs;
      std::vector<float> moms;
      std::vector<float> vtx_x;
      std::vector<float> vtx_y;
      std::vector<float> vtx_z;
      std::vector<float> vtx_ux;
      std::vector<float> vtx_uy;
      std::vector<float> vtx_uz;
      for (int ip=0; ip<mcTruth->at(0).NParticles(); ++ip) {
	decltype(auto) gp = mcTruth->at(0).GetParticle(ip);
	if (gp.StatusCode()!=1) continue;
	if (std::abs(gp.PdgCode())!=13 && std::abs(gp.PdgCode())!=221 && std::abs(gp.PdgCode())!=321 && std::abs(gp.PdgCode())!=2212) continue;
	ngp++;
	pdgs.push_back(gp.PdgCode());
	moms.push_back(gp.P());
	vtx_x.push_back(gp.Vx()+xOffset);
	vtx_y.push_back(gp.Vy()+yOffset);
	vtx_z.push_back(gp.Vz()+zOffset);
	vtx_ux.push_back(gp.Px()/gp.P());
	vtx_uy.push_back(gp.Py()/gp.P());
	vtx_uz.push_back(gp.Pz()/gp.P());
	if (0) std::cout << "part=" << ip << " pdg=" << gp.PdgCode() << " status=" << gp.StatusCode() << " mom=" << gp.Mother() << " process=" << gp.Process() 
			 << " E=" << gp.E() << " P=" << gp.P() << " pos=(" << gp.Vx() << "," << gp.Vy() << "," << gp.Vz() << ")"
			 << " dir=(" << gp.Px()/gp.P() << "," << gp.Py()/gp.P() << "," << gp.Pz()/gp.P() << ")"
			 << std::endl;
      }
      if (0) std::cout << "number of daughters=" << ngp << std::endl;
      mcntks = ngp;
      mcpos = {float(nuvtx.X()), float(nuvtx.Y()), float(nuvtx.Z())};
      mcpos_unc = { float(nuvtxbare.X()), float(nuvtxbare.Y()), float(nuvtxbare.Z()) };
      mc_vtx_x.assign( vtx_x.begin(), vtx_x.end() );
      mc_vtx_y.assign( vtx_y.begin(), vtx_y.end() );
      mc_vtx_z.assign( vtx_z.begin(), vtx_z.end() );
      mc_vtx_ux.assign( vtx_ux.begin(), vtx_ux.end() );
      mc_vtx_uy.assign( vtx_uy.begin(), vtx_uy.end() );
      mc_vtx_uz.assign( vtx_uz.begin(), vtx_uz.end() );
      mcpdg.assign( pdgs.begin(), pdgs.end() );
      mcmom.assign( moms.begin(), moms.end() );
    }
    //
    tree->Fill();
    //
  }
  //
}

DEFINE_ART_MODULE(SplitVertexNtuplizer)
