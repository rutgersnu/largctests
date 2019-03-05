////////////////////////////////////////////////////////////////////////
// Class:       TestBackTracking
// Plugin Type: analyzer (art v3_01_02)
// File:        TestBackTracking_module.cc
//
// Generated at Tue Mar  5 11:09:28 2019 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_05_01.
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

#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"

class TestBackTracking;


class TestBackTracking : public art::EDAnalyzer {
public:
  explicit TestBackTracking(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TestBackTracking(TestBackTracking const&) = delete;
  TestBackTracking(TestBackTracking&&) = delete;
  TestBackTracking& operator=(TestBackTracking const&) = delete;
  TestBackTracking& operator=(TestBackTracking&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginJob() override;

  void resetTree();

private:

  // Declare member data here.
  TTree* tree;
  //
  int run, subrun, eventid;
  //
  int plane;
  int wire;
  int channel;
  float peakTime;
  float sigmaPeakTime;
  float rms;
  float peakAmplitude;
  float sigmaPeakAmplitude;
  float summedADC;
  float integral;
  float sigmaIntegral;
  int multiplicity;
  int localIndex;
  float goodnessOfFit;
  int degreesOfFreedom;
  int startTick;
  int endTick;
  //
  int nmcp;
  std::vector<int> mcp_id;
  std::vector<int> mcp_pdg;
  std::vector<int> mcp_mother;
  std::vector<float> mcp_p;
  std::vector<float> mcp_ux;
  std::vector<float> mcp_uy;
  std::vector<float> mcp_uz;
  std::vector<float> mcp_x;
  std::vector<float> mcp_y;
  std::vector<float> mcp_z;
  std::vector<float> mcp_ideFraction;
  std::vector<int> mcp_isMaxIDE;
  std::vector<float> mcp_ideNFraction;
  std::vector<int> mcp_isMaxIDEN;
  std::vector<float> mcp_numElectrons;
  std::vector<float> mcp_energy;

};


TestBackTracking::TestBackTracking(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void TestBackTracking::resetTree() {
  run = -999;
  subrun = -999;
  eventid = -999;
  //
  plane = -999;
  wire = -999;
  channel = -999;
  peakTime = -999;
  sigmaPeakTime = -999;
  rms = -999;
  peakAmplitude = -999;
  sigmaPeakAmplitude = -999;
  summedADC = -999;
  integral = -999;
  sigmaIntegral = -999;
  multiplicity = -999;
  localIndex = -999;
  goodnessOfFit = -999;
  degreesOfFreedom = -999;
  startTick = -999;
  endTick = -999;
  //
  nmcp = -999;
  mcp_id.clear();
  mcp_pdg.clear();
  mcp_mother.clear();
  mcp_p.clear();
  mcp_ux.clear();
  mcp_uy.clear();
  mcp_uz.clear();
  mcp_x.clear();
  mcp_y.clear();
  mcp_z.clear();
  mcp_ideFraction.clear();
  mcp_isMaxIDE.clear();
  mcp_ideNFraction.clear();
  mcp_isMaxIDEN.clear();
  mcp_numElectrons.clear();
  mcp_energy.clear();
  //
}

void TestBackTracking::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  //
  tree = tfs->make<TTree>("tree", "tree");
  //
  tree->Branch("run", &run,"run/I");
  tree->Branch("subrun", &subrun, "subrun/I");
  tree->Branch("eventid", &eventid, "eventid/I");
  //
  tree->Branch("plane", &plane, "plane/I");
  tree->Branch("wire", &wire, "wire/I");
  tree->Branch("channel", &channel, "channel/I");
  tree->Branch("peakTime", &peakTime, "peakTime/F");
  tree->Branch("sigmaPeakTime", &sigmaPeakTime, "sigmaPeakTime/F");
  tree->Branch("rms", &rms, "rms/F");
  tree->Branch("peakAmplitude", &peakAmplitude, "peakAmplitude/F");
  tree->Branch("sigmaPeakAmplitude",&sigmaPeakAmplitude,"sigmaPeakAmplitude/F");
  tree->Branch("summedADC", &summedADC, "summedADC/F");
  tree->Branch("integral", &integral, "integral/F");
  tree->Branch("sigmaIntegral", &sigmaIntegral, "sigmaIntegral/F");
  tree->Branch("multiplicity", &multiplicity, "multiplicity/I");
  tree->Branch("localIndex", &localIndex, "localIndex/I");
  tree->Branch("goodnessOfFit", &goodnessOfFit, "goodnessOfFit/F");
  tree->Branch("degreesOfFreedom  ",&degreesOfFreedom,"degreesOfFreedom/I");
  tree->Branch("startTick", &startTick, "startTick/I");
  tree->Branch("endTick",&endTick,"endTick/I");
  //
  tree->Branch("nmcp", &nmcp, "nmcp/I");
  tree->Branch("mcp_id", &mcp_id);
  tree->Branch("mcp_pdg", &mcp_pdg);
  tree->Branch("mcp_mother", &mcp_mother);
  tree->Branch("mcp_p", &mcp_p);
  tree->Branch("mcp_ux", &mcp_ux);
  tree->Branch("mcp_uy", &mcp_uy);
  tree->Branch("mcp_uz", &mcp_uz);
  tree->Branch("mcp_x", &mcp_x);
  tree->Branch("mcp_y", &mcp_y);
  tree->Branch("mcp_z", &mcp_z);
  tree->Branch("mcp_ideFraction", &mcp_ideFraction);
  tree->Branch("mcp_isMaxIDE", &mcp_isMaxIDE);
  tree->Branch("mcp_ideNFraction", &mcp_ideNFraction);
  tree->Branch("mcp_isMaxIDEN", &mcp_isMaxIDEN);
  tree->Branch("mcp_numElectrons", &mcp_numElectrons);
  tree->Branch("mcp_energy", &mcp_energy);
}

void TestBackTracking::analyze(art::Event const& e)
{
  // gaushit,  recob::Hits
  art::InputTag HitInputTag("gaushit");
  art::InputTag HitTruthInputTag("gaushitTruthMatch");
  auto hits = proxy::getCollection<std::vector<recob::Hit> > (e, HitInputTag, 
							      proxy::withAssociatedMeta<simb::MCParticle,anab::BackTrackerHitMatchingData>(HitTruthInputTag));
  //
  for (auto h : hits) {
    const auto& assocMCP = h.get<simb::MCParticle>();
    //std::cout << "hit time=" << h->PeakTime() << " nAssocMcp=" << assocMCP.size() << std::endl;
    resetTree();
    //
    run = e.run();
    subrun = e.subRun();
    eventid = e.event();
    //
    plane = h->WireID().Plane;
    wire = h->WireID().Wire;
    channel = h->Channel();
    peakTime = h->PeakTime();
    sigmaPeakTime = h->SigmaPeakTime();
    rms = h->RMS();
    peakAmplitude = h->PeakAmplitude();
    sigmaPeakAmplitude = h->SigmaPeakAmplitude();
    summedADC = h->SummedADC();
    integral = h->Integral();
    sigmaIntegral = h->SigmaIntegral();
    multiplicity = h->Multiplicity();
    localIndex = h->LocalIndex();
    goodnessOfFit = h->GoodnessOfFit();
    degreesOfFreedom = h->DegreesOfFreedom();
    startTick = h->StartTick();
    endTick = h->EndTick();
    //
    nmcp = assocMCP.size();
    for (auto mcp : assocMCP) {
      mcp_id.push_back(mcp->TrackId());
      mcp_pdg.push_back(mcp->PdgCode());
      mcp_mother.push_back(mcp->Mother());
      // std::cout << mcp->Process() << " " << mcp->Mother() << " " << mcp->StatusCode() << std::endl;
      mcp_p.push_back(mcp->P());
      mcp_ux.push_back(mcp->Px()/mcp->P());
      mcp_uy.push_back(mcp->Py()/mcp->P());
      mcp_uz.push_back(mcp->Pz()/mcp->P());
      mcp_x.push_back(mcp->Vx());
      mcp_y.push_back(mcp->Vy());
      mcp_z.push_back(mcp->Vz());
      auto md = mcp.data();
      mcp_ideFraction.push_back(md.ideFraction);
      mcp_isMaxIDE.push_back(md.isMaxIDE);
      mcp_ideNFraction.push_back(md.ideNFraction);
      mcp_isMaxIDEN.push_back(md.isMaxIDEN);
      mcp_numElectrons.push_back(md.numElectrons);
      mcp_energy.push_back(md.energy);
    }
    //
    tree->Fill();
  }
}

DEFINE_ART_MODULE(TestBackTracking)
