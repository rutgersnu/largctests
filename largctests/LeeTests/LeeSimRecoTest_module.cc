////////////////////////////////////////////////////////////////////////
// Class:       LeeSimRecoTest
// Plugin Type: analyzer (art v2_05_01)
// File:        LeeSimRecoTest_module.cc
//
// Generated at Tue Jul 17 13:22:48 2018 by Giuseppe Cerati using cetskelgen
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

#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"

/*

TODO:
- check that the hit navigation is correct (make sure the hits associated from the gaushit collection are the same as those from the PFP) - DONE
- split PFParticles into the different hierarchies, down to secondaries - DONE
- make an ntuple
- run over the full dev sample

 */

class LeeSimRecoTest;

class LeeSimRecoTest : public art::EDAnalyzer {
public:
  explicit LeeSimRecoTest(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LeeSimRecoTest(LeeSimRecoTest const &) = delete;
  LeeSimRecoTest(LeeSimRecoTest &&) = delete;
  LeeSimRecoTest & operator = (LeeSimRecoTest const &) = delete;
  LeeSimRecoTest & operator = (LeeSimRecoTest &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  void beginJob() override;

  void resetTree();

private:

  TTree* tree;
  //
  int run, subrun, eventid;
  //
  float nue, nup, nux, nuy, nuz, scnux, scnuy, scnuz;
  int infidvol;
  //
  float ele, elp;
  std::vector<double> proe;
  std::vector<double> prop;
  std::vector<double> pioe;
  std::vector<double> piop;
  std::vector<double> neue;
  std::vector<double> neup;
  //
  int tothits;
  int nelhits;
  int ntotprhits;
  int ntotpihits;
  std::vector<int> nprhits;
  std::vector<int> npihits;
  //
  int tothitscr;
  int nelhitscr;
  int ntotprhitscr;
  int ntotpihitscr;
  std::vector<int> nprhitscr;
  std::vector<int> npihitscr;
  //
  int pftothits;
  int npfnus;
  std::vector<int> nupfpdg;
  std::vector<int> nupftothits;
  std::vector<float> nuvtxx;
  std::vector<float> nuvtxy;
  std::vector<float> nuvtxz;
  std::vector< std::vector<int> > pfpdg;
  std::vector< std::vector<int> > pfhits;
  std::vector< std::vector<int> > pfelhits;
  std::vector< std::vector<int> > pftotprhits;
  std::vector< std::vector<int> > pftotpihits;
  std::vector< std::vector<int> > pfprhits;
  std::vector< std::vector<int> > pfpihits;
  std::vector< std::vector<int> > pfprid;
  std::vector< std::vector<int> > pfpiid;
  std::vector< std::vector<int> > pfnpr;
  std::vector< std::vector<int> > pfnpi;
  std::vector< std::vector<float> > pfvtxx;
  std::vector< std::vector<float> > pfvtxy;
  std::vector< std::vector<float> > pfvtxz;
  std::vector<float> allvtxx;
  std::vector<float> allvtxy;
  std::vector<float> allvtxz;
  //
};


LeeSimRecoTest::LeeSimRecoTest(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void LeeSimRecoTest::resetTree() {
  run = -999;
  subrun = -999;
  eventid = -999;
  //
  nue = -999;
  nup = -999;
  nux = -999;
  nuy = -999;
  nuz = -999;
  scnux = -999;
  scnuy = -999;
  scnuz = -999;
  infidvol = -999;
  //
  ele = -999;
  elp = -999;
  proe.clear();
  prop.clear();
  pioe.clear();
  piop.clear();
  neue.clear();
  neup.clear();
  //
  tothits = -999;
  nelhits = -999;
  ntotprhits = -999;
  ntotpihits = -999;
  nprhits.clear();
  npihits.clear();
  //
  tothitscr = -999;
  nelhitscr = -999;
  ntotprhitscr = -999;
  ntotpihitscr = -999;
  nprhitscr.clear();
  npihitscr.clear();
  //
  pftothits = -999;
  npfnus = -999;
  nupfpdg.clear();
  nupftothits.clear();
  nuvtxx.clear();
  nuvtxy.clear();
  nuvtxz.clear();
  pfpdg.clear();
  pfhits.clear();
  pfelhits.clear();
  pftotprhits.clear();
  pftotpihits.clear();
  pfprhits.clear();
  pfpihits.clear();
  pfprid.clear();
  pfpiid.clear();
  pfnpr.clear();
  pfnpi.clear();
  pfvtxx.clear();
  pfvtxy.clear();
  pfvtxz.clear();
  allvtxx.clear();
  allvtxy.clear();
  allvtxz.clear();
}

void LeeSimRecoTest::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  //
  tree = tfs->make<TTree>("tree", "tree");
  //
  tree->Branch("run", &run,"run/I");
  tree->Branch("subrun", &subrun, "subrun/I");
  tree->Branch("eventid", &eventid, "eventid/I");
  //
  tree->Branch("nue", &nue, "nue/F");
  tree->Branch("nup", &nup, "nup/F");
  tree->Branch("nux", &nux, "nux/F");
  tree->Branch("nuy", &nuy, "nuy/F");
  tree->Branch("nuz", &nuz, "nuz/F");
  tree->Branch("scnux", &scnux, "scnux/F");
  tree->Branch("scnuy", &scnuy, "scnuy/F");
  tree->Branch("scnuz", &scnuz, "scnuz/F");
  tree->Branch("infidvol", &infidvol, "infidvol/I");
  //
  tree->Branch("ele", &ele, "ele/F");
  tree->Branch("elp", &elp, "elp/F");
  tree->Branch("proe", &proe);
  tree->Branch("prop", &prop);
  tree->Branch("pioe", &pioe);
  tree->Branch("piop", &piop);
  tree->Branch("neue", &neue);
  tree->Branch("neup", &neup);
  //
  tree->Branch("tothits", &tothits, "tothits/I");
  tree->Branch("nelhits", &nelhits, "nelhits/I");
  tree->Branch("ntotprhits", &ntotprhits, "ntotprhits/I");
  tree->Branch("ntotpihits", &ntotpihits, "ntotpihits/I");
  tree->Branch("nprhits", &nprhits);
  tree->Branch("npihits", &npihits);
  //
  tree->Branch("tothitscr", &tothitscr, "tothitscr/I");
  tree->Branch("nelhitscr", &nelhitscr, "nelhitscr/I");
  tree->Branch("ntotprhitscr", &ntotprhitscr, "ntotprhitscr/I");
  tree->Branch("ntotpihitscr", &ntotpihitscr, "ntotpihitscr/I");
  tree->Branch("nprhitscr", &nprhitscr);
  tree->Branch("npihitscr", &npihitscr);
  //
  tree->Branch("pftothits", &pftothits, "pftothits/I");
  tree->Branch("npfnus", &npfnus, "npfnus/I");
  tree->Branch("nupfpdg", &nupfpdg);
  tree->Branch("nuvtxx", &nuvtxx);
  tree->Branch("nuvtxy", &nuvtxy);
  tree->Branch("nuvtxz", &nuvtxz);
  tree->Branch("nupftothits", &nupftothits);
  tree->Branch("pfpdg", &pfpdg);
  tree->Branch("pfhits", &pfhits);
  tree->Branch("pfelhits", &pfelhits);
  tree->Branch("pftotprhits", &pftotprhits);
  tree->Branch("pftotpihits", &pftotpihits);
  tree->Branch("pfprhits", &pfprhits);
  tree->Branch("pfpihits", &pfpihits);
  tree->Branch("pfprid", &pfprid);
  tree->Branch("pfpiid", &pfpiid);
  tree->Branch("pfnpr", &pfnpr);
  tree->Branch("pfnpi", &pfnpi);
  tree->Branch("pfvtxx", &pfvtxx);
  tree->Branch("pfvtxy", &pfvtxy);
  tree->Branch("pfvtxz", &pfvtxz);
  //
  tree->Branch("allvtxx", &allvtxx);
  tree->Branch("allvtxy", &allvtxy);
  tree->Branch("allvtxz", &allvtxz);
  //
  // tree->Branch("", &, "/I");
  // tree->Branch("", &, "/F");
  // tree->Branch("", &);
}

void LeeSimRecoTest::analyze(art::Event const & e)
{
  // Implementation of required member function here.

  resetTree();
  //
  run = e.run();
  subrun = e.subRun();
  eventid = e.event();

  int pdglep = -1;
  double elep = -1.;
  double plep = -1.;
  bool isinfidvol = false;
  std::vector<double> epro;
  std::vector<double> epio;
  std::vector<double> eneu;
  std::vector<double> ppro;
  std::vector<double> ppio;
  std::vector<double> pneu;

  detinfo::DetectorProperties const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
  detinfo::DetectorClocks const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>(); 

  // generator, simb::MCTruth
  art::InputTag GenInputTag("generator");
  art::ValidHandle<std::vector<simb::MCTruth> > inputMCTruth = e.getValidHandle<std::vector<simb::MCTruth> >(GenInputTag);
  /// std::cout << "N mct=" << inputMCTruth->size() << std::endl;
  if (inputMCTruth->size()<1) return;
  const auto& mct = inputMCTruth->at(0);
  pdglep = mct.GetNeutrino().Lepton().PdgCode();
  elep = mct.GetNeutrino().Lepton().E();
  plep = mct.GetNeutrino().Lepton().P();
  /// std::cout << "Neutrino E=" << mct.GetNeutrino().Nu().E() << " vtx=" << mct.GetNeutrino().Lepton().Vx() << ", " << mct.GetNeutrino().Lepton().Vy() << ", " << mct.GetNeutrino().Lepton().Vz() << " eleid=" << mct.GetNeutrino().Lepton().TrackId() << std::endl;
  // infidvol = (mct.GetNeutrino().Lepton().Vx()>0. && mct.GetNeutrino().Lepton().Vx()<260 && mct.GetNeutrino().Lepton().Vy()>-115. && mct.GetNeutrino().Lepton().Vy()<115. && mct.GetNeutrino().Lepton().Vz()>0. && mct.GetNeutrino().Lepton().Vz()<1040.);
  isinfidvol = (mct.GetNeutrino().Lepton().Vx()>30. && mct.GetNeutrino().Lepton().Vx()<220 && mct.GetNeutrino().Lepton().Vy()>-85. && mct.GetNeutrino().Lepton().Vy()<85. && mct.GetNeutrino().Lepton().Vz()>30. && mct.GetNeutrino().Lepton().Vz()<1010.);
  for ( int im=0; im<mct.NParticles(); ++im ) {
    const auto& mcp = mct.GetParticle(im);
    if (mcp.StatusCode()==1) {
       // std::cout << "id=" << mcp.TrackId() << " status=" << mcp.StatusCode() << " pdg=" << mcp.PdgCode() << " mother= " << mcp.Mother() << " E=" << mcp.E() << " P=" << mcp.P() << std::endl;
      if (std::abs(mcp.PdgCode())==211) {
	epio.push_back(mcp.E());
	ppio.push_back(mcp.P());
      }
      if (std::abs(mcp.PdgCode())==2212) {
	epro.push_back(mcp.E());
	ppro.push_back(mcp.P());
      }
      if (std::abs(mcp.PdgCode())==2112) {
	eneu.push_back(mcp.E());
	pneu.push_back(mcp.P());
      }
    }
  }
  if (abs(pdglep)!=11 && abs(pdglep)!=13) return;

  if (!isinfidvol) return;

  //
  nue = mct.GetNeutrino().Nu().E();
  nup = mct.GetNeutrino().Nu().P();
  nux = mct.GetNeutrino().Lepton().Vx();
  nuy = mct.GetNeutrino().Lepton().Vy();
  nuz = mct.GetNeutrino().Lepton().Vz();
  auto scecorr = SCE->GetPosOffsets({nux,nuy,nuz});
  double g4Ticks = detClocks->TPCG4Time2Tick(mct.GetNeutrino().Lepton().T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
  double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
  double yOffset = scecorr.Y();
  double zOffset = scecorr.Z();
  recob::tracking::Point_t mcpos(nux+xOffset,nuy+yOffset,nuz+zOffset);
  scnux = mcpos.X();
  scnuy = mcpos.Y();
  scnuz = mcpos.Z();
  infidvol = isinfidvol;
  //
  ele = elep;
  elp = plep;
  proe = epro;
  prop = ppro;
  pioe = epio;
  piop = ppio;
  neue = eneu;
  neup = pneu;
  //

  //need to find the id of the primary MCParticles in the right collection to be matched with the backtracker map
  int eleid = -1;
  std::vector<int> swrid;//needed to capture particles from the electron shower
  std::vector<int> proid;
  std::vector<int> pioid;
  std::vector<int> neuid;
  // largeant, simb::MCParticle
  art::InputTag MCInputTag("largeant");
  art::ValidHandle<std::vector<simb::MCParticle> > inputMCParticle = e.getValidHandle<std::vector<simb::MCParticle> >(MCInputTag);
  /// std::cout << "N mcp=" << inputMCParticle->size() << std::endl;
  for ( unsigned int im=0; im<inputMCParticle->size(); ++im ) {
    const auto& mcp = inputMCParticle->at(im);
    if (mcp.StatusCode()==1) {
      // std::cout << "id=" << mcp.TrackId() << " status=" << mcp.StatusCode()<< " pdg=" << mcp.PdgCode() << " mother= " << mcp.Mother() << " E=" << mcp.E() << std::endl;
      if (mcp.Mother()==0 && mcp.StatusCode()==1) {
	// /// std::cout << "id=" << mcp.TrackId() << " status=" << mcp.StatusCode()<< " pdg=" << mcp.PdgCode() << " mother= " << mcp.Mother() << " E=" << mcp.E() << std::endl;
	if (abs(mcp.PdgCode())==11 && std::abs(elep-mcp.E())<0.00001) { 
	  eleid = mcp.TrackId();
	  swrid.push_back(eleid);
	} else if (abs(mcp.PdgCode())==2212) {
	  for (unsigned int ip=0; ip<ppro.size(); ip++) {
	    if (std::abs(ppro[ip]-mcp.P())<0.00000001) {
	      //  /// std::cout << "ip=" << ip << " im=" << im << " epro[ip]=" << epro[ip] << " mcp.E()=" << mcp.E() << " mcp.P()=" << mcp.P() << " mother= " << mcp.Mother() << std::endl;
	      proid.push_back(mcp.TrackId());

	    }
	  }
	} else if (abs(mcp.PdgCode())==211) {
	  for (unsigned int ip=0; ip<ppio.size(); ip++) {
	    if (std::abs(ppio[ip]-mcp.P())<0.00000001) pioid.push_back(mcp.TrackId());
	  }
	} else if (abs(mcp.PdgCode())==2112) {
	  for (unsigned int ip=0; ip<pneu.size(); ip++) {
	    if (std::abs(pneu[ip]-mcp.P())<0.00000001) neuid.push_back(mcp.TrackId());
	  }
	}
      }
      if (std::find(swrid.begin(), swrid.end(), mcp.Mother())!=swrid.end()) {
	swrid.push_back(mcp.TrackId());
	// std::cout << "adding swr particle with mother=" << mcp.Mother() << " and id=" << mcp.TrackId() << std::endl;
      }
    }
  }

  /// std::cout << "pdglep=" << pdglep << " elep=" << elep << " plep=" << plep << " eleid=" << eleid << " npr=" << ppro.size() << " npi=" << ppio.size() << " nne=" << pneu.size() << std::endl;

  //  /// std::cout << "ppro.size()=" << ppro.size() << " proid.size()=" << proid.size() << std::endl;

  // gaushit,  recob::Hits
  art::InputTag HitInputTag("gaushit");
  art::InputTag HitTruthInputTag("gaushitTruthMatch");
  art::ValidHandle<std::vector<recob::Hit> > inputHits = e.getValidHandle<std::vector<recob::Hit> >(HitInputTag);
  std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> > assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>(inputHits, e, HitTruthInputTag));
   /// std::cout << "N hits=" << inputHits->size() << std::endl;
  int nelehits = 0;
  std::vector<int> nprohits(ppro.size(),0);
  std::vector<int> npiohits(ppio.size(),0);
  for ( unsigned int ih=0; ih<inputHits->size(); ih++) {
    // const auto& hit = inputHits->at(ih);
    auto assmcp = assocMCPart->at(ih);
    // if (hit.Channel()==1038)  /// std::cout << "gh hit id=" << ih << " channel=" << hit.Channel() << " time[s,p,e]=" << hit.StartTick() << ", " << hit.PeakTime() << ", " << hit.EndTick() 
    // 				       << " rms=" << hit.RMS() << " mult=" << hit.Multiplicity() << " nAssocMCP=" << assmcp.size() << std::endl;
    for (unsigned int ia=0; ia<assmcp.size(); ++ia) {
      auto mcp = assmcp[ia];
      // if (hit.Channel()==1038 && mcp->StatusCode()==1)  /// std::cout << "mcp#" << ia << "id=" << mcp->TrackId() << " pdg=" << mcp->PdgCode() << " mother= " << mcp->Mother() << " E=" << mcp->E()  << std::endl;
      // std::cout << "mcp#" << ia << "id=" << mcp->TrackId() << " pdg=" << mcp->PdgCode() << " mother= " << mcp->Mother() << " E=" << mcp->E()  << std::endl;
      // if (eleid == mcp->TrackId()) {
      if (std::find(swrid.begin(),swrid.end(),mcp->TrackId())!=swrid.end()) {
	// const auto& hit = inputHits->at(ih);
	//  /// std::cout << "hit id=" << ih << " channel=" << hit.Channel() << " time=" << hit.PeakTime() << " nAssocMCP=" << assmcp.size() << std::endl;
	//  /// std::cout << "id=" << mcp->TrackId() << " pdg=" << mcp->PdgCode() << " mother= " << mcp->Mother() << " E=" << mcp->E()  << std::endl;
	nelehits++;
	break;
      } else if (std::find(proid.begin(),proid.end(),mcp->TrackId())!=proid.end()) {
	nprohits[std::find(proid.begin(),proid.end(),mcp->TrackId())-proid.begin()]++;
	break;
      } else if (std::find(pioid.begin(),pioid.end(),mcp->TrackId())!=pioid.end()) {
	npiohits[std::find(pioid.begin(),pioid.end(),mcp->TrackId())-pioid.begin()]++;
	break;
      }
    }
  }
  int countprohits = 0;
  for (auto n : nprohits) countprohits+=n;
  int countpiohits = 0;
  for (auto n : npiohits) countpiohits+=n;
  //std::cout << "N electron hits=" << nelehits << " proton=" << countprohits << " pion=" << countpiohits << std::endl;

  //
  tothits = inputHits->size();
  nelhits = nelehits;
  ntotprhits = countprohits;
  ntotpihits = countpiohits;
  nprhits = nprohits;
  npihits = npiohits;
  //

  // pandoraCosmicHitRemoval,  recob::Hits
  art::InputTag CRHitInputTag("pandoraCosmicHitRemoval");
  art::InputTag CRHitTruthInputTag("crHitRemovalTruthMatch");
  art::ValidHandle<std::vector<recob::Hit> > inputCRHits = e.getValidHandle<std::vector<recob::Hit> >(CRHitInputTag);
  //map the crhits to gaushits, the BackTrackerHitMatchingData is buggy!
  std::vector<int> crindex(inputCRHits->size());
  for ( unsigned int ic=0; ic<inputCRHits->size(); ic++) {
    const auto& crhit = inputCRHits->at(ic);
    for ( unsigned int ig=0; ig<inputHits->size(); ig++) {
      const auto& gahit = inputHits->at(ig);
      if (crhit.Channel()==gahit.Channel() && std::abs(crhit.PeakTime()-gahit.PeakTime())<0.000001) crindex[ic] = ig;
    }
  }
  // std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> > assocMCPartCR = std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>(inputCRHits, e, CRHitTruthInputTag));
   /// std::cout << "N CR hits=" << inputCRHits->size() << std::endl;
  int nelehitsCR = 0;
  std::vector<int> nprohitsCR(ppro.size(),0);
  std::vector<int> npiohitsCR(ppio.size(),0);
  for ( unsigned int ih=0; ih<inputCRHits->size(); ih++) {
    // const auto& hit = inputCRHits->at(ih);
    // auto assmcp = assocMCPartCR->at(ih);
    auto assmcp = assocMCPart->at(crindex[ih]);
    // if (hit.Channel()==1038)  /// std::cout << "cr hit id=" << ih << " channel=" << hit.Channel() << " time[s,p,e]=" << hit.StartTick() << ", " << hit.PeakTime() << ", " << hit.EndTick() 
    // 				       << " rms=" << hit.RMS() << " mult=" << hit.Multiplicity() << " nAssocMCP=" << assmcp.size() << std::endl;
    for (unsigned int ia=0; ia<assmcp.size(); ++ia) {
      auto mcp = assmcp[ia];
      // if (hit.Channel()==1038 && mcp->StatusCode()==1)  /// std::cout << "mcp#" << ia << "id=" << mcp->TrackId() << " pdg=" << mcp->PdgCode() << " mother= " << mcp->Mother() << " E=" << mcp->E()  << std::endl;
      // if (eleid == mcp->TrackId()) {
      if (std::find(swrid.begin(),swrid.end(),mcp->TrackId())!=swrid.end()) {
	//  /// std::cout << "cr hit id=" << ih << " channel=" << hit.Channel() << " time=" << hit.PeakTime() << " nAssocMCP=" << assmcp.size() << std::endl;
	//  /// std::cout << "id=" << mcp->TrackId() << " pdg=" << mcp->PdgCode() << " mother= " << mcp->Mother() << " E=" << mcp->E()  << std::endl;
	nelehitsCR++;
	break;
      } else if (std::find(proid.begin(),proid.end(),mcp->TrackId())!=proid.end()) {
	bool test = ( int(nprohitsCR.size())>int(std::find(proid.begin(),proid.end(),mcp->TrackId())-proid.begin()) );
	if (!test) {
	  std::cout << __FILE__ << " " << __LINE__ << std::endl;
	  exit(1);
	}
	nprohitsCR[std::find(proid.begin(),proid.end(),mcp->TrackId())-proid.begin()]++;
	break;
      } else if (std::find(pioid.begin(),pioid.end(),mcp->TrackId())!=pioid.end()) {
	bool test = ( int(npiohitsCR.size()) > int(std::find(pioid.begin(),pioid.end(),mcp->TrackId())-pioid.begin()) );
	if (!test) {
	  std::cout << __FILE__ << " " << __LINE__ << std::endl;
	  exit(1);
	}
	npiohitsCR[std::find(pioid.begin(),pioid.end(),mcp->TrackId())-pioid.begin()]++;
	break;
      }
    }
  }
  int countprohitsCR = 0;
  for (auto n : nprohitsCR) countprohitsCR+=n;
  int countpiohitsCR = 0;
  for (auto n : npiohitsCR) countpiohitsCR+=n;
   /// std::cout << "N hits after CR electron=" << nelehitsCR << " proton=" << countprohitsCR << " pion=" << countpiohitsCR << std::endl;

  //
  tothitscr = inputCRHits->size();
  nelhitscr = nelehitsCR;
  ntotprhitscr = countprohitsCR;
  ntotpihitscr = countpiohitsCR;
  nprhitscr = nprohitsCR;
  npihitscr = npiohitsCR;
  //

  // pandoraNu PFParticles - will have to redo the MC truth matching by hand
  //
  int mypftothits = 0;
  int mynpfnus = 0;
  std::vector<int> mynupfpdg;
  std::vector<int> mynupftothits;
  std::vector<float> mynuvtxx;
  std::vector<float> mynuvtxy;
  std::vector<float> mynuvtxz;
  std::vector< std::vector<int> > mypfpdg;
  std::vector< std::vector<int> > mypfhits;
  std::vector< std::vector<int> > mypfelhits;
  std::vector< std::vector<int> > mypftotprhits;
  std::vector< std::vector<int> > mypftotpihits;
  std::vector< std::vector<int> > mypfprhits;
  std::vector< std::vector<int> > mypfpihits;
  std::vector< std::vector<int> > mypfprid;
  std::vector< std::vector<int> > mypfpiid;
  std::vector< std::vector<int> > mypfnpr;
  std::vector< std::vector<int> > mypfnpi;
  std::vector< std::vector<float> > mypfvtxx;
  std::vector< std::vector<float> > mypfvtxy;
  std::vector< std::vector<float> > mypfvtxz;
  art::InputTag PfInputTag("pandoraNu");
  art::ValidHandle<std::vector<recob::PFParticle> > inputPfParticle = e.getValidHandle<std::vector<recob::PFParticle> >(PfInputTag);
  art::ValidHandle<std::vector<recob::Cluster> > inputCluster = e.getValidHandle<std::vector<recob::Cluster> >(PfInputTag);
  auto assocCluster = std::unique_ptr<art::FindManyP<recob::Cluster> >(new art::FindManyP<recob::Cluster>(inputPfParticle, e, PfInputTag));
  auto assocHit = std::unique_ptr<art::FindManyP<recob::Hit> >(new art::FindManyP<recob::Hit>(inputCluster, e, PfInputTag));
  auto assocVertex = std::unique_ptr<art::FindManyP<recob::Vertex> >(new art::FindManyP<recob::Vertex>(inputPfParticle, e, PfInputTag));
  /// std::cout << "N pfps=" << inputPfParticle->size() << std::endl;
  for (unsigned int inpf=0; inpf<inputPfParticle->size(); ++inpf) {
    art::Ptr<recob::PFParticle> npfp(inputPfParticle,inpf);
    if (std::abs(npfp->PdgCode())!=12 && std::abs(npfp->PdgCode())!=14) continue;
    mynpfnus++;
    mynupfpdg.push_back(npfp->PdgCode());
    mynupftothits.push_back(0);
    const std::vector<art::Ptr<recob::Vertex> >& vertexVec = assocVertex->at(npfp.key());
    double xyz[3];
    vertexVec.front()->XYZ(xyz);
    mynuvtxx.push_back(xyz[0]);
    mynuvtxy.push_back(xyz[1]);
    mynuvtxz.push_back(xyz[2]);
    mypfpdg.push_back( std::vector<int>() );
    mypfhits.push_back( std::vector<int>() );
    mypfelhits.push_back( std::vector<int>() );
    mypftotprhits.push_back( std::vector<int>() );
    mypftotpihits.push_back( std::vector<int>() );
    mypfprhits.push_back( std::vector<int>() );
    mypfpihits.push_back( std::vector<int>() );
    mypfprid.push_back( std::vector<int>() );
    mypfpiid.push_back( std::vector<int>() );
    mypfnpr.push_back( std::vector<int>() );
    mypfnpi.push_back( std::vector<int>() );
    mypfvtxx.push_back( std::vector<float>() );
    mypfvtxy.push_back( std::vector<float>() );
    mypfvtxz.push_back( std::vector<float>() );
    /// std::cout << "neutrino pfp pdg=" << npfp->PdgCode() << " and nDaughters=" << npfp->NumDaughters() << std::endl;
    std::vector<size_t> dpfps = npfp->Daughters();
    // include also pfps one more layer down in the hierarchy (secondaries)
    for (unsigned int ipf : dpfps) {
      art::Ptr<recob::PFParticle> ppfp(inputPfParticle,ipf);
      const std::vector<size_t>& dpfps_tmp = ppfp->Daughters();
      for (auto k : dpfps_tmp) dpfps.push_back(k);
    }
    // now let's go over primaries and secondaries
    for (unsigned int ipf : dpfps) {
      int npfhits = 0;
      int nelehitsPF = 0;
      std::vector<int> nprohitsPF(ppro.size(),0);
      std::vector<int> npiohitsPF(ppio.size(),0);
      art::Ptr<recob::PFParticle> ppfp(inputPfParticle,ipf);
      const std::vector<art::Ptr<recob::Cluster> >& clustVec = assocCluster->at(ppfp.key());
      for (unsigned int icl=0; icl<clustVec.size(); icl++) {
	art::Ptr<recob::Cluster> cl = clustVec[icl];
	const std::vector<art::Ptr<recob::Hit> >& hitVec = assocHit->at(cl.key());
	npfhits += hitVec.size();
	for (unsigned int ih=0; ih<hitVec.size(); ih++) {
	  art::Ptr<recob::Hit> hitp = hitVec[ih];
	  auto assmcp = assocMCPart->at(crindex[hitp.key()]);
	  for (unsigned int ia=0; ia<assmcp.size(); ++ia) {
	    auto mcp = assmcp[ia];
	    // if (eleid == mcp->TrackId()) {
	    if (std::find(swrid.begin(),swrid.end(),mcp->TrackId())!=swrid.end()) {
	      //  /// std::cout << "hit id=" << hitp.key() << " channel=" << hitp->Channel() << " time=" << hitp->PeakTime() << " nAssocMCP=" << assmcp.size() << std::endl;
	      nelehitsPF++;
	      break;
	    } else if (std::find(proid.begin(),proid.end(),mcp->TrackId())!=proid.end()) {
	      //  /// std::cout << "nprohitsPF.size()=" << nprohitsPF.size() << " std::find(proid.begin(),proid.end(),mcp->TrackId())-proid.begin()=" << std::find(proid.begin(),proid.end(),mcp->TrackId())-proid.begin() << std::endl;
	      bool test = ( int(nprohitsPF.size()) > int(std::find(proid.begin(),proid.end(),mcp->TrackId())-proid.begin()) );
	      if (!test) {
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		exit(1);
	      }
	      nprohitsPF[std::find(proid.begin(),proid.end(),mcp->TrackId())-proid.begin()]++;
	      break;
	    } else if (std::find(pioid.begin(),pioid.end(),mcp->TrackId())!=pioid.end()) {
	      //  /// std::cout << "npiohitsPF.size()=" << npiohitsPF.size() << " std::find(pioid.begin(),pioid.end(),mcp->TrackId())-pioid.begin()=" << std::find(pioid.begin(),pioid.end(),mcp->TrackId())-pioid.begin() << std::endl;
	      bool test = ( int(npiohitsPF.size()) > int(std::find(pioid.begin(),pioid.end(),mcp->TrackId())-pioid.begin()) );
	      if (!test) {
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		exit(1);
	      }
	      npiohitsPF[std::find(pioid.begin(),pioid.end(),mcp->TrackId())-pioid.begin()]++;
	      break;
	    }
	  }
	}
      }
      auto pfp = inputPfParticle->at(ipf);
      int countprohitsPF = 0;
      int inpr = 0;
      int inprn = 0;
      for (auto n : nprohitsPF) {
	countprohitsPF+=n;
	if (n>0) {
	  mypfprhits.back().push_back(n);
	  mypfprid.back().push_back(inpr);
	  inprn++;
	}
	inpr++;
      }
      int countpiohitsPF = 0;
      int inpi = 0;
      int inpin = 0;
      for (auto n : npiohitsPF) {
	countpiohitsPF+=n;
	if (n>0) {
	  mypfpihits.back().push_back(n);
	  mypfpiid.back().push_back(inpi);
	  inpin++;
	}
	inpi++;
      }
      /// std::cout << "pf #" << ppfp.key() << " pdg=" << pfp.PdgCode() << " has nhits=" << npfhits << " ele=" << nelehitsPF << " pro=" << countprohitsPF << " pio=" << countpiohitsPF << std::endl;
      mypftothits+=npfhits;
      mynupftothits.back()+=npfhits;
      mypfpdg.back().push_back(pfp.PdgCode());
      mypfhits.back().push_back(npfhits);
      mypfelhits.back().push_back(nelehitsPF);
      mypftotprhits.back().push_back(countprohitsPF);
      mypftotpihits.back().push_back(countpiohitsPF);
      mypfnpr.back().push_back(inprn);
      mypfnpi.back().push_back(inpin);
      const std::vector<art::Ptr<recob::Vertex> >& vertexVec = assocVertex->at(ppfp.key());
      /// std::cout << "pfp nvtx=" << vertexVec.size() << std::endl;
      double xyz[3] = {-999,-999,-999};
      if (vertexVec.size()>0) {
	vertexVec.front()->XYZ(xyz);
	mypfvtxx.back().push_back(xyz[0]);
	mypfvtxy.back().push_back(xyz[1]);
	mypfvtxz.back().push_back(xyz[2]);
      }
      /// std::cout << "end pfp" << std::endl;
    }
  }
  //
  pftothits = mypftothits;
  npfnus = mynpfnus;
  nupfpdg = mynupfpdg;
  nuvtxx = mynuvtxx;
  nuvtxy = mynuvtxy;
  nuvtxz = mynuvtxz;
  nupftothits = mynupftothits;
  pfpdg = mypfpdg;
  pfhits = mypfhits;
  pfelhits = mypfelhits;
  pftotprhits = mypftotprhits;
  pftotpihits = mypftotpihits;
  pfprhits = mypfprhits;
  pfpihits = mypfpihits;
  pfprid = mypfprid;
  pfpiid = mypfpiid;
  pfnpr = mypfnpr;
  pfnpi = mypfnpi;
  pfvtxx = mypfvtxx;
  pfvtxy = mypfvtxy;
  pfvtxz = mypfvtxz;
  //
  art::ValidHandle<std::vector<recob::Vertex> > inputVertex = e.getValidHandle<std::vector<recob::Vertex> >(PfInputTag);
  for (auto vtx : *inputVertex) {
    double xyz[3] = {-999,-999,-999};
    vtx.XYZ(xyz);
    allvtxx.push_back(xyz[0]);
    allvtxy.push_back(xyz[1]);
    allvtxz.push_back(xyz[2]);
  }
  //
  tree->Fill();
  //
}

DEFINE_ART_MODULE(LeeSimRecoTest)
