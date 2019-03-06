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
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"

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

  std::vector<size_t> pfpDaughtersKey(const art::Ptr<recob::PFParticle> ppfp,
				      art::ValidHandle<std::vector<recob::PFParticle> > inputPfParticle);

  void incrementCounts(const art::Ptr<recob::PFParticle> ppfp,
		       const std::unique_ptr<art::FindManyP<recob::Cluster> >& assocCluster,
		       const std::unique_ptr<art::FindManyP<recob::Hit> >& assocHit,
		       const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& assocMCPart,
		       const std::vector<int>& swrid, const std::vector<int>& proid, const std::vector<int>& pioid, const std::vector<int>& npiid, const std::vector<int>& neuid,
		       int& npfhits, int& nelehitsPF, std::vector<int>& nprohitsPF, std::vector<int>& npiohitsPF, int& nnpihitsPF, int& nneuhitsPF) const;

  void countSliceHits(const std::vector<art::Ptr<recob::Hit> >& hits,
		       const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& assocMCPart,
		       const std::vector<int>& swrid, const std::vector<int>& proid, const std::vector<int>& pioid, const std::vector<int>& npiid, const std::vector<int>& neuid,
		       int& nelehitsPF, int& nprohitsPF, int& npiohitsPF, int& nnpihitsPF, int& nneuhitsPF) const;

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
  std::vector<double> npie;
  std::vector<double> npip;
  std::vector<double> neue;
  std::vector<double> neup;
  //
  int tothits;
  int nelhits;
  int ntotprhits;
  int ntotpihits;
  int ntotnphits;
  int ntotnehits;
  std::vector<int> nprhits;
  std::vector<int> npihits;
  //
  int tothitscr;
  int nelhitscr;
  int ntotprhitscr;
  int ntotpihitscr;
  int ntotnphitscr;
  int ntotnehitscr;
  std::vector<int> nprhitscr;
  std::vector<int> npihitscr;
  //
  int pftothits;
  int npfnus;
  std::vector<int> nupfpdg;
  std::vector<int> nupftothits;
  std::vector<int> nusltothits;
  std::vector<int> nuslelhits;
  std::vector<int> nusltotprhits;
  std::vector<int> nusltotpihits;
  std::vector<int> nusltotnphits;
  std::vector<int> nusltotnehits;
  std::vector<float> nupfscore;
  std::vector<float> nuvtxx;
  std::vector<float> nuvtxy;
  std::vector<float> nuvtxz;
  std::vector< std::vector<int> > pfpdg;
  std::vector< std::vector<int> > pfhits;
  std::vector< std::vector<int> > pfelhits;
  std::vector< std::vector<int> > pftotprhits;
  std::vector< std::vector<int> > pftotpihits;
  std::vector< std::vector<int> > pftotnphits;
  std::vector< std::vector<int> > pftotnehits;
  std::vector< std::vector<int> > pfprhits;
  std::vector< std::vector<int> > pfpihits;
  std::vector< std::vector<int> > pfprid;
  std::vector< std::vector<int> > pfpiid;
  std::vector< std::vector<int> > pfnpr;
  std::vector< std::vector<int> > pfnpi;
  std::vector< std::vector<float> > pfscore;
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
  npie.clear();
  npip.clear();
  neue.clear();
  neup.clear();
  //
  tothits = -999;
  nelhits = -999;
  ntotprhits = -999;
  ntotpihits = -999;
  ntotnphits = -999;
  ntotnehits = -999;
  nprhits.clear();
  npihits.clear();
  //
  tothitscr = -999;
  nelhitscr = -999;
  ntotprhitscr = -999;
  ntotpihitscr = -999;
  ntotnphitscr = -999;
  ntotnehitscr = -999;
  nprhitscr.clear();
  npihitscr.clear();
  //
  pftothits = -999;
  npfnus = -999;
  nupfpdg.clear();
  nupftothits.clear();
  nusltothits.clear();
  nuslelhits.clear();
  nusltotprhits.clear();
  nusltotpihits.clear();
  nusltotnphits.clear();
  nusltotnehits.clear();
  nupfscore.clear();
  nuvtxx.clear();
  nuvtxy.clear();
  nuvtxz.clear();
  pfpdg.clear();
  pfhits.clear();
  pfelhits.clear();
  pftotprhits.clear();
  pftotpihits.clear();
  pftotnphits.clear();
  pftotnehits.clear();
  pfprhits.clear();
  pfpihits.clear();
  pfprid.clear();
  pfpiid.clear();
  pfnpr.clear();
  pfnpi.clear();
  pfscore.clear();
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
  tree->Branch("npie", &npie);
  tree->Branch("npip", &npip);
  tree->Branch("neue", &neue);
  tree->Branch("neup", &neup);
  //
  tree->Branch("tothits", &tothits, "tothits/I");
  tree->Branch("nelhits", &nelhits, "nelhits/I");
  tree->Branch("ntotprhits", &ntotprhits, "ntotprhits/I");
  tree->Branch("ntotpihits", &ntotpihits, "ntotpihits/I");
  tree->Branch("ntotnphits", &ntotnphits, "ntotnphits/I");
  tree->Branch("ntotnehits", &ntotnehits, "ntotnehits/I");
  tree->Branch("nprhits", &nprhits);
  tree->Branch("npihits", &npihits);
  //
  tree->Branch("tothitscr", &tothitscr, "tothitscr/I");
  tree->Branch("nelhitscr", &nelhitscr, "nelhitscr/I");
  tree->Branch("ntotprhitscr", &ntotprhitscr, "ntotprhitscr/I");
  tree->Branch("ntotpihitscr", &ntotpihitscr, "ntotpihitscr/I");
  tree->Branch("ntotnphitscr", &ntotnphitscr, "ntotnphitscr/I");
  tree->Branch("ntotnehitscr", &ntotnehitscr, "ntotnehitscr/I");
  tree->Branch("nprhitscr", &nprhitscr);
  tree->Branch("npihitscr", &npihitscr);
  //
  tree->Branch("pftothits", &pftothits, "pftothits/I");
  tree->Branch("npfnus", &npfnus, "npfnus/I");
  tree->Branch("nupfpdg", &nupfpdg);
  tree->Branch("nupfscore", &nupfscore);
  tree->Branch("nuvtxx", &nuvtxx);
  tree->Branch("nuvtxy", &nuvtxy);
  tree->Branch("nuvtxz", &nuvtxz);
  tree->Branch("nupftothits", &nupftothits);
  tree->Branch("nusltothits", &nusltothits);
  tree->Branch("nuslelhits", &nuslelhits);
  tree->Branch("nusltotprhits", &nusltotprhits);
  tree->Branch("nusltotpihits", &nusltotpihits);
  tree->Branch("nusltotnphits", &nusltotnphits);
  tree->Branch("nusltotnehits", &nusltotnehits);
  tree->Branch("pfpdg", &pfpdg);
  tree->Branch("pfhits", &pfhits);
  tree->Branch("pfelhits", &pfelhits);
  tree->Branch("pftotprhits", &pftotprhits);
  tree->Branch("pftotpihits", &pftotpihits);
  tree->Branch("pftotnphits", &pftotnphits);
  tree->Branch("pftotnehits", &pftotnehits);
  tree->Branch("pfprhits", &pfprhits);
  tree->Branch("pfpihits", &pfpihits);
  tree->Branch("pfprid", &pfprid);
  tree->Branch("pfpiid", &pfpiid);
  tree->Branch("pfnpr", &pfnpr);
  tree->Branch("pfnpi", &pfnpi);
  tree->Branch("pfscore", &pfscore);
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
  std::vector<double> enpi;
  std::vector<double> eneu;
  std::vector<double> ppro;
  std::vector<double> ppio;
  std::vector<double> pnpi;
  std::vector<double> pneu;

  detinfo::DetectorProperties const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
  detinfo::DetectorClocks const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>(); 

  // generator, simb::MCTruth
  art::InputTag GenInputTag("generator");
  art::ValidHandle<std::vector<simb::MCTruth> > inputMCTruth = e.getValidHandle<std::vector<simb::MCTruth> >(GenInputTag);
  //std::cout << "N mct=" << inputMCTruth->size() << std::endl;
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
      if (std::abs(mcp.PdgCode())==111) {
	enpi.push_back(mcp.E());
	pnpi.push_back(mcp.P());
      }
      if (std::abs(mcp.PdgCode())==2212) {
	epro.push_back(mcp.E());
	ppro.push_back(mcp.P());
        // std::cout << "id=" << mcp.TrackId() << " status=" << mcp.StatusCode() << " pdg=" << mcp.PdgCode() << " mother= " << mcp.Mother() << " E=" << mcp.E() << " P=" << mcp.P() << std::endl;//HERE
      }
      if (std::abs(mcp.PdgCode())==2112) {
	eneu.push_back(mcp.E());
	pneu.push_back(mcp.P());
      }
    }
  }
  if (abs(pdglep)!=11 && abs(pdglep)!=13) return;

  if (!isinfidvol) return;

  // fill tree variables with MCTruth info
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
  npie = enpi;
  npip = pnpi;
  neue = eneu;
  neup = pneu;
  //

  //need to find the id of the primary MCParticles in the right collection to be matched with the backtracker map
  int eleid = -1;
  std::vector<int> swrid;//needed to capture particles from the electron shower
  std::vector<int> proid;
  std::vector<int> pioid;
  std::vector<int> neuid;
  std::vector<int> npiid;
  // largeant, simb::MCParticle
  art::InputTag MCInputTag("largeant");
  art::ValidHandle<std::vector<simb::MCParticle> > inputMCParticle = e.getValidHandle<std::vector<simb::MCParticle> >(MCInputTag);
  /// std::cout << "N mcp=" << inputMCParticle->size() << std::endl;
  for ( unsigned int im=0; im<inputMCParticle->size(); ++im ) {
    const auto& mcp = inputMCParticle->at(im);
    if (mcp.StatusCode()==1) {
      // std::cout << "id=" << mcp.TrackId() << " status=" << mcp.StatusCode()<< " pdg=" << mcp.PdgCode() << " mother= " << mcp.Mother() << " E=" << mcp.E() << std::endl;
      // if (abs(mcp.PdgCode())==2212) std::cout << "id=" << mcp.TrackId() << " status=" << mcp.StatusCode()<< " pdg=" << mcp.PdgCode() << " mother= " << mcp.Mother() << " E=" << mcp.E() << std::endl;//HERE
      if (mcp.Mother()==0 && mcp.StatusCode()==1) {
	// /// std::cout << "id=" << mcp.TrackId() << " status=" << mcp.StatusCode()<< " pdg=" << mcp.PdgCode() << " mother= " << mcp.Mother() << " E=" << mcp.E() << std::endl;
	if (abs(mcp.PdgCode())==pdglep && std::abs(elep-mcp.E())<0.00001) {
	  eleid = mcp.TrackId();
	  swrid.push_back(eleid);
	} else if (abs(mcp.PdgCode())==2212) {
	  for (unsigned int ip=0; ip<ppro.size(); ip++) {
	    if (std::abs(ppro[ip]-mcp.P())<0.00000001) {
	      // std::cout << "ip=" << ip << " im=" << im << " epro[ip]=" << epro[ip] << " mcp.E()=" << mcp.E() << " mcp.P()=" << mcp.P() << " mother= " << mcp.Mother() << std::endl;//HERE
	      proid.push_back(mcp.TrackId());

	    }
	  }
	} else if (abs(mcp.PdgCode())==211) {
	  for (unsigned int ip=0; ip<ppio.size(); ip++) {
	    if (std::abs(ppio[ip]-mcp.P())<0.00000001) pioid.push_back(mcp.TrackId());
	  }
	} else if (abs(mcp.PdgCode())==111) {
	  for (unsigned int ip=0; ip<pnpi.size(); ip++) {
	    if (std::abs(pnpi[ip]-mcp.P())<0.00000001) npiid.push_back(mcp.TrackId());
	  }
	} else if (abs(mcp.PdgCode())==2112) {
	  for (unsigned int ip=0; ip<pneu.size(); ip++) {
	    if (std::abs(pneu[ip]-mcp.P())<0.00000001) neuid.push_back(mcp.TrackId());
	  }
	}
      }
      // navigate shower history
      if (std::find(swrid.begin(), swrid.end(), mcp.Mother())!=swrid.end()) {
	swrid.push_back(mcp.TrackId());
	// std::cout << "adding swr particle with mother=" << mcp.Mother() << " and id=" << mcp.TrackId() << std::endl;
      }
      if (std::find(npiid.begin(), npiid.end(), mcp.Mother())!=npiid.end()) {
	npiid.push_back(mcp.TrackId());
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
  int nnpihits = 0;
  int nneuhits = 0;
  for ( unsigned int ih=0; ih<inputHits->size(); ih++) {
    auto assmcp = assocMCPart->at(ih);
    auto assmdt = assocMCPart->data(ih);
    // const auto& hit = inputHits->at(ih);
    // if (hit.Channel()==1038)  
    // if (assmcp.size()>0)
    // std::cout << "gh hit id=" << ih << " view=" << hit.View() << " channel=" << hit.Channel() << " time[s,mr,p,pr,e]=" << hit.StartTick() << ", " << hit.PeakTimeMinusRMS() << ", " << hit.PeakTime() << ", " << hit.PeakTimePlusRMS() << ", " << hit.EndTick() 
    // 	      << " rms=" << hit.RMS() << " mult=" << hit.Multiplicity() << " nAssocMCP=" << assmcp.size() << std::endl;
    for (unsigned int ia=0; ia<assmcp.size(); ++ia) {
      auto mcp = assmcp[ia];
      auto amd = assmdt[ia];
      if (amd->isMaxIDE!=1) continue;
      // if (hit.Channel()==1038 && mcp->StatusCode()==1)  /// std::cout << "mcp#" << ia << "id=" << mcp->TrackId() << " pdg=" << mcp->PdgCode() << " mother= " << mcp->Mother() << " E=" << mcp->E()  << std::endl;
      // if ( hit.Channel()==3818 /*mcp->PdgCode()==2212*/) {
      // if (mcp->PdgCode()==13)
      // std::cout << "mcp#" << ia << " id=" << mcp->TrackId() << " view=" << hit.View() << " channel=" << hit.Channel() << " ideFrac=" << amd->ideFraction << " max=" << amd->isMaxIDE << " numEl=" << amd->numElectrons << " e=" << amd->energy << " pdg=" << mcp->PdgCode() << " mother= " << mcp->Mother() << " E=" << mcp->E()  << std::endl;
      // 	std::cout << "gh hit id=" << ih << " channel=" << hit.Channel() << " time[s,p,e]=" << hit.StartTick() << ", " << hit.PeakTime() << ", " << hit.EndTick() 
      // 		  << " rms=" << hit.RMS() << " mult=" << hit.Multiplicity() << " nAssocMCP=" << assmcp.size() << std::endl;
      // }
      // if (eleid == mcp->TrackId()) {
      if (std::find(swrid.begin(),swrid.end(),mcp->TrackId())!=swrid.end()) {
	// const auto& hit = inputHits->at(ih);
	//  /// std::cout << "hit id=" << ih << " channel=" << hit.Channel() << " time=" << hit.PeakTime() << " nAssocMCP=" << assmcp.size() << std::endl;
	//  /// std::cout << "id=" << mcp->TrackId() << " pdg=" << mcp->PdgCode() << " mother= " << mcp->Mother() << " E=" << mcp->E()  << std::endl;
	nelehits++;
	//break;//not needed with isMaxIDE
      } else if (std::find(proid.begin(),proid.end(),mcp->TrackId())!=proid.end()) {
	nprohits[std::find(proid.begin(),proid.end(),mcp->TrackId())-proid.begin()]++;
	//break;//not needed with isMaxIDE
      } else if (std::find(pioid.begin(),pioid.end(),mcp->TrackId())!=pioid.end()) {
	npiohits[std::find(pioid.begin(),pioid.end(),mcp->TrackId())-pioid.begin()]++;
	//break;//not needed with isMaxIDE
      } else if (std::find(npiid.begin(),npiid.end(),mcp->TrackId())!=npiid.end()) {
	nnpihits++;
	//break;//not needed with isMaxIDE
      } else if (std::find(neuid.begin(),neuid.end(),mcp->TrackId())!=neuid.end()) {
	nneuhits++;
	//break;//not needed with isMaxIDE
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
  ntotnphits = nnpihits;
  ntotnehits = nneuhits;
  nprhits = nprohits;
  npihits = npiohits;
  //

  // pandora PFParticles
  //
  int mytothitscr = 0;
  int mynelhitscr = 0;
  int myntotprhitscr = 0;
  int myntotpihitscr = 0;
  int myntotnphitscr = 0;
  int myntotnehitscr = 0;
  std::vector<int> mynprhitscr(ppro.size(),0);
  std::vector<int> mynpihitscr(ppio.size(),0);
  //
  int mypftothits = 0;
  int mynpfnus = 0;
  std::vector<int> mynupfpdg;
  std::vector<int> mynupftothits;
  std::vector<int> mynusltothits;
  std::vector<int> mynuslelhits;
  std::vector<int> mynusltotprhits;
  std::vector<int> mynusltotpihits;
  std::vector<int> mynusltotnphits;
  std::vector<int> mynusltotnehits;
  std::vector<float> mynupfscore;
  std::vector<float> mynuvtxx;
  std::vector<float> mynuvtxy;
  std::vector<float> mynuvtxz;
  std::vector< std::vector<int> > mypfpdg;
  std::vector< std::vector<int> > mypfhits;
  std::vector< std::vector<int> > mypfelhits;
  std::vector< std::vector<int> > mypftotprhits;
  std::vector< std::vector<int> > mypftotpihits;
  std::vector< std::vector<int> > mypftotnphits;
  std::vector< std::vector<int> > mypftotnehits;
  std::vector< std::vector<int> > mypfprhits;
  std::vector< std::vector<int> > mypfpihits;
  std::vector< std::vector<int> > mypfprid;
  std::vector< std::vector<int> > mypfpiid;
  std::vector< std::vector<int> > mypfnpr;
  std::vector< std::vector<int> > mypfnpi;
  std::vector< std::vector<float> > mypfscore;
  std::vector< std::vector<float> > mypfvtxx;
  std::vector< std::vector<float> > mypfvtxy;
  std::vector< std::vector<float> > mypfvtxz;
  //
  art::InputTag PfInputTag("pandora");
  art::ValidHandle<std::vector<recob::PFParticle> > inputPfParticle = e.getValidHandle<std::vector<recob::PFParticle> >(PfInputTag);
  art::ValidHandle<std::vector<recob::Cluster> > inputCluster = e.getValidHandle<std::vector<recob::Cluster> >(PfInputTag);
  auto assocCluster = std::unique_ptr<art::FindManyP<recob::Cluster> >(new art::FindManyP<recob::Cluster>(inputPfParticle, e, PfInputTag));
  auto assocHit = std::unique_ptr<art::FindManyP<recob::Hit> >(new art::FindManyP<recob::Hit>(inputCluster, e, PfInputTag));
  auto assocVertex = std::unique_ptr<art::FindManyP<recob::Vertex> >(new art::FindManyP<recob::Vertex>(inputPfParticle, e, PfInputTag));
  art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(inputPfParticle, e, PfInputTag);

  auto assocSlice = std::unique_ptr<art::FindManyP<recob::Slice> >(new art::FindManyP<recob::Slice>(inputPfParticle, e, PfInputTag));
  art::ValidHandle<std::vector<recob::Slice> > inputSlice = e.getValidHandle<std::vector<recob::Slice> >(PfInputTag);
  auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit> >(new art::FindManyP<recob::Hit>(inputSlice, e, PfInputTag));

  /// std::cout << "N pfps=" << inputPfParticle->size() << std::endl;
  for (unsigned int inpf=0; inpf<inputPfParticle->size(); ++inpf) {
    art::Ptr<recob::PFParticle> npfp(inputPfParticle,inpf);
    // std::cout << "inpf=" << inpf << " pdg=" << npfp->PdgCode() << " primary=" << npfp->IsPrimary() << " parent=" << npfp->Parent() << std::endl;
    //
    /*
      Here's the structure of these Metadata
      Primary PfParticles are either
      1) IsClearCosmic = 1 (for unambiguous cosmics)
      2) NuScore = 0.108586, SliceIndex = 1 (for cosmic slices)
      3) IsNeutrino = 1, NuScore = 0.170914, SliceIndex = 2 (for the nu slice)
      Then, for PfParticles that are daughter of the nu, the track score is saved, e.g.:
      4) TrackScore = 0.671488
      PfParticles that are not primary and that are not daugthers of the neutrino have empty Metadata
    */
    bool isClearCosmic = false;
    bool isTheNeutrino = false;
    double nuScore = -1;
    const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(npfp.key()));
    if (!pfParticleMetadataList.empty()) {
      for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j) {
	const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
	const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
	for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	  if (it->first=="IsClearCosmic" && it->second==1) isClearCosmic = true;
	  if (it->first=="NuScore") nuScore = it->second;
	  if (it->first=="IsNeutrino" && it->second==1) isTheNeutrino = true;
	}
      }
    }
    if (npfp->IsPrimary()==false) continue; 
    //std::cout << "Primary pfp with isClearCosmic=" << isClearCosmic << " isTheNeutrino=" << isTheNeutrino << " nuScore=" << nuScore << std::endl;

    //  recob::Slices are unique for the neutrino slice and for clear cosmics. Ambiguous cosmics may have >1 primary PFParticle per slice
    // (i.e. the slice is not unique in terms of primary PFParticles)
    auto slices = assocSlice->at(npfp.key());
    if (slices.size()!=1) {
      std::cout << "WRONG!!! n slices = " << slices.size() << std::endl;
    }
    auto slicehit = assocSliceHit->at(slices[0].key());

    if (isClearCosmic) {
      // this is the vector of pfps in this slice 
      std::vector<size_t> spfps = pfpDaughtersKey(npfp,inputPfParticle);//initialize with the daughters
      // include also pfps one more layer down in the hierarchy (secondaries)
      for (unsigned int ipf : spfps) {
	art::Ptr<recob::PFParticle> ppfp(inputPfParticle,ipf);
	const std::vector<size_t>& dpfps_tmp = pfpDaughtersKey(ppfp,inputPfParticle);
	for (auto k : dpfps_tmp) spfps.push_back(k);
      }
      //finally add the primary as well, since it's not the neutrino
      spfps.push_back(npfp.key());
      int npfhits = 0;
      int nelehitsCR = 0;
      std::vector<int> nprohitsCR(ppro.size(),0);
      std::vector<int> npiohitsCR(ppio.size(),0);
      int nnpihitsCR = 0;
      int nneuhitsCR = 0;
      for (unsigned int ipf : spfps) {
	art::Ptr<recob::PFParticle> ppfp(inputPfParticle,ipf);
	incrementCounts(ppfp, assocCluster, assocHit, assocMCPart, swrid, proid, pioid, npiid, neuid, npfhits, nelehitsCR, nprohitsCR, npiohitsCR, nnpihitsCR, nneuhitsCR);
      }
      //
      int countprohitsCR = 0;
      for (auto n : nprohitsCR) countprohitsCR+=n;
      int countpiohitsCR = 0;
      for (auto n : npiohitsCR) countpiohitsCR+=n;
      /// std::cout << "N hits after CR electron=" << nelehitsCR << " proton=" << countprohitsCR << " pion=" << countpiohitsCR << std::endl;
      mytothitscr += npfhits;
      mynelhitscr += nelehitsCR;
      myntotprhitscr += countprohitsCR;
      myntotpihitscr += countpiohitsCR;
      myntotnphitscr += nnpihitsCR;
      myntotnehitscr += nneuhitsCR;
      for (unsigned int k=0;k<nprohitsCR.size();k++) mynprhitscr[k] += nprohitsCR[k];
      for (unsigned int k=0;k<npiohitsCR.size();k++) mynpihitscr[k] += npiohitsCR[k];
      continue;
    } 

    mynpfnus++;
    mynupfpdg.push_back(npfp->PdgCode());
    mynupfscore.push_back(nuScore);
    mynupftothits.push_back(0);
    int nslelhits=0;
    int nsltotprhits=0;
    int nsltotpihits=0;
    int nsltotnphits=0;
    int nsltotnehits=0;
    countSliceHits(slicehit,assocMCPart,swrid,proid,pioid,npiid,neuid,nslelhits,nsltotprhits,nsltotpihits,nsltotnphits,nsltotnehits);
    mynusltothits.push_back(slicehit.size());
    mynuslelhits.push_back(nslelhits);
    mynusltotprhits.push_back(nsltotprhits);
    mynusltotpihits.push_back(nsltotpihits);
    mynusltotnphits.push_back(nsltotnphits);
    mynusltotnehits.push_back(nsltotnehits);
    const std::vector<art::Ptr<recob::Vertex> >& vertexVec = assocVertex->at(npfp.key());
    if (vertexVec.size()>0){
      double xyz[3];
      vertexVec.front()->XYZ(xyz);
      mynuvtxx.push_back(xyz[0]);
      mynuvtxy.push_back(xyz[1]);
      mynuvtxz.push_back(xyz[2]);
    } else {
      mynuvtxx.push_back(-999.);
      mynuvtxy.push_back(-999.);
      mynuvtxz.push_back(-999.);
    }
    mypfpdg.push_back( std::vector<int>() );
    mypfhits.push_back( std::vector<int>() );
    mypfelhits.push_back( std::vector<int>() );
    mypftotprhits.push_back( std::vector<int>() );
    mypftotpihits.push_back( std::vector<int>() );
    mypftotnphits.push_back( std::vector<int>() );
    mypftotnehits.push_back( std::vector<int>() );
    mypfprhits.push_back( std::vector<int>() );
    mypfpihits.push_back( std::vector<int>() );
    mypfprid.push_back( std::vector<int>() );
    mypfpiid.push_back( std::vector<int>() );
    mypfnpr.push_back( std::vector<int>() );
    mypfnpi.push_back( std::vector<int>() );
    mypfscore.push_back( std::vector<float>() );
    mypfvtxx.push_back( std::vector<float>() );
    mypfvtxy.push_back( std::vector<float>() );
    mypfvtxz.push_back( std::vector<float>() );
    /// std::cout << "neutrino pfp pdg=" << npfp->PdgCode() << " and nDaughters=" << npfp->NumDaughters() << std::endl;
    std::vector<size_t> dpfps = pfpDaughtersKey(npfp,inputPfParticle);
    // include also pfps one more layer down in the hierarchy (secondaries)
    for (unsigned int ipf : dpfps) {
      art::Ptr<recob::PFParticle> ppfp(inputPfParticle,ipf);
      const std::vector<size_t>& dpfps_tmp = pfpDaughtersKey(ppfp,inputPfParticle);
      for (auto k : dpfps_tmp) dpfps.push_back(k);
    }
    // if it's not THE neutrino, need to add the primary as well
    if (isTheNeutrino==false) dpfps.push_back(npfp.key());
    // now let's go over primaries and secondaries
    for (unsigned int ipf : dpfps) {
      int npfhits = 0;
      int nelehitsPF = 0;
      std::vector<int> nprohitsPF(ppro.size(),0);
      std::vector<int> npiohitsPF(ppio.size(),0);
      int nnpihitsPF = 0;
      int nneuhitsPF = 0;
      art::Ptr<recob::PFParticle> ppfp(inputPfParticle,ipf);
      incrementCounts(ppfp, assocCluster, assocHit, assocMCPart, swrid, proid, pioid, npiid, neuid, npfhits, nelehitsPF, nprohitsPF, npiohitsPF, nnpihitsPF, nneuhitsPF);
      //
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
      mypfpdg.back().push_back(ppfp->PdgCode());
      mypfhits.back().push_back(npfhits);
      mypfelhits.back().push_back(nelehitsPF);
      mypftotprhits.back().push_back(countprohitsPF);
      mypftotpihits.back().push_back(countpiohitsPF);
      mypftotnphits.back().push_back(nnpihitsPF);
      mypftotnehits.back().push_back(nneuhitsPF);
      mypfnpr.back().push_back(inprn);
      mypfnpi.back().push_back(inpin);
      //
      double trackScore = -1;
      const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(ipf));
      if (!pfParticleMetadataList.empty()) {
	for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j) {
	  const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
	  const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
	  for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	    if (it->first=="TrackScore") trackScore = it->second;
	  }
	}
      }
      mypfscore.back().push_back(trackScore);
      //
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
  tothitscr = mytothitscr;
  nelhitscr = mynelhitscr;
  ntotprhitscr = myntotprhitscr;
  ntotpihitscr = myntotpihitscr;
  ntotnphitscr = myntotnphitscr;
  ntotnehitscr = myntotnehitscr;
  nprhitscr = mynprhitscr;
  npihitscr = mynpihitscr;
  //
  pftothits = mypftothits;
  npfnus = mynpfnus;
  nupfpdg = mynupfpdg;
  nupfscore = mynupfscore;
  nuvtxx = mynuvtxx;
  nuvtxy = mynuvtxy;
  nuvtxz = mynuvtxz;
  nupftothits = mynupftothits;
  nusltothits = mynusltothits;
  nuslelhits = mynuslelhits;
  nusltotprhits = mynusltotprhits;
  nusltotpihits = mynusltotpihits;
  nusltotnphits = mynusltotnphits;
  nusltotnehits = mynusltotnehits;
  pfpdg = mypfpdg;
  pfhits = mypfhits;
  pfelhits = mypfelhits;
  pftotprhits = mypftotprhits;
  pftotpihits = mypftotpihits;
  pftotnphits = mypftotnphits;
  pftotnehits = mypftotnehits;
  pfprhits = mypfprhits;
  pfpihits = mypfpihits;
  pfprid = mypfprid;
  pfpiid = mypfpiid;
  pfnpr = mypfnpr;
  pfnpi = mypfnpi;
  pfscore = mypfscore;
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

void LeeSimRecoTest::incrementCounts(const art::Ptr<recob::PFParticle> ppfp,
				     const std::unique_ptr<art::FindManyP<recob::Cluster> >& assocCluster,
				     const std::unique_ptr<art::FindManyP<recob::Hit> >& assocHit,
				     const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& assocMCPart,
				     const std::vector<int>& swrid, const std::vector<int>& proid, const std::vector<int>& pioid, const std::vector<int>& npiid, const std::vector<int>& neuid,
				     int& npfhits, int& nelehitsPF, std::vector<int>& nprohitsPF, std::vector<int>& npiohitsPF, int& nnpihitsPF, int& nneuhitsPF) const
{
  const std::vector<art::Ptr<recob::Cluster> >& clustVec = assocCluster->at(ppfp.key());
  for (unsigned int icl=0; icl<clustVec.size(); icl++) {
    art::Ptr<recob::Cluster> cl = clustVec[icl];
    const std::vector<art::Ptr<recob::Hit> >& hitVec = assocHit->at(cl.key());
    npfhits += hitVec.size();
    for (unsigned int ih=0; ih<hitVec.size(); ih++) {
      art::Ptr<recob::Hit> hitp = hitVec[ih];
      auto assmcp = assocMCPart->at(hitp.key());
      auto assmdt = assocMCPart->data(hitp.key());
      for (unsigned int ia=0; ia<assmcp.size(); ++ia) {
	auto mcp = assmcp[ia];
	auto amd = assmdt[ia];
	if (amd->isMaxIDE!=1) continue;
	if (std::find(swrid.begin(),swrid.end(),mcp->TrackId())!=swrid.end()) {
	  nelehitsPF++;
	  break;
	} else if (std::find(proid.begin(),proid.end(),mcp->TrackId())!=proid.end()) {
	  bool test = ( int(nprohitsPF.size()) > int(std::find(proid.begin(),proid.end(),mcp->TrackId())-proid.begin()) );
	  if (!test) {
	    std::cout << __FILE__ << " " << __LINE__ << std::endl;
	    exit(1);
	  }
	  nprohitsPF[std::find(proid.begin(),proid.end(),mcp->TrackId())-proid.begin()]++;
	  break;
	} else if (std::find(pioid.begin(),pioid.end(),mcp->TrackId())!=pioid.end()) {
	  bool test = ( int(npiohitsPF.size()) > int(std::find(pioid.begin(),pioid.end(),mcp->TrackId())-pioid.begin()) );
	  if (!test) {
	    std::cout << __FILE__ << " " << __LINE__ << std::endl;
	    exit(1);
	  }
	  npiohitsPF[std::find(pioid.begin(),pioid.end(),mcp->TrackId())-pioid.begin()]++;
	  break;
	} else if (std::find(npiid.begin(),npiid.end(),mcp->TrackId())!=npiid.end()) {
	  nnpihitsPF++;
	  break;
	} else if (std::find(neuid.begin(),neuid.end(),mcp->TrackId())!=neuid.end()) {
	  nneuhitsPF++;
	  break;
	}
      }
    }
  }
}

void LeeSimRecoTest::countSliceHits(const std::vector<art::Ptr<recob::Hit> >& hits,
				    const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& assocMCPart,
				    const std::vector<int>& swrid, const std::vector<int>& proid, const std::vector<int>& pioid, const std::vector<int>& npiid, const std::vector<int>& neuid,
				    int& nelehitsPF, int& nprohitsPF, int& npiohitsPF, int& nnpihitsPF, int& nneuhitsPF) const
{
  for (unsigned int ih=0; ih<hits.size(); ih++) {
    art::Ptr<recob::Hit> hitp = hits[ih];
    auto assmcp = assocMCPart->at(hitp.key());
    auto assmdt = assocMCPart->data(hitp.key());
    for (unsigned int ia=0; ia<assmcp.size(); ++ia) {
      auto mcp = assmcp[ia];
      auto amd = assmdt[ia];
      if (amd->isMaxIDE!=1) continue;
      if (std::find(swrid.begin(),swrid.end(),mcp->TrackId())!=swrid.end()) {
	nelehitsPF++;
	break;
      } else if (std::find(proid.begin(),proid.end(),mcp->TrackId())!=proid.end()) {
	nprohitsPF++;
	break;
      } else if (std::find(pioid.begin(),pioid.end(),mcp->TrackId())!=pioid.end()) {
	npiohitsPF++;
	break;
      } else if (std::find(npiid.begin(),npiid.end(),mcp->TrackId())!=npiid.end()) {
	nnpihitsPF++;
      } else if (std::find(neuid.begin(),neuid.end(),mcp->TrackId())!=neuid.end()) {
	nneuhitsPF++;
	break;
      }
    }
  }
}

std::vector<size_t> LeeSimRecoTest::pfpDaughtersKey(const art::Ptr<recob::PFParticle> ppfp,
						    art::ValidHandle<std::vector<recob::PFParticle> > inputPfParticle) {
  std::vector<size_t> result;
  // Daugthers returns the id as in pfp->Self() not the key
  // so need to find the key for the pfp with Self==ipfd
  auto& pfd = ppfp->Daughters();
  for (auto ipfd : pfd) {
    for (size_t jPF = 0; jPF < inputPfParticle->size(); ++jPF) {
      art::Ptr<recob::PFParticle> pfpd(inputPfParticle, jPF);
      if (pfpd->Self()!=ipfd) continue;
      result.push_back(pfpd.key());
      break;
    }
  }
  return result;
}

DEFINE_ART_MODULE(LeeSimRecoTest)
