////////////////////////////////////////////////////////////////////////
// Class:       VectorHitFinder
// Plugin Type: analyzer (art v2_09_06)
// File:        VectorHitFinder_module.cc
//
// Generated at Mon Apr  9 15:11:04 2018 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_01_03.
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
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/TrackingPlaneHelper.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/RecoObjects/TrackStatePropagator.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/MCBase/MCHitCollection.h"

#include <iterator>
#include <iostream>
#include <fstream>

#include "Event.h"

class VectorHitFinder;

class VectorHitFinder : public art::EDAnalyzer {
public:
  explicit VectorHitFinder(fhicl::ParameterSet const & p);
  ~VectorHitFinder();
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  VectorHitFinder(VectorHitFinder const &) = delete;
  VectorHitFinder(VectorHitFinder &&) = delete;
  VectorHitFinder & operator = (VectorHitFinder const &) = delete;
  VectorHitFinder & operator = (VectorHitFinder &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.
  art::ServiceHandle<geo::Geometry> geom;
  const detinfo::DetectorProperties* detprop;
  std::ofstream outfile;
  gshf::DataFile data_file;
  long long int savedEvents = 0;
};


VectorHitFinder::VectorHitFinder(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
  outfile.open("hitfinder.txt");
  data_file.OpenWrite("hitfinder.bin", 1);
}

VectorHitFinder::~VectorHitFinder()
{
  outfile.close();
  data_file.CloseWrite(savedEvents);
  printf("\nSaved %lli events\n\n",savedEvents);
}

void VectorHitFinder::analyze(art::Event const & e) {

  gshf::Event EE(e.event());
  std::vector<gshf::wiredata> &wd_vec_ = EE.wd_vec_;
  std::vector<gshf::refdata>  &rd_vec_ = EE.rd_vec_;

  detinfo::DetectorClocks const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  trkf::TrackStatePropagator prop(1.0, 0.1, 10, 10., 0.01, false);

  const auto& mcps = e.getValidHandle<std::vector<simb::MCParticle> >("largeant");
  if (0) std::cout << "nmcp=" << mcps->size() << std::endl;
  const simb::MCParticle* m = nullptr;
  for (const auto& mcp : (*mcps)) {
    if (mcp.Mother()!=0) continue;
    if (0) std::cout << "mcp np=" << mcp.NumberTrajectoryPoints()
	      << " status=" << mcp.StatusCode()
	      << " pdg=" << mcp.PdgCode()
	      << " mother=" << mcp.Mother()
	      << std::endl;
    m = &mcp;
    break;
  }

  const auto& mcts = e.getValidHandle<std::vector<sim::MCTrack> >("mcreco");
  for (const auto& mct : (*mcts)) {
    if (0) std::cout << "mct np=" << mct.size()
	      << " pdg=" << mct.PdgCode()
	      << std::endl;
  }
  
  auto const& hitsp = proxy::getCollection<std::vector<recob::Hit> >(e,"gaushit",proxy::withAssociated<recob::Wire>(),proxy::withAssociated<raw::RawDigit>(),proxy::withAssociatedMeta<simb::MCParticle,anab::BackTrackerHitMatchingData>(art::InputTag("pandoraNuTruthMatch")));
  size_t maxwf = 0;
  for (const auto& h : hitsp) {
    double t = h->PeakTime();
    double x = detprop->ConvertTicksToX(t, h->WireID().Plane, h->WireID().TPC, h->WireID().Cryostat);
    //
    /*
    std::cout << "hit view=" << h->View() << " wire=" << h->Channel() << " mult=" << h->Multiplicity() << " start=" << h->StartTick() << " end=" << h->EndTick()
	      << " peakTick=" << h->PeakTime() << " ampl=" << h->PeakAmplitude() << " RMS=" << h->RMS() << " sigma=" << h->SigmaPeakTime()
	      << " x=" << x
	      << std::endl;
    */
    const auto& mps = h.get<simb::MCParticle>();
    size_t nmps = 0;
    bool foundmcp = false;
    for (const auto& im : mps) {
      const auto& bthmd = im.data();
      if (0) std::cout << "mcp pdg=" << im->PdgCode() << " ideFraction=" << bthmd.ideFraction << " isMaxIDE=" << bthmd.isMaxIDE
		<< " ideNFraction=" << bthmd.ideNFraction << " isMaxIDEN=" << bthmd.isMaxIDEN
		<< " numElectrons=" << bthmd.numElectrons << " energy=" << bthmd.energy
		<< std::endl;
      if (m == &*im && bthmd.isMaxIDE && bthmd.ideFraction>0.9) {
	foundmcp = true;
      }
      nmps++;
    }
    //std::cout << "foundmcp=" << foundmcp << std::endl;
    if (!foundmcp) continue;
    const auto& w = h.get<recob::Wire>().front();
    const auto& rd = h.get<raw::RawDigit>().front();
    const auto& sigv = w->Signal();
    const auto& radcs = rd->ADCs();
    size_t nonzeroadcs = 0;
    for (size_t t=0;t<sigv.size();++t) {
      if (std::abs(sigv[t])>0) {
	if (0) std::cout << "tick=" << t << " ADC=" << sigv[t] << std::endl;
	nonzeroadcs++;
      }
    }
    size_t nonzerorawadcs = 0;
    for (size_t t=0;t<radcs.size();++t) {
      if (radcs[t]>0) {
	if (0) std::cout << "tick=" << t << " ADC=" << radcs[t] << std::endl;
	nonzerorawadcs++;
      }
    }
    if (0) std::cout << "wire view=" << w->View() << " wire=" << w->Channel()
	      << " nsig=" << w->NSignal() << " nROI=" << w->SignalROI().size() << " nTicks=" << w->Signal().size()
	      << " nonZeroADCs=" << nonzeroadcs
	      << " nRawADCs=" << radcs.size() << " nonZeroRaw=" << nonzerorawadcs
	      << std::endl;
    auto plane = recob::tracking::makePlane( geom->WireIDToWireGeo(h->WireID()) );
    for (size_t ip=0;ip<m->NumberTrajectoryPoints()-1;++ip) {
      double g4Ticks0 = detClocks->TPCG4Time2Tick(m->T(ip))+detprop->GetXTicksOffset(h->WireID().planeID())-detprop->TriggerOffset();
      double xOffset0 = detprop->ConvertTicksToX(g4Ticks0, h->WireID().planeID());
      double g4Ticks1 = detClocks->TPCG4Time2Tick(m->T(ip+1))+detprop->GetXTicksOffset(h->WireID().planeID())-detprop->TriggerOffset();
      double xOffset1 = detprop->ConvertTicksToX(g4Ticks1, h->WireID().planeID());
      if (0) {
	std::cout << "m->T(ip)=" << m->T(ip) << std::endl;
	std::cout << "g4Time=" << detClocks->G4ToElecTime(m->T(ip)) << std::endl;
	std::cout << "detClocks->TPCG4Time2Tick(m->T(ip))=" << detClocks->TPCG4Time2Tick(m->T(ip)) << std::endl;
	std::cout << "detprop->TriggerOffset()=" << detprop->TriggerOffset() << std::endl;
	std::cout << "detClocks->TriggerOffsetTPC()=" << detClocks->TriggerOffsetTPC() << std::endl;
	std::cout << "detClocks->TriggerTime()=" << detClocks->TriggerTime() << std::endl;
	std::cout << "detprop->GetXTicksOffset(h->WireID().planeID())=" << detprop->GetXTicksOffset(h->WireID().planeID()) << std::endl;
	std::cout << "t0=" << detprop->GetXTicksOffset(h->WireID().planeID())-detprop->TriggerOffset() << std::endl;
      }
      // std::cout << "g4TOffset0=" << detClocks->TPCG4Time2Tick(m->T(ip)) << " g4TOffset1=" << detClocks->TPCG4Time2Tick(m->T(ip+1)) << " ticksTOffset=" << detprop->GetXTicksOffset(h->WireID().planeID())
      // 		<< " trigTOffset=" << detprop->TriggerOffset() << " totTOffset0=" << g4Ticks0 << " xOffset0=" << xOffset0 << " totTOffset1=" << g4Ticks1 << " xOffset1=" << xOffset1  << std::endl;
      recob::tracking::Point_t p0(m->Vx(ip)+xOffset0,m->Vy(ip),m->Vz(ip));
      recob::tracking::Point_t p1(m->Vx(ip+1)+xOffset1,m->Vy(ip+1),m->Vz(ip+1));
      recob::tracking::Vector_t d0(m->Px(ip)/m->P(ip),m->Py(ip)/m->P(ip),m->Pz(ip)/m->P(ip));
      recob::tracking::Vector_t d1(m->Px(ip+1)/m->P(ip+1),m->Py(ip+1)/m->P(ip+1),m->Pz(ip+1)/m->P(ip+1));
      //
      bool success = true;
      recob::tracking::Point_t px = prop.propagatedPosByDistance(p0,d0, prop.distanceToPlane(success, p0, d0, plane) );
      //
      auto simtick = detprop->ConvertXToTicks(px.X(),h->WireID().planeID());
      //
      if ( plane.direction().Dot(p0-plane.position())*plane.direction().Dot(p1-plane.position())<0. && simtick>=h->StartTick() && simtick<=h->EndTick() ) {
	if (0) std::cout << "hit view=" << h->View() << " wire=" << h->Channel() << " mult=" << h->Multiplicity() << " nmps=" << nmps
		  << " simX=" << px.X() << " ticksX=" << detprop->ConvertXToTicks(px.X(),h->WireID().planeID())
		  << " peakTick=" << h->PeakTime() << " RMS=" << h->RMS() << " sigma=" << h->SigmaPeakTime() << " recX=" << x
		  << " deltaX=" << x-px.X() << " deltaTicks=" << h->PeakTime()-detprop->ConvertXToTicks(px.X(),h->WireID().planeID())
		  << " nonZeroADCs=" << nonzeroadcs << " nADCs=" << h->EndTick()-h->StartTick() << " nRawADCs=" << radcs.size()
		  << " xOffset0=" << xOffset0
		  << std::endl;
      outfile << "hit view=" << h->View() << " wire=" << h->Channel() << " mult=" << h->Multiplicity() << " multMC=" << nmps
	      << " simX=" << px.X() << " simTick=" << simtick
	      << " recX=" << x << " recTick=" << h->PeakTime() << " RMS=" << h->RMS() << " sigma=" << h->SigmaPeakTime()
	      << " startTick=" << h->StartTick() << " endTick=" << h->EndTick()
	      << " deltaX=" << x-px.X() << " deltaTicks=" << h->PeakTime()-detprop->ConvertXToTicks(px.X(),h->WireID().planeID())
	      << " nRoiADCs=" << nonzeroadcs << " nHitADCs=" << h->EndTick()-h->StartTick() << " nRawADCs=" << radcs.size();
      gshf::refdata  rd;
      rd.simtck = simtick;
      rd.rectck = h->PeakTime();
      rd.rms = h->RMS();
      rd_vec_.push_back(rd);
      gshf::wiredata wd;
      wd.vw = h->View();
      wd.ntck = 0;
      try
	{
          // ####################################
          // ### Using BackTrackerService HitCheater ###
          // ####################################
          art::ServiceHandle<cheat::BackTrackerService> bt_serv;
          //trackides = bt_serv->HitToTrackIDEs(*h);
          std::vector<double> xyz = bt_serv->HitToXYZ(*h);
	  if (0) std::cout << "from backtracker x=" << xyz[0]+xOffset0 /*<< " y=" << xyz[1] << " z=" << xyz[2]*/ << std::endl;
	  outfile << " backtrackerX=" << xyz[0]+xOffset0;
	}
      catch(cet::exception e)
	{
	  mf::LogWarning("GausHitFinderAna") << "BackTrackerService Failed";
	  continue;
	}

      outfile << "\n";
      for (size_t t=0;t<sigv.size();++t) {
	if (std::abs(sigv[t])>0) {
	  outfile << "tick=" << t << " ADC=" << sigv[t] << "\n";
	  gshf::waveform wf;
	  wf.tck = t;
	  wf.adc = sigv[t];
	  wd.wv.push_back(wf);
	  wd.ntck++;
	}
      }
      wd_vec_.push_back(wd);
      if (maxwf <= wd.ntck) maxwf = wd.ntck;
      /*
	std::cout << "\tplane pos=" << plane.position() << " mc p0=" << p0 << " p1=" << p1 << " px=" << px
	<< " d0=" << d0 << " d1=" << d1
	<< " ticks0=" << detprop->ConvertXToTicks(p0.X(),h->WireID().planeID()) //<< " xOffset0=" << xOffset0
	<< " ticks1=" << detprop->ConvertXToTicks(p1.X(),h->WireID().planeID())
	<< " ticksX=" << detprop->ConvertXToTicks(px.X(),h->WireID().planeID())
	// << " dp0=" << p0-plane.position() << " dp1=" << p1-plane.position()
	// << " dot=" << plane.direction().Dot(p0-plane.position())*plane.direction().Dot(p1-plane.position())
	<< std::endl;
      */
      break;
      }
    }
  }
  
  // const auto& wires = e.getValidHandle<std::vector<recob::Wire> >("caldata");
  // for (const auto& w : (*wires)) {
  //   std::cout << "wire view=" << w.View() << " wire=" << w.Channel()
  // 	      << " nsig=" << w.NSignal() << " nROI=" << w.SignalROI().size() << " nTicks=" << w.Signal().size()
  // 	      << std::endl;
  //   const auto& sigv = w.Signal();
  //   for (size_t t=0;t<sigv.size();++t) {
  //     if (sigv[t]>0) std::cout << "tick=" << t << " ADC=" << sigv[t] << std::endl;
  //   }
  // }

  std::cout << "writing event with wd_vec.size=" << EE.wd_vec_.size() << " and rd_vec.size=" << EE.rd_vec_.size() << " and mxwf=" << maxwf << std::endl;

  
  EE.write_out(data_file);
  savedEvents++;
}

DEFINE_ART_MODULE(VectorHitFinder)
