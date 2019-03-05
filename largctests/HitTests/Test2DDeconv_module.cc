////////////////////////////////////////////////////////////////////////
// Class:       Test2DDeconv
// Plugin Type: analyzer (art v2_11_03)
// File:        Test2DDeconv_module.cc
//
// Generated at Wed Aug 29 12:37:40 2018 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_03_01.
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

#include "lardata/RecoBaseProxy/Track.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TH1.h"
#include "TH2.h"

class Test2DDeconv;

class Test2DDeconv : public art::EDAnalyzer {
public:
  explicit Test2DDeconv(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Test2DDeconv(Test2DDeconv const &) = delete;
  Test2DDeconv(Test2DDeconv &&) = delete;
  Test2DDeconv & operator = (Test2DDeconv const &) = delete;
  Test2DDeconv & operator = (Test2DDeconv &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  void beginJob() override;
private:

  TH1F* h_ntracks;
  TH1F* h_ntrackhits;
  TH1F* h_ntracks50;
  TH1F* h_ntrackhits50;
  //
  TH1F* h_dt_3m10;
  TH1F* h_dx_3m10;
  //
  TH1F* h_dt_3m12;
  TH1F* h_dx_3m12;
  //
  TH1F* h_dt_3m02;
  TH1F* h_dx_3m02;
  //
  TH1F* h_dt_loose_2m12;
  TH1F* h_dx_loose_2m12;
  TH1F* h_dt_loose_2m02;
  TH1F* h_dx_loose_2m02;
  TH1F* h_dt_loose_2m10;
  TH1F* h_dx_loose_2m10;
  //
  TH1F* h_dt_tight_2m12;
  TH1F* h_dx_tight_2m12;
  TH1F* h_dt_tight_2m02;
  TH1F* h_dx_tight_2m02;
  TH1F* h_dt_tight_2m10;
  TH1F* h_dx_tight_2m10;
  //
  TH2F* h_dyz_vs_dt_2m12;
  TH2F* h_dyz_vs_dt_2m02;
  TH2F* h_dyz_vs_dt_2m10;
  TH2F* h_dyz_vs_dt_2m01;
  TH2F* h_dyz_vs_dt_2m20;
  //
  TH2F* h_wire0_vs_dt_2m10;
  TH2F* h_wire1_vs_dt_2m10;
  //
  TH2F* h_mult0_vs_dt_2m10;
  TH2F* h_mult1_vs_dt_2m10;
  //
  TH2F* h_theta_vs_dt_2m10;
  TH2F* h_phi_vs_dt_2m10;
  TH2F* h_thetax_vs_dt_2m10;
  TH2F* h_thetay_vs_dt_2m10;
};


Test2DDeconv::Test2DDeconv(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void Test2DDeconv::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  //
  h_ntracks = tfs->make<TH1F>("h_ntracks","h_ntracks",50, 0, 100);
  h_ntrackhits = tfs->make<TH1F>("h_ntrackhits","h_ntrackhits",50, 0, 2500);
  h_ntracks50 = tfs->make<TH1F>("h_ntracks50","h_ntracks50",50, 0, 50);
  h_ntrackhits50 = tfs->make<TH1F>("h_ntrackhits50","h_ntrackhits50",50, 0, 2500);
  //
  h_dt_3m10 = tfs->make<TH1F>("h_dt_3m10","h_dt_3m10", 100, -20, 20);
  h_dx_3m10 = tfs->make<TH1F>("h_dx_3m10","h_dx_3m10", 40, -0.5, 0.5);
  //
  h_dt_3m12 = tfs->make<TH1F>("h_dt_3m12","h_dt_3m12", 100, -20, 20);
  h_dx_3m12 = tfs->make<TH1F>("h_dx_3m12","h_dx_3m12", 40, -0.5, 0.5);
  //
  h_dt_3m02 = tfs->make<TH1F>("h_dt_3m02","h_dt_3m02", 100, -20, 20);
  h_dx_3m02 = tfs->make<TH1F>("h_dx_3m02","h_dx_3m02", 40, -0.5, 0.5);
  //
  h_dt_loose_2m12 = tfs->make<TH1F>("h_dt_loose_2m12","h_dt_loose_2m12", 40, -20, 20);
  h_dx_loose_2m12 = tfs->make<TH1F>("h_dx_loose_2m12","h_dx_loose_2m12", 40, -0.1, 0.1);
  h_dt_loose_2m02 = tfs->make<TH1F>("h_dt_loose_2m02","h_dt_loose_2m02", 40, -20, 20);
  h_dx_loose_2m02 = tfs->make<TH1F>("h_dx_loose_2m02","h_dx_loose_2m02", 40, -0.1, 0.1);
  h_dt_loose_2m10 = tfs->make<TH1F>("h_dt_loose_2m10","h_dt_loose_2m10", 40, -20, 20);
  h_dx_loose_2m10 = tfs->make<TH1F>("h_dx_loose_2m10","h_dx_loose_2m10", 40, -0.1, 0.1);
  //
  h_dt_tight_2m12 = tfs->make<TH1F>("h_dt_tight_2m12","h_dt_tight_2m12", 40, -20, 20);
  h_dx_tight_2m12 = tfs->make<TH1F>("h_dx_tight_2m12","h_dx_tight_2m12", 40, -0.1, 0.1);
  h_dt_tight_2m02 = tfs->make<TH1F>("h_dt_tight_2m02","h_dt_tight_2m02", 40, -20, 20);
  h_dx_tight_2m02 = tfs->make<TH1F>("h_dx_tight_2m02","h_dx_tight_2m02", 40, -0.1, 0.1);
  h_dt_tight_2m10 = tfs->make<TH1F>("h_dt_tight_2m10","h_dt_tight_2m10", 40, -20, 20);
  h_dx_tight_2m10 = tfs->make<TH1F>("h_dx_tight_2m10","h_dx_tight_2m10", 40, -0.1, 0.1);
  //
  h_dyz_vs_dt_2m12 = tfs->make<TH2F>("h_dyz_vs_dt_2m12","h_dyz_vs_dt_2m12", 40, -20, 20, 40, 0, 2);
  h_dyz_vs_dt_2m02 = tfs->make<TH2F>("h_dyz_vs_dt_2m02","h_dyz_vs_dt_2m02", 40, -20, 20, 40, 0, 2);
  h_dyz_vs_dt_2m10 = tfs->make<TH2F>("h_dyz_vs_dt_2m10","h_dyz_vs_dt_2m10", 40, -20, 20, 40, 0, 2);
  h_dyz_vs_dt_2m01 = tfs->make<TH2F>("h_dyz_vs_dt_2m01","h_dyz_vs_dt_2m01", 40, -20, 20, 40, 0, 2);
  h_dyz_vs_dt_2m20 = tfs->make<TH2F>("h_dyz_vs_dt_2m20","h_dyz_vs_dt_2m20", 40, -20, 20, 40, 0, 2);
  //
  h_wire0_vs_dt_2m10 = tfs->make<TH2F>("h_wire0_vs_dt_2m10","h_wire0_vs_dt_2m10", 40, -20, 20, 500, 0, 2500);
  h_wire1_vs_dt_2m10 = tfs->make<TH2F>("h_wire1_vs_dt_2m10","h_wire1_vs_dt_2m10", 40, -20, 20, 500, 0, 2500);
  //
  h_mult0_vs_dt_2m10 = tfs->make<TH2F>("h_mult0_vs_dt_2m10","h_mult0_vs_dt_2m10", 40, -20, 20, 30, 0, 30);
  h_mult1_vs_dt_2m10 = tfs->make<TH2F>("h_mult1_vs_dt_2m10","h_mult1_vs_dt_2m10", 40, -20, 20, 30, 0, 30);
  //
  h_theta_vs_dt_2m10 = tfs->make<TH2F>("h_theta_vs_dt_2m10","h_theta_vs_dt_2m10", 40, -20, 20, 100, 0, 3.15);
  h_phi_vs_dt_2m10 = tfs->make<TH2F>("h_phi_vs_dt_2m10","h_phi_vs_dt_2m10", 40, -20, 20, 50, -3.15, 0);
  h_thetax_vs_dt_2m10 = tfs->make<TH2F>("h_thetax_vs_dt_2m10","h_thetax_vs_dt_2m10", 40, -20, 20, 100, -3.15, 3.15);
  h_thetay_vs_dt_2m10 = tfs->make<TH2F>("h_thetay_vs_dt_2m10","h_thetay_vs_dt_2m10", 40, -20, 20, 50, -3.15, 0);
}

void Test2DDeconv::analyze(art::Event const & e)
{

  constexpr double dyzcut_loose_2m = 1.0;
  constexpr double dyzcut_tight_2m = 0.1;
  constexpr double dyzcut_3m = 1.0;
  constexpr double dtimecut_3m = 10;

  art::InputTag trkTag("pandora");
  auto const& tracks   = proxy::getCollection<proxy::Tracks>(e,trkTag);
  h_ntracks->Fill(tracks.size());
  int ntracksall = 0;
  int ntracks50 = 0;
  for (const auto& t : tracks) {
    ntracksall++;
    std::cout << "track has npoints=" << t->NumberTrajectoryPoints() << " valid=" << t->CountValidPoints() << " and nhits=" << t.nHits() << std::endl;
    h_ntrackhits->Fill(t->CountValidPoints());
    if (t->CountValidPoints()<50) continue;
    ntracks50++;
    h_ntrackhits50->Fill(t->CountValidPoints());
    const auto& tkhits = t.hits();
    //
    // Here we first loop over hits in plane2, and then find those in plane1 and plane0
    //
    for (size_t ih = 0; ih < t->NumberTrajectoryPoints(); ih++) {
      if (t->HasValidPoint(ih)==0) continue;
      if (tkhits[ih]->View()!=2) continue;

      //std::cout << " view2 pos=" << t->Trajectory().LocationAtPoint(ih) << std::endl;

      size_t ih1 = t->NumberTrajectoryPoints();
      float dyz_min1 = 9999.; 
      for (size_t jh = std::max(size_t(0),ih-10); jh < std::min(t->NumberTrajectoryPoints(),ih+11); jh++) {
	if (t->HasValidPoint(jh)==0) continue;
	if (tkhits[jh]->View()!=1) continue;
	float dyz = sqrt( pow(t->Trajectory().LocationAtPoint(ih).Y()-t->Trajectory().LocationAtPoint(jh).Y(),2) +
			  pow(t->Trajectory().LocationAtPoint(ih).Z()-t->Trajectory().LocationAtPoint(jh).Z(),2) );
	//std::cout << " view1 pos=" << t->Trajectory().LocationAtPoint(jh) << " dyz=" << dyz << std::endl;
	//if (dyz>1.) continue;
	if (dyz<dyz_min1) {
	  dyz_min1 = dyz;
	  ih1 = jh;
	}
	//std::cout << " view1 pos=" << t->Trajectory().LocationAtPoint(jh) << std::endl;
      }
      
      size_t ih0 = t->NumberTrajectoryPoints();
      float dyz_min0 = 9999.; 
      for (size_t jh = std::max(size_t(0),(ih-10)); jh < std::min(t->NumberTrajectoryPoints(),ih+11); jh++) {
	if (t->HasValidPoint(jh)==0) continue;
	if (tkhits[jh]->View()!=0) continue;
	float dyz = sqrt( pow(t->Trajectory().LocationAtPoint(ih).Y()-t->Trajectory().LocationAtPoint(jh).Y(),2) +
			  pow(t->Trajectory().LocationAtPoint(ih).Z()-t->Trajectory().LocationAtPoint(jh).Z(),2) );
	//std::cout << " view0 pos=" << t->Trajectory().LocationAtPoint(jh) << " dyz=" << dyz << std::endl;
	//if (dyz>1.) continue;
	if (dyz<dyz_min0) {
	  dyz_min0 = dyz;
	  ih0 = jh;
	}
	//std::cout << " view0 pos=" << t->Trajectory().LocationAtPoint(jh) << std::endl;
      }

      // std::cout << "dyz_min1=" << dyz_min1 << " dyz_min0=" << dyz_min0 << std::endl;
      if (ih1<t->NumberTrajectoryPoints() && ih0<t->NumberTrajectoryPoints()) {
	// if (dyz_min1<1 && dyz_min0<1) std::cout << "BINGO" << std::endl; 
	// std::cout << " view2 pos=" << t->Trajectory().LocationAtPoint(ih) << std::endl;
	// std::cout << " view1 pos=" << t->Trajectory().LocationAtPoint(ih1) << std::endl;
	// std::cout << " view0 pos=" << t->Trajectory().LocationAtPoint(ih0) << std::endl;
	// std::cout << std::endl;
	const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
	std::cout << "offset0=" << detprop->GetXTicksOffset(tkhits[ih0]->WireID().Plane,tkhits[ih0]->WireID().TPC,tkhits[ih0]->WireID().Cryostat) << std::endl;
	std::cout << "offset1=" << detprop->GetXTicksOffset(tkhits[ih1]->WireID().Plane,tkhits[ih1]->WireID().TPC,tkhits[ih1]->WireID().Cryostat) << std::endl;
	std::cout << "offset2=" << detprop->GetXTicksOffset(tkhits[ih]->WireID().Plane,tkhits[ih]->WireID().TPC,tkhits[ih]->WireID().Cryostat) << std::endl;
	if (dyz_min1<dyzcut_3m && dyz_min0<dyzcut_3m) {
	  h_dt_3m10->Fill(tkhits[ih1]->PeakTime()-tkhits[ih0]->PeakTime());
	  h_dt_3m12->Fill(tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime());
	  h_dt_3m02->Fill(tkhits[ih0]->PeakTime()-tkhits[ih]->PeakTime());
	  if (std::abs(tkhits[ih1]->PeakTime()-tkhits[ih0]->PeakTime())<dtimecut_3m) {
	    h_dx_3m10->Fill(t->Trajectory().LocationAtPoint(ih1).X()-t->Trajectory().LocationAtPoint(ih0).X());
	    h_dx_3m12->Fill(t->Trajectory().LocationAtPoint(ih1).X()-t->Trajectory().LocationAtPoint(ih).X());
	    h_dx_3m02->Fill(t->Trajectory().LocationAtPoint(ih0).X()-t->Trajectory().LocationAtPoint(ih).X());
	  }
	}
      }
      //
      if (ih1<t->NumberTrajectoryPoints()) {
	h_dyz_vs_dt_2m12->Fill( tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime(), dyz_min1);
	if (dyz_min1>=dyzcut_loose_2m) continue;
	h_dt_loose_2m12->Fill(tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime());
	h_dx_loose_2m12->Fill( t->Trajectory().LocationAtPoint(ih1).X()-t->Trajectory().LocationAtPoint(ih).X() );
	if (dyz_min1<dyzcut_tight_2m) {
	  h_dt_tight_2m12->Fill(tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime());
	  h_dx_tight_2m12->Fill( t->Trajectory().LocationAtPoint(ih1).X()-t->Trajectory().LocationAtPoint(ih).X() );
	}
	// if (std::abs(tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime())>200) {
	//   std::cout << "ARGH" << std::endl;
	//   std::cout << "ih=" << ih << " time=" << tkhits[ih]->PeakTime() << " ih1=" << ih1 << " time=" << tkhits[ih1]->PeakTime() << std::endl;
	//   std::cout << " view2 pos=" << t->Trajectory().LocationAtPoint(ih) << std::endl;
	//   std::cout << " view1 pos=" << t->Trajectory().LocationAtPoint(ih1) << std::endl;
	//   for (const art::Ptr<recob::Hit>& h : t.hits()) {
	//     std::cout << "\t\thit wire=" << h->WireID() << " peak time=" << h->PeakTime() << std::endl;;
	//   }
	// }
      }
      //
      if (ih0<t->NumberTrajectoryPoints()) {
	h_dyz_vs_dt_2m02->Fill( tkhits[ih0]->PeakTime()-tkhits[ih]->PeakTime(), dyz_min0);
	if (dyz_min0>=dyzcut_loose_2m) continue;
	h_dt_loose_2m02->Fill(tkhits[ih0]->PeakTime()-tkhits[ih]->PeakTime());
	h_dx_loose_2m02->Fill( t->Trajectory().LocationAtPoint(ih0).X()-t->Trajectory().LocationAtPoint(ih).X() );
	if (dyz_min0<dyzcut_tight_2m) {
	  h_dt_tight_2m02->Fill(tkhits[ih0]->PeakTime()-tkhits[ih]->PeakTime());
	  h_dx_tight_2m02->Fill( t->Trajectory().LocationAtPoint(ih0).X()-t->Trajectory().LocationAtPoint(ih).X() );
	}
      }
    }
    //
    // Here we first loop over hits in plane0, and then find those in plane2 and plane1
    //
    for (size_t ih = 0; ih < t->NumberTrajectoryPoints(); ih++) {
      if (t->HasValidPoint(ih)==0) continue;
      if (tkhits[ih]->View()!=0) continue;

      //std::cout << " view2 pos=" << t->Trajectory().LocationAtPoint(ih) << std::endl;

      size_t ih1 = t->NumberTrajectoryPoints();
      float dyz_min1 = 9999.; 
      for (size_t jh = std::max(size_t(0),ih-10); jh < std::min(t->NumberTrajectoryPoints(),ih+11); jh++) {
	if (t->HasValidPoint(jh)==0) continue;
	if (tkhits[jh]->View()!=1) continue;
	float dyz = sqrt( pow(t->Trajectory().LocationAtPoint(ih).Y()-t->Trajectory().LocationAtPoint(jh).Y(),2) +
			  pow(t->Trajectory().LocationAtPoint(ih).Z()-t->Trajectory().LocationAtPoint(jh).Z(),2) );
	//std::cout << " view1 pos=" << t->Trajectory().LocationAtPoint(jh) << " dyz=" << dyz << std::endl;
	//if (dyz>1.) continue;
	if (dyz<dyz_min1) {
	  dyz_min1 = dyz;
	  ih1 = jh;
	}
	//std::cout << " view1 pos=" << t->Trajectory().LocationAtPoint(jh) << std::endl;
      }
      //
      size_t ih2 = t->NumberTrajectoryPoints();
      float dyz_min2 = 9999.; 
      for (size_t jh = std::max(size_t(0),ih-10); jh < std::min(t->NumberTrajectoryPoints(),ih+11); jh++) {
	if (t->HasValidPoint(jh)==0) continue;
	if (tkhits[jh]->View()!=2) continue;
	float dyz = sqrt( pow(t->Trajectory().LocationAtPoint(ih).Y()-t->Trajectory().LocationAtPoint(jh).Y(),2) +
			  pow(t->Trajectory().LocationAtPoint(ih).Z()-t->Trajectory().LocationAtPoint(jh).Z(),2) );
	//std::cout << " view1 pos=" << t->Trajectory().LocationAtPoint(jh) << " dyz=" << dyz << std::endl;
	//if (dyz>1.) continue;
	if (dyz<dyz_min2) {
	  dyz_min2 = dyz;
	  ih2 = jh;
	}
	//std::cout << " view1 pos=" << t->Trajectory().LocationAtPoint(jh) << std::endl;
      }
      //
      if (ih1<t->NumberTrajectoryPoints()) {
	h_dyz_vs_dt_2m10->Fill( tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime(), dyz_min1);
	if (dyz_min1>=dyzcut_loose_2m) continue;
	h_dt_loose_2m10->Fill(tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime());
	h_dx_loose_2m10->Fill( t->Trajectory().LocationAtPoint(ih1).X()-t->Trajectory().LocationAtPoint(ih).X() );
	h_wire0_vs_dt_2m10->Fill( tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime(), tkhits[ih]->WireID().Wire);
	h_wire1_vs_dt_2m10->Fill( tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime(), tkhits[ih1]->WireID().Wire);
	h_mult0_vs_dt_2m10->Fill( tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime(), tkhits[ih]->Multiplicity());
	h_mult1_vs_dt_2m10->Fill( tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime(), tkhits[ih1]->Multiplicity());
	h_theta_vs_dt_2m10->Fill( tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime(), t->Theta());
	h_phi_vs_dt_2m10->Fill( tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime(), t->Phi());
	h_thetax_vs_dt_2m10->Fill( tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime(), atan2(t->StartMomentumVector().X(), t->StartMomentumVector().Z()));
	h_thetay_vs_dt_2m10->Fill( tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime(), atan2(t->StartMomentumVector().Y(), t->StartMomentumVector().Z()));
	//
	if ( (tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime())>0 && (tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime())<1 &&
	     atan2(t->StartMomentumVector().Y(), t->StartMomentumVector().Z())>-1.4 && atan2(t->StartMomentumVector().Y(), t->StartMomentumVector().Z())<-1.2 ) {
	  //
	  std::cout << Form("event=(%i, %i, %i) track#%i nhits=%i phi=%f theta=%f plane0 wire=%i hit time=%f plane1 wire=%i hit time=%f",
			    e.run(),e.subRun(),e.event(),ntracksall,int(t.nHits()),t->Phi(),t->Theta(),tkhits[ih]->WireID().Wire,tkhits[ih]->PeakTime(),tkhits[ih1]->WireID().Wire,tkhits[ih1]->PeakTime())
		    << std::endl;
	  //
	}
	//
	if (dyz_min1<dyzcut_tight_2m) {
	  h_dt_tight_2m10->Fill(tkhits[ih1]->PeakTime()-tkhits[ih]->PeakTime());
	  h_dx_tight_2m10->Fill( t->Trajectory().LocationAtPoint(ih1).X()-t->Trajectory().LocationAtPoint(ih).X() );
	}
      }
      //
      if (ih2<t->NumberTrajectoryPoints()) {
	h_dyz_vs_dt_2m20->Fill( tkhits[ih]->PeakTime()-tkhits[ih2]->PeakTime(), dyz_min2);
	//if (dyz_min2>=dyzcut_loose_2m) continue;
      }
    }
    //
    for (size_t ih = 0; ih < t->NumberTrajectoryPoints(); ih++) {
      if (t->HasValidPoint(ih)==0) continue;
      if (tkhits[ih]->View()!=1) continue;
      //
      size_t ih0 = t->NumberTrajectoryPoints();
      float dyz_min0 = 9999.; 
      for (size_t jh = std::max(size_t(0),ih-10); jh < std::min(t->NumberTrajectoryPoints(),ih+11); jh++) {
	if (t->HasValidPoint(jh)==0) continue;
	if (tkhits[jh]->View()!=0) continue;
	float dyz = sqrt( pow(t->Trajectory().LocationAtPoint(ih).Y()-t->Trajectory().LocationAtPoint(jh).Y(),2) +
			  pow(t->Trajectory().LocationAtPoint(ih).Z()-t->Trajectory().LocationAtPoint(jh).Z(),2) );
	//std::cout << " view1 pos=" << t->Trajectory().LocationAtPoint(jh) << " dyz=" << dyz << std::endl;
	//if (dyz>1.) continue;
	if (dyz<dyz_min0) {
	  dyz_min0 = dyz;
	  ih0 = jh;
	}
	//std::cout << " view1 pos=" << t->Trajectory().LocationAtPoint(jh) << std::endl;
      }
      //
      if (ih0<t->NumberTrajectoryPoints()) {
	h_dyz_vs_dt_2m01->Fill( tkhits[ih]->PeakTime()-tkhits[ih0]->PeakTime(), dyz_min0);
	//if (dyz_min0>=dyzcut_loose_2m) continue;
      }
      //     
    }
  }
  h_ntracks50->Fill(ntracks50);
}

DEFINE_ART_MODULE(Test2DDeconv)
