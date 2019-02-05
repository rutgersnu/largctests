////////////////////////////////////////////////////////////////////////
// Class:       TestHitSorting
// Plugin Type: analyzer (art v3_00_00)
// File:        TestHitSorting_module.cc
//
// Generated at Thu Dec 27 14:18:37 2018 by Giuseppe Cerati using cetskelgen
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

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/ForEachAssociatedGroup.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"

#include "lardata/RecoBaseProxy/Track.h"

class TestHitSorting;


class TestHitSorting : public art::EDAnalyzer {
public:
  explicit TestHitSorting(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TestHitSorting(TestHitSorting const&) = delete;
  TestHitSorting(TestHitSorting&&) = delete;
  TestHitSorting& operator=(TestHitSorting const&) = delete;
  TestHitSorting& operator=(TestHitSorting&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

};


TestHitSorting::TestHitSorting(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void TestHitSorting::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  using namespace std;

  string label = "pandoraTrack";
  
  const auto& inputTracks = e.getValidHandle<vector<recob::Track> >(label);

  art::Assns<recob::Track, recob::Hit> tkHitsAssn = *e.getValidHandle<art::Assns<recob::Track, recob::Hit> >(label);
  const auto& trackHitsGroups = util::associated_groups(tkHitsAssn);

  art::FindMany<recob::Hit> assocHits(inputTracks, e, label);
  
  for (size_t itk = 0; itk<inputTracks->size();++itk) {
    art::Ptr<recob::Track> track(inputTracks,itk);

    cout << "track points=" << track->NumberTrajectoryPoints() << " valid=" << track->CountValidPoints() << " length=" << track->Length() << std::endl;
    
    std::vector<art::Ptr<recob::Hit> > inHits;
    decltype(auto) hitsRange = util::groupByIndex(trackHitsGroups, track.key());
    for (art::Ptr<recob::Hit> const& hit: hitsRange) inHits.push_back(hit);
    std::cout << "groupByIndex n hits=" << inHits.size() << std::endl;
    for (auto h : inHits) {
      cout << "p=" << h->WireID().Plane << " w=" << h->WireID().Wire << " t=" << h->PeakTime() << endl;
    }

    const std::vector<recob::Hit const*>& hits = assocHits.at(track.key());
    std::cout << "FindMany n hits=" << hits.size() << std::endl;
    for (auto h : hits) {
      cout << "p=" << h->WireID().Plane << " w=" << h->WireID().Wire << " t=" << h->PeakTime() << endl;
    }
  }

  auto const& tracks = proxy::getCollection<proxy::Tracks>(e,label);
  for (const auto& t : tracks) {
    std::cout << "proxy n hits=" << t.nHits() << std::endl;
    for (const art::Ptr<recob::Hit>& h : t.hits()) {
      cout << "p=" << h->WireID().Plane << " w=" << h->WireID().Wire << " t=" << h->PeakTime() << endl;
    }    
  }
}

DEFINE_ART_MODULE(TestHitSorting)
