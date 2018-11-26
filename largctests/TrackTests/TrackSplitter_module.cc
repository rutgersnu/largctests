////////////////////////////////////////////////////////////////////////
// Class:       TrackSplitter
// Plugin Type: producer (art v2_06_03)
// File:        TrackSplitter_module.cc
//
// Generated at Thu May 25 08:40:39 2017 by Giuseppe Cerati using cetskelgen
// from cetlib version v2_03_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "larreco/ClusterFinder/ClusterCreator.h"
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"

#include <memory>

class TrackSplitter;


class TrackSplitter : public art::EDProducer {
public:

  explicit TrackSplitter(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  TrackSplitter(TrackSplitter const &) = delete;
  TrackSplitter(TrackSplitter &&) = delete;
  TrackSplitter & operator = (TrackSplitter const &) = delete;
  TrackSplitter & operator = (TrackSplitter &&) = delete;

  recob::Cluster BuildCluster(const int id, const art::PtrVector<recob::Hit> &hitVector, cluster::ClusterParamsAlgBase &algo);

  // Required functions.
  void produce(art::Event & e) override;

private:
  art::InputTag pfParticleInputTag;
  bool reverseHits;
  int hitSubsetMode;
  int midPointHalfGap;
  int interleaveStep;
};


TrackSplitter::TrackSplitter(fhicl::ParameterSet const & p)
  : pfParticleInputTag(p.get<std::string>("pfParticleInputTag")),
    reverseHits(p.get<bool>("reverseHits")),
    hitSubsetMode(p.get<int>("hitSubsetMode")),
    midPointHalfGap(p.get<int>("midPointHalfGap")),
    interleaveStep(p.get<int>("interleaveStep"))
{
  produces<std::vector<recob::Track> >();
  produces<std::vector<recob::Cluster> >();
  produces<std::vector<recob::PFParticle> >();
  produces<art::Assns<recob::Track, recob::Hit> >();
  produces<art::Assns<recob::PFParticle, recob::Track> >();
  produces<art::Assns<recob::PFParticle, recob::Cluster> >();
  produces<art::Assns<recob::PFParticle, recob::Vertex> >();
  produces<art::Assns<recob::Cluster, recob::Hit> >();
}

void TrackSplitter::produce(art::Event & e)
{
  // Implementation of required member function here.
  auto outputClusters = std::make_unique<std::vector<recob::Cluster> >();
  auto outputTracks   = std::make_unique<std::vector<recob::Track> >();
  auto outputPFs      = std::make_unique<std::vector<recob::PFParticle> >();
  auto outputHitsAssn = std::make_unique<art::Assns<recob::Track, recob::Hit> >();
  auto outputPFTkAssn = std::make_unique<art::Assns<recob::PFParticle, recob::Track> >();
  auto outputPFClAssn = std::make_unique<art::Assns<recob::PFParticle, recob::Cluster> >();
  auto outputPFVxAssn = std::make_unique<art::Assns<recob::PFParticle, recob::Vertex> >();
  auto outputClHtAssn = std::make_unique<art::Assns<recob::Cluster, recob::Hit> >();
  //
  art::PtrMaker<recob::Cluster> makeClusterPtr(e, *this);
  art::PtrMaker<recob::Track> makeTrackPtr(e, *this);
  art::PtrMaker<recob::PFParticle> makePFParticlePtr(e, *this);
  //
  art::ValidHandle<std::vector<recob::PFParticle> > inputPFParticle = e.getValidHandle<std::vector<recob::PFParticle> >(pfParticleInputTag);
  std::unique_ptr<art::FindManyP<recob::Track> >  assocTracks = std::unique_ptr<art::FindManyP<recob::Track> >(new art::FindManyP<recob::Track>(inputPFParticle, e, pfParticleInputTag));
  std::unique_ptr<art::FindManyP<recob::Vertex> > assocVertices = std::unique_ptr<art::FindManyP<recob::Vertex> >(new art::FindManyP<recob::Vertex>(inputPFParticle, e, pfParticleInputTag));
  auto const& tkHitsAssn = *e.getValidHandle<art::Assns<recob::Track, recob::Hit> >(pfParticleInputTag);
  //
  for (unsigned int iPF = 0; iPF < inputPFParticle->size(); ++iPF) {
    //
    const recob::PFParticle& pf = inputPFParticle->at(iPF);
    recob::PFParticle outPF(pf);//copy PFParticle hierarchy, but need to remake tracks and clusters
    //
    outputPFs->emplace_back(std::move(outPF));
    art::Ptr<recob::PFParticle> aptrpf = makePFParticlePtr(outputPFs->size()-1);
    //
    const std::vector<art::Ptr<recob::Track> >& tracks = assocTracks->at(iPF);
    //
    for (unsigned int iTrack = 0; iTrack < tracks.size(); ++iTrack) {
      //
      const recob::Track& track = *tracks[iTrack];
      art::Ptr<recob::Track> ptrack = tracks[iTrack];
      //
      // prepare the hits - this is not computationally optimal, but at least preserves the order unlike FindManyP
      std::vector<art::Ptr<recob::Hit> > inHits;
      for (auto it = tkHitsAssn.begin(); it!=tkHitsAssn.end(); ++it) {
	if (it->first == ptrack) inHits.push_back(it->second);
	else if (inHits.size()>0) break;
      }
      //
      // here we do the splitting
      //
      recob::tracking::Positions_t    positions;
      recob::tracking::Momenta_t      momenta;
      recob::TrackTrajectory::Flags_t flags;
      art::PtrVector<recob::Hit> outHits;
      const int beg = (reverseHits ? inHits.size()-1 : 0);
      const int end = (reverseHits ? -1 : inHits.size());
      int tknpoints = track.NPoints();
      for (int ihit = beg; ihit!=end; (reverseHits ? ihit-- : ihit++)) {
	auto hit = inHits[ihit];
	//
	if (hitSubsetMode==1 && hit->WireID().Wire % 2 == 1) continue;
	if (hitSubsetMode==2 && hit->WireID().Wire % 2 == 0) continue;
	//
	if (hitSubsetMode==3 && ihit>(inHits.size()*0.5-midPointHalfGap)) continue;
	if (hitSubsetMode==4 && ihit<(inHits.size()*0.5+midPointHalfGap)) continue;
	//
	int tc4r = (ihit % (4*interleaveStep))/interleaveStep;
	assert(tc4r<4);
	if (hitSubsetMode==5 && tc4r!=0) continue;
	if (hitSubsetMode==6 && tc4r!=2) continue;
	//
	int tc6r = (ihit % (6*interleaveStep))/interleaveStep;
	assert(tc6r<6);
	if (hitSubsetMode==7 && tc6r!=0 && tc6r!=1) continue;
	if (hitSubsetMode==8 && tc6r!=3 && tc6r!=4) continue;
	//
	//
	outHits.push_back(hit);
	if (ihit>=tknpoints) continue;
	positions.push_back(track.Trajectory().LocationAtPoint(ihit));
	momenta.push_back(track.Trajectory().MomentumVectorAtPoint(ihit));
	flags.push_back(track.Trajectory().FlagsAtPoint(ihit));
      }
      //
      // if track splitting fails make a dummy track with 2 hits, which will fails cuts later on
      // but at least does not break the PFP hierarchy
      //
      int nvalid = 0;
      for (auto f : flags) {
	if (f.isPointValid()) nvalid++;
      }
      if (nvalid<2) {
	positions.clear();
	momenta.clear();
	flags.clear();
	outHits.clear();
	int count = 0;
	for (int ihit = beg; ihit!=end; (reverseHits ? ihit-- : ihit++)) {
	  if (count>1) break;
	  auto hit = inHits[ihit];
	  outHits.push_back(hit);
	  positions.push_back(track.Trajectory().LocationAtPoint(ihit));
	  momenta.push_back(track.Trajectory().MomentumVectorAtPoint(ihit));
	  flags.push_back(recob::TrackTrajectory::PointFlags_t());
	  count++;
	}
      }
      //
      // make clusters for split track 
      //
      cluster::StandardClusterParamsAlg ClusterParamAlgo;
      for (unsigned int plane=0;plane<3;plane++) {
	//
	art::PtrVector<recob::Hit> viewHits;
	for (auto hit : outHits) {
	  //fixme does not work for multi-TPC
	  if (hit->WireID().planeID().Plane!=plane) continue;
	  viewHits.push_back(hit);
	}
	if (viewHits.empty()) continue;
	recob::Cluster outCluster = BuildCluster(outputClusters->size(), viewHits, ClusterParamAlgo);
	outputClusters->emplace_back(std::move(outCluster));
	art::Ptr<recob::Cluster> aptrcl = makeClusterPtr(outputClusters->size()-1);
	outputPFClAssn->addSingle(aptrpf, aptrcl);
	for (auto phit : viewHits)  outputClHtAssn->addSingle(aptrcl, phit);
      }
      //
      // make split track 
      //
      recob::Track outTrack = recob::Track(recob::TrackTrajectory(std::move(positions),std::move(momenta),std::move(flags),track.HasMomentum()),
					   0,-1.,0,recob::tracking::SMatrixSym55(),recob::tracking::SMatrixSym55(),track.ID());
      //
      // std::cout << "orgnl track start=" << track.Trajectory().Start() << " dir" <<  track.Trajectory().StartDirection() << " nh=" << track.CountValidPoints() << " id=" << track.ID() << std::endl;
      // std::cout << "split track start=" << outTrack.Trajectory().Start() << " dir" <<  outTrack.Trajectory().StartDirection() << " nh=" << outTrack.CountValidPoints() << " id=" << outTrack.ID() << std::endl;
      // std::cout << std::endl;
      //
      outputTracks->emplace_back(std::move(outTrack));
      art::Ptr<recob::Track> aptrtk= makeTrackPtr(outputTracks->size()-1);
      for (auto const& trhit: outHits) {
	outputHitsAssn->addSingle(aptrtk, trhit);
      }
      outputPFTkAssn->addSingle(aptrpf, aptrtk);
    }
    //
    // vertices are not remade at this stage, but need the association
    //
    const std::vector<art::Ptr<recob::Vertex> >& vertices = assocVertices->at(iPF);
    for (unsigned int iVertex = 0; iVertex < vertices.size(); ++iVertex) {
      art::Ptr<recob::Vertex> pvertex = vertices[iVertex];
      outputPFVxAssn->addSingle(aptrpf, pvertex);
    }
    //
  }
  e.put(std::move(outputTracks));
  e.put(std::move(outputClusters));
  e.put(std::move(outputPFs));
  e.put(std::move(outputHitsAssn));
  e.put(std::move(outputPFTkAssn));
  e.put(std::move(outputPFClAssn));
  e.put(std::move(outputPFVxAssn));
  e.put(std::move(outputClHtAssn));
}

recob::Cluster TrackSplitter::BuildCluster(const int id, const art::PtrVector<recob::Hit> &hitVector, cluster::ClusterParamsAlgBase &algo) {
  //
  // Fill list of cluster properties
  geo::View_t view(geo::kUnknown);
  geo::PlaneID planeID;
  //
  double startWire(+std::numeric_limits<float>::max()), sigmaStartWire(0.0);
  double startTime(+std::numeric_limits<float>::max()), sigmaStartTime(0.0);
  double endWire(-std::numeric_limits<float>::max()), sigmaEndWire(0.0);
  double endTime(-std::numeric_limits<float>::max()), sigmaEndTime(0.0);
  //
  std::vector<recob::Hit const*> hits_for_params;
  hits_for_params.reserve(hitVector.size());
  //
  for (const art::Ptr<recob::Hit> &hit : hitVector) {
    //
    const double thisWire(hit->WireID().Wire);
    const double thisWireSigma(0.5);
    const double thisTime(hit->PeakTime());
    const double thisTimeSigma(double(2.*hit->RMS()));
    const geo::View_t thisView(hit->View());
    const geo::PlaneID thisPlaneID(hit->WireID().planeID());
    //
    if (geo::kUnknown == view) {
      view = thisView;
      planeID = thisPlaneID;
    }
    //
    hits_for_params.push_back(&*hit);
    //
    if (thisWire < startWire || (thisWire == startWire && thisTime < startTime)) {
      startWire = thisWire;
      sigmaStartWire = thisWireSigma;
      startTime = thisTime;
      sigmaStartTime = thisTimeSigma;
    }
    //
    if (thisWire > endWire || (thisWire == endWire && thisTime > endTime)) {
      endWire = thisWire;
      sigmaEndWire = thisWireSigma;
      endTime = thisTime;
      sigmaEndTime = thisTimeSigma;
    }
  }
  //
  // feed the algorithm with all the cluster hits
  algo.SetHits(hits_for_params);
  //
  // create the recob::Cluster directly in the vector
  return cluster::ClusterCreator(algo,startWire,sigmaStartWire,startTime,sigmaStartTime,
				 endWire,sigmaEndWire,endTime,sigmaEndTime,id,view,planeID,recob::Cluster::Sentry).move();
}

DEFINE_ART_MODULE(TrackSplitter)
