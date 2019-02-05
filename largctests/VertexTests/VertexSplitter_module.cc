////////////////////////////////////////////////////////////////////////
// Class:       VertexSplitter
// Plugin Type: producer (art v2_08_04)
// File:        VertexSplitter_module.cc
//
// Generated at Fri Dec 29 11:09:13 2017 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_01_01.
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

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "larreco/RecoAlg/Geometric3DVertexFitter.h"

#include <memory>

class VertexSplitter;

class VertexSplitter : public art::EDProducer {
public:

  struct Inputs {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<art::InputTag> inputPFParticleLabel {
       Name("inputPFParticleLabel"),
       Comment("Label of recob::PFParticle Collection to be fit")
    };
    fhicl::Atom<art::InputTag> inputVertexLabel {
       Name("inputVertexLabel"),
       Comment("Label of recob::Vertex Collection associated to PFParticles")
    };
  };

  struct Config {
    using Name = fhicl::Name;
    fhicl::Table<VertexSplitter::Inputs> inputs {
      Name("inputs"),
    };
    fhicl::Table<trkf::Geometric3DVertexFitter::Config> geom3dvtxfit {
      Name("geom3dvtxfit")
    };
    fhicl::Table<trkf::TrackStatePropagator::Config> propagator {
      Name("propagator")
    };
  };
  using Parameters = art::EDProducer::Table<Config>;

  explicit VertexSplitter(Parameters const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  VertexSplitter(VertexSplitter const &) = delete;
  VertexSplitter(VertexSplitter &&) = delete;
  VertexSplitter & operator = (VertexSplitter const &) = delete;
  VertexSplitter & operator = (VertexSplitter &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:
  art::InputTag pfParticleInputTag;
  art::InputTag vertexInputTag;
  trkf::Geometric3DVertexFitter fitter;

};


VertexSplitter::VertexSplitter(Parameters const & p)
  : pfParticleInputTag(p().inputs().inputPFParticleLabel())
  , vertexInputTag(p().inputs().inputVertexLabel())
  , fitter(p().geom3dvtxfit,p().propagator)
{
  produces<std::vector<recob::Vertex> >("even");
  produces<art::Assns<recob::PFParticle, recob::Vertex> >("even");
  produces<art::Assns<recob::Vertex, recob::Track, recob::VertexAssnMeta> >("even");
  produces<std::vector<recob::Vertex> >("odd");
  produces<art::Assns<recob::PFParticle, recob::Vertex> >("odd");
  produces<art::Assns<recob::Vertex, recob::Track, recob::VertexAssnMeta> >("odd");
}

void VertexSplitter::produce(art::Event & e)
{
  // Implementation of required member function here.

  using namespace std;
  using namespace trkf;

  auto outputVertices1   = make_unique<vector<recob::Vertex> >();
  auto outputPFVxAssn1   = make_unique<art::Assns<recob::PFParticle, recob::Vertex> >();
  auto outputVxTkMtAssn1 = make_unique<art::Assns<recob::Vertex, recob::Track, recob::VertexAssnMeta> >();
  //
  auto outputVertices2   = make_unique<vector<recob::Vertex> >();
  auto outputPFVxAssn2   = make_unique<art::Assns<recob::PFParticle, recob::Vertex> >();
  auto outputVxTkMtAssn2 = make_unique<art::Assns<recob::Vertex, recob::Track, recob::VertexAssnMeta> >();
  //
  const auto& inputPFParticle = e.getValidHandle<vector<recob::PFParticle> >(pfParticleInputTag);
  auto assocVertex = unique_ptr<art::FindManyP<recob::Vertex> >(new art::FindManyP<recob::Vertex>(inputPFParticle, e, vertexInputTag));
  const auto& inputVertexColl = e.getValidHandle<vector<recob::Vertex> >(vertexInputTag);
  auto assocTracks = unique_ptr<art::FindManyP<recob::Track, recob::VertexAssnMeta> >(new art::FindManyP<recob::Track, recob::VertexAssnMeta>(inputVertexColl, e, vertexInputTag));

  art::PtrMaker<recob::Vertex> vtxPtrMaker1(e, "even");
  art::PtrMaker<recob::Vertex> vtxPtrMaker2(e, "odd");

  for (size_t iPF = 0; iPF < inputPFParticle->size(); ++iPF) {
    //
    art::Ptr<recob::PFParticle> pfp(inputPFParticle, iPF);
    if (pfp->IsPrimary()==false || pfp->NumDaughters()<2) continue;

    auto vertices = assocVertex->at(pfp.key());
    if (vertices.size()!=1) continue;
    auto vertex = vertices[0];
    auto tracks = assocTracks->at(vertex.key());
    auto metas  = assocTracks->data(vertex.key());
    //
    vector< art::Ptr<recob::Track> > usedTracks;
    for (size_t it=0; it<tracks.size(); ++it) {
      if (metas[it]->status()==recob::VertexAssnMeta::IncludedInFit) usedTracks.push_back(tracks[it]);
    }
    if (usedTracks.size()!=4) continue;
    //
    std::cout << "Found vertex with ntracks=4!!!!!" << std::endl;
    std::cout << "original pos=" << vertex->position() << std::endl;
    //
    vector< art::Ptr<recob::Track> > tracks1;
    tracks1.push_back(usedTracks[0]);
    tracks1.push_back(usedTracks[2]);
    VertexWrapper vtx1 = fitter.fitTracks(tracks1);
    if (vtx1.isValid()==false) continue;
    vtx1.setVertexId(vertex->ID());
    std::cout << "split1 pos=" << vtx1.position() << std::endl;
    auto meta1 = fitter.computeMeta(vtx1, tracks1);
    outputVertices1->emplace_back(vtx1.vertex());
    const art::Ptr<recob::Vertex> aptr1 = vtxPtrMaker1(outputVertices1->size()-1);
    outputPFVxAssn1->addSingle( art::Ptr<recob::PFParticle>(inputPFParticle, iPF), aptr1);
    size_t itt1 = 0;
    for (auto t1 : tracks1) {
      outputVxTkMtAssn1->addSingle(aptr1, t1, meta1[itt1]);
      itt1++;
    }
    //
    vector< art::Ptr<recob::Track> > tracks2;
    tracks2.push_back(usedTracks[1]);
    tracks2.push_back(usedTracks[3]);
    VertexWrapper vtx2 = fitter.fitTracks(tracks2);
    if (vtx2.isValid()==false) continue;
    std::cout << "split2 pos=" << vtx2.position() << std::endl;
    vtx2.setVertexId(vertex->ID());
    auto meta2 = fitter.computeMeta(vtx2, tracks2);
    outputVertices2->emplace_back(vtx2.vertex());
    const art::Ptr<recob::Vertex> aptr2 = vtxPtrMaker2(outputVertices2->size()-1);
    outputPFVxAssn2->addSingle( art::Ptr<recob::PFParticle>(inputPFParticle, iPF), aptr2);
    size_t itt2 = 0;
    for (auto t2 : tracks2) {
      outputVxTkMtAssn2->addSingle(aptr2, t2, meta2[itt2]);
      itt2++;
    }
  }
  //
  e.put(std::move(outputVertices1),"even");
  e.put(std::move(outputPFVxAssn1),"even");
  e.put(std::move(outputVxTkMtAssn1),"even");
  e.put(std::move(outputVertices2),"odd");
  e.put(std::move(outputPFVxAssn2),"odd");
  e.put(std::move(outputVxTkMtAssn2),"odd");
}

DEFINE_ART_MODULE(VertexSplitter)
