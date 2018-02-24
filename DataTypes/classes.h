
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"

#include "uboone/UBXSec/DataTypes/FlashMatch.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "uboone/UBXSec/DataTypes/MCGhost.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "uboone/UBXSec/DataTypes/UBXSecEvent.h"

#include "uboone/UBXSec/DataTypes/SelectionResult.h"


#include <vector>

template class art::Assns<anab::FlashMatch,recob::PFParticle>;
template class art::Assns<recob::PFParticle,anab::FlashMatch>;

template class art::Wrapper<art::Assns<anab::FlashMatch,  recob::PFParticle>>;
template class art::Wrapper<art::Assns<recob::PFParticle, anab::FlashMatch >>;


template class std::vector<ubana::FlashMatch>;

template class art::Assns<ubana::FlashMatch,recob::PFParticle>;
template class art::Assns<recob::PFParticle,ubana::FlashMatch>;
template class art::Assns<ubana::FlashMatch,recob::Track,void>;
template class art::Assns<recob::Track,ubana::FlashMatch,void>;

template class art::Wrapper<art::Assns<ubana::FlashMatch,  recob::PFParticle>>;
template class art::Wrapper<art::Assns<recob::PFParticle, ubana::FlashMatch >>;
template class art::Wrapper<art::Assns<ubana::FlashMatch,recob::Track,void> >;
template class art::Wrapper<art::Assns<recob::Track,ubana::FlashMatch,void> >;




template class std::vector<ubana::TPCObject>;

template class art::Assns<ubana::TPCObject,recob::PFParticle>;
template class art::Assns<recob::PFParticle,ubana::TPCObject>;
template class art::Assns<ubana::TPCObject,recob::Track,void>;
template class art::Assns<recob::Track,ubana::TPCObject,void>;
template class art::Assns<ubana::TPCObject,recob::Shower,void>;
template class art::Assns<recob::Shower,ubana::TPCObject,void>;
template class art::Assns<ubana::TPCObject,recob::Vertex,void>;
template class art::Assns<recob::Vertex,ubana::TPCObject,void>;
template class art::Assns<ubana::TPCObject,ubana::FlashMatch,void>;
template class art::Assns<ubana::FlashMatch,ubana::TPCObject,void>;
template class art::Assns<ubana::TPCObject,anab::CosmicTag,void>;
template class art::Assns<anab::CosmicTag,ubana::TPCObject,void>;

template class art::Wrapper<art::Assns<ubana::TPCObject,  recob::PFParticle>>;
template class art::Wrapper<art::Assns<recob::PFParticle, ubana::TPCObject >>;
template class art::Wrapper<art::Assns<recob::Track,ubana::TPCObject,void> >;
template class art::Wrapper<art::Assns<ubana::TPCObject,recob::Track,void> >;
template class art::Wrapper<art::Assns<recob::Shower,ubana::TPCObject,void> >;
template class art::Wrapper<art::Assns<ubana::TPCObject,recob::Shower,void> >;
template class art::Wrapper<art::Assns<recob::Vertex,ubana::TPCObject,void> >;
template class art::Wrapper<art::Assns<ubana::TPCObject,recob::Vertex,void> >;
template class art::Wrapper<art::Assns<ubana::TPCObject,ubana::FlashMatch,void> >;
template class art::Wrapper<art::Assns<ubana::FlashMatch,ubana::TPCObject,void> >;
template class art::Wrapper<art::Assns<ubana::TPCObject,anab::CosmicTag,void> >;
template class art::Wrapper<art::Assns<anab::CosmicTag,ubana::TPCObject,void> >;


template class std::vector<ubana::MCGhost>;

template class art::Assns<simb::MCParticle,ubana::MCGhost>;
template class art::Assns<ubana::MCGhost,simb::MCParticle>;
template class art::Assns<recob::PFParticle,ubana::MCGhost>;
template class art::Assns<ubana::MCGhost,recob::PFParticle>;
template class art::Wrapper<art::Assns<simb::MCParticle,ubana::MCGhost>>;
template class art::Wrapper<art::Assns<ubana::MCGhost,simb::MCParticle>>;
template class art::Wrapper<art::Assns<recob::PFParticle,ubana::MCGhost>>;
template class art::Wrapper<art::Assns<ubana::MCGhost,recob::PFParticle>>;


template class std::vector<ubana::SelectionResult>;
template class art::Wrapper<art::Assns<ubana::TPCObject,ubana::SelectionResult,void> >;
template class art::Wrapper<art::Assns<ubana::SelectionResult,ubana::TPCObject,void> >;

