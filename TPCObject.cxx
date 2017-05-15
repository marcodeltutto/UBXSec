#include "TPCObject.h"
#include <vector>

namespace ubana {

  TPCObject::TPCObject() {
    fOrigin      = ubana::kUnknown;
  }

  TPCObject::~TPCObject(){
  }

  // Setter methoths
  void TPCObject::SetTracks (std::vector<recob::Track> tracks)  { this->fTracks = tracks; }
  void TPCObject::SetVertex (recob::Vertex vertex)              { this->fVertex = vertex; }       
  void TPCObject::SetOrigin (ubana::TPCObjectOrigin origin )    { this->fOrigin = origin; }

  // Getter methods
  const std::vector<recob::Track>  & TPCObject::GetTracks() const { return this->fTracks; }
  const recob::Vertex              & TPCObject::GetVertex() const { return this->fVertex; }
  const ubana::TPCObjectOrigin     & TPCObject::GetOrigin() const { return this->fOrigin; }

}


