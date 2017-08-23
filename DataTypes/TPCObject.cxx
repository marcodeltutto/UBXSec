#include "TPCObject.h"
#include <vector>

namespace ubana {

  TPCObject::TPCObject() {
    fOrigin      = ubana::kUnknown;
  }

  TPCObject::~TPCObject(){
  }

  // Setter methoths
  void TPCObject::SetTracks (std::vector<recob::Track> tracks)    { this->fTracks = tracks;    }
  void TPCObject::SetPFPs   (std::vector<recob::PFParticle> pfps) { this->fPFParticles = pfps; }
  void TPCObject::SetVertex (recob::Vertex vertex)                { this->fVertex = vertex;    }       
  void TPCObject::SetOrigin (ubana::TPCObjectOrigin origin )      { this->fOrigin = origin;    }
  void TPCObject::SetMultiplicity(int p, int t, int s)            { this->fPfpMult = p; this->fTrackMult = t; this->fShowerMult = s;}

  // Getter methods
  const std::vector<recob::Track>       & TPCObject::GetTracks() const { return this->fTracks;      }
  const std::vector<recob::PFParticle>  & TPCObject::GetPFPs()   const { return this->fPFParticles; }
  const recob::Vertex                   & TPCObject::GetVertex() const { return this->fVertex;      }
  const ubana::TPCObjectOrigin          & TPCObject::GetOrigin() const { return this->fOrigin;      }

  const size_t TPCObject::GetNTracks()  const { return (this->fTracks).size();      }
  const size_t TPCObject::GetNPFP()     const { return (this->fPFParticles).size(); }
  const size_t TPCObject::GetNShowers() const { return (this->fShowers).size();     }

  const void TPCObject::GetMultiplicity(int &p, int &t, int &s) const { p = this->fPfpMult; t = this->fTrackMult; s = this->fShowerMult;}

}


