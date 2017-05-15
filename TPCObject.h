/**
 * \class ubana::TPCObject
 *
 * \ingroup UBXSec
 *
 * \brief Data product to store a flash matching results
 * 
 *
 * \author $Author: Marco Del Tutto<marco.deltutto@physics.ox.ac.uk> $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2017/03/02 $
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Friday, May 15, 2017 at 16:34:34
 *
 */

#ifndef TPCObject_h
#define TPCObject_h

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include <vector>

namespace ubana {
  enum TPCObjectOrigin{
    kUnknown = -1,          // -1           
    kBeamNeutrino = 0,      // 0
    kCosmicRay,             // 1
    kMixed,                 // 2
  };
}

namespace ubana {

  class TPCObject {

  public:

    TPCObject();
    virtual ~TPCObject();

    // Setter methods
    void SetTracks(std::vector<recob::Track>);
    void SetVertex(recob::Vertex);
    void SetOrigin(ubana::TPCObjectOrigin);

    // Getter methods
    const std::vector<recob::Track>  & GetTracks()  const;
    const recob::Vertex              & GetVertex()  const;
    const ubana::TPCObjectOrigin     & GetOrigin()  const;

  private:

    std::vector<recob::Track>  fTracks;
    recob::Vertex              fVertex;
    ubana::TPCObjectOrigin     fOrigin;

 };
}

#endif /* TPCObject_h */
