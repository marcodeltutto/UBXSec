/**
 * \class ubana::MCGhost
 *
 * \ingroup UBXSec
 *
 * \brief Data product to store a TPC Object
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

#ifndef MCGHOST_H
#define MCGHOST_H

#include <string>

namespace ubana {

  class MCGhost {

  public:

    MCGhost();
    virtual ~MCGhost();

    // Setter methods
    void SetMode(std::string);

    // Getter methods
    const std::string & GetMode() const;

  private:

    std::string _mode;
 };
}

#endif /* MCGHOST_H */
