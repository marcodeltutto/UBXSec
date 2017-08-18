#include "MCGhost.h"

#include <string>

namespace ubana {

  MCGhost::MCGhost() {
    _mode = "unknown";
  }

  MCGhost::~MCGhost(){
  }

  // Setter methoths
  void MCGhost::SetMode (std::string mode) { this->_mode = mode;    }

  // Getter methods
  const std::string & MCGhost::GetMode() const { return this->_mode;      }
}


