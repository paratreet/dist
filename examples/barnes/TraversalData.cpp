#include "TraversalData.h"

std::ostream &operator<<(std::ostream &out, const TraversalData &data){
  out << "[open: " << data.open();
  out << " pn: " << data.pn();
  out << " pp: " << data.pp();
  out << "]";
  return out;
}
