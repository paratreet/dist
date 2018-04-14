#ifndef LB_ON_OFF_GROUP_H
#define LB_ON_OFF_GROUP_H

#include "barnes.decl.h"

class LbOnOffGroup : public CBase_LbOnOffGroup {
  public:
  LbOnOffGroup();

  void on(const CkCallback &cb);
  void off(const CkCallback &cb);
};

#endif // SPLITTER_GROUP_H
