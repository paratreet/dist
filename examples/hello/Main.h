#ifndef MAIN_H
#define MAIN_H

#include "pgm.decl.h"

class Main : public CBase_Main {
  public:
  Main(CkArgMsg *msg);
  void done(void);
};

#endif // MAIN_H
