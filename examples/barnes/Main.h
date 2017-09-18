#ifndef MAIN_H
#define MAIN_H

#include "barnes.decl.h"
#include "TemplateInstanceTypedefs.h"

class Main : public CBase_Main {
  CkCallback traversalCallback_;
  int traversalsDone_;
  int numConcurrentTraversals_;

  private:
  void setTraversalCallback(const CkCallback &cb);
  void checkTraversalsDone();


  public:
  Main(CkArgMsg *msg);
  void commence();

  void gravityDone(CkReductionMsg *m);
  void sphDone(CkReductionMsg *m);
};

#endif // MAIN_H
