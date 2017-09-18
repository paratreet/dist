#include <sstream>

#include "Main.h"
#include "mdt.h"
using namespace MultiphaseDistree;

Main::Main(CkArgMsg *m) {
    int nPieces=3*CkNumPes();
    CkArrayOptions opts(nPieces);
    
    CkPrintf("Main handing off work to %d tree pieces\n",nPieces);
    CProxy_TreePiece treePieceProxy = CProxy_TreePiece::ckNew(opts);
    
    CkCallback cb(CkReductionTarget(Main,done),thisProxy);
    treePieceProxy.runTreeWork(cb);
}

void Main::done() {
    CkPrintf("Main is now done running...\n");
    CkExit();
}


class TreePiece : public CBase_TreePiece 
{
public:
    TreePiece() {
        
    }
    TreePiece(CkMigrateMessage *m) {}
    
    /* entry */
    void runTreeWork(CkCallback &cb) {
        CkPrintf("Hello from index %d\n",thisIndex);
        contribute(cb);
    }
};



namespace MultiphaseDistree {
    void userInitializeReducers(void) {
    } 
};


#include "pgm.def.h"


