#include "Decomposer.h"
#include "ParamStorage.h"
#include "DoubleBufferedVec.h"
#include "Utilities.h"

extern ParamStorage parameters;
extern CProxy_TreePiece treePieceProxy;
extern CProxy_SplitterGroup splitterGroupProxy;

Decomposer::Decomposer(){
}

void Decomposer::decompose(BoundingBox &box, const CkCallback &cb){
  DoubleBufferedVec<Key> keys;
  keys.add(Key(1));
  keys.add(~Key(0));
  keys.sync();

  CkReductionMsg *msg;
  treePieceProxy.count(keys.get(), CkCallbackResumeThread((void *&)msg));

  for(int iteration = 0; ; iteration++){
    int *counts = (int *) msg->getData();
    int nCounts = msg->getSize()/sizeof(int);
    
    processCounts(counts, nCounts, keys); 
    keys.sync();

    delete msg;

    if(keys.size() == 0) break;

    treePieceProxy.count(keys.get(), CkCallbackResumeThread((void *&) msg));
  }
  
  int nPieces = splitters_.size();
  CkPrintf("[decomposer] nPieces %d\n", nPieces);

  if(nPieces > parameters.nPieces){
    CkPrintf("[decomposer] too few pieces; try with --nPieces %d\n", nPieces);
    CkAbort("too few pieces\n");
  }

  splitters_.quickSort();
  splitterGroupProxy.splitters(splitters_, CkCallbackResumeThread());
  treePieceProxy.flush(CkCallbackResumeThread());
  splitters_.resize(0);

  cb.send();
}

void Decomposer::processCounts(int *counts, int nCounts, DoubleBufferedVec<Key> &keys){
  CkAssert(2 * nCounts == keys.size());

  Real threshold = (DECOMP_TOLERANCE * Real(parameters.maxppc));
  for(int i = 0; i < nCounts; i++){
    Key from = keys.get(2*i);
    Key to = keys.get(2*i+1);

    int np = counts[i];
    if((Real) np > threshold){
      keys.add(from << 1);
      keys.add((from << 1)+1);

      keys.add((from << 1)+1);
      if(to == (~Key(0))){
        keys.add(~Key(0));
      }
      else{
        keys.add(to << 1);
      }
    }
    else{
      splitters_.push_back(Splitter(Utilities::toSplitter(from), Utilities::toSplitter(to), from, np));
    }
  }
}

void Decomposer::printKeys(DoubleBufferedVec<Key> &keys){
  ostringstream oss;
  for(int i = 0; i < keys.size(); i += 2){
    oss << "[" << keys.get(i) << ",";
    oss << keys.get(i+1) << ")";
    oss << ", ";
  }

  CkPrintf("decomposer keys: %s\n", oss.str().c_str());
}


