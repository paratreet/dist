#include "LbOnOffGroup.h"

LbOnOffGroup::LbOnOffGroup(){
}

void LbOnOffGroup::on(const CkCallback &cb){
  LBTurnInstrumentOn();
  contribute(cb);
}

void LbOnOffGroup::off(const CkCallback &cb){
  LBTurnInstrumentOff();
  contribute(cb);
}
