#ifndef PARAMS_H
#define PARAMS_H

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "Environment.h"

using namespace std;

class Param {
  bool hasValue_;
  bool extractedValue_;
  bool compulsory_;

  public:
  Param(bool hasValue, bool compulsory) : 
    hasValue_(hasValue),
    extractedValue_(false),
    compulsory_(compulsory)
  {}


  virtual void extract(const char *val) = 0;
  virtual void push() = 0;

  bool hasValue() const {
    return hasValue_;
  }

  void hasExtractedValue(){
    extractedValue_ = true;
  }

  bool extractedValue() const {
    return extractedValue_;
  }

  bool isCompulsory() const {
    return compulsory_;
  }
};

template<typename ParamType>
class TypedParam : public Param {
  ParamType *saveTo_;
  ParamType parsedValue_;
  ParamType defaultValue_;

  public:
  TypedParam(ParamType *saveTo, ParamType defaultValue, bool hasValue, bool compulsory) : 
    Param(hasValue, compulsory), 
    saveTo_(saveTo),
    defaultValue_(defaultValue)
  {}

  TypedParam(ParamType *saveTo, bool hasValue, bool compulsory) : 
    Param(hasValue, compulsory), 
    saveTo_(saveTo)
  {}

  void extract(const char *val) {}

  void push(){
    if(extractedValue()) *saveTo_ = parsedValue_;
    else *saveTo_ = defaultValue_;
  }
};

#include <stdlib.h>
#include <assert.h>

template<>
inline void TypedParam<bool>::extract(const char *val){
  string valString(val);
  if(valString == "true"){
    parsedValue_ = true;  
  }
  else{
    parsedValue_ = false;
  }
}

template<>
inline void TypedParam<int>::extract(const char *val){
  parsedValue_ = atoi(val);
}

template<>
inline void TypedParam<unsigned int>::extract(const char *val){
  parsedValue_ = strtoul(val, NULL, 0);
}

template<>
inline void TypedParam<float>::extract(const char *val){
  parsedValue_ = atof(val);
}

template<>
inline void TypedParam<double>::extract(const char *val){
  parsedValue_ = strtod(val, NULL);
}

template<>
inline void TypedParam<string>::extract(const char *val){
  parsedValue_ = string(val);
}

class ParamCollection {
  map<string, Param *> params_;

  private:
  // these are 'visitor' classes for param list
  // this one extracts all values required by the
  // parallel driver of moose
  class ParamExtractionWorker {
    public:
    void matchValue(int argi, int argc, char **argv, Param *param){
      assert(argi + 1 < argc);
      param->extract(argv[argi + 1]);
    }

    void matchNoValue(int argi, int argc, char **argv, Param *param){
      param->extract("true");
    }

    void noMatch(int argi, int argc, char **argv){
      cerr << "unrecognized parameter `" << argv[argi] << "'; skipping" << endl;
    }
  };
  
  // this one "purges" the arg list for moose-core,
  // i.e. it removes all args required only by the 
  // parallel driver, since those are not understood
  // by the core
  class ParamPurgeWorker {
    vector< string > *purged_;

    public:
    ParamPurgeWorker(vector< string > *purged) : 
      purged_(purged)
    {
    }

    void matchValue(int argi, int argc, char **argv, Param *param){
    }

    void matchNoValue(int argi, int argc, char **argv, Param *param){
    }

    void noMatch(int argi, int argc, char **argv){
      (*purged_).push_back(argv[argi]);
    }
  };
 


  public:
  template<typename ParamType>
  void add(const string &description, ParamType *saveTo, ParamType defaultValue, bool hasValue){
    // since this parameter has a default value, it is
    // not compulsory
    params_[Environment::FlagPrefix_ + description] = new TypedParam<ParamType>(saveTo, defaultValue, hasValue, false);
  }

  template<typename ParamType>
  void add(const string &description, ParamType *saveTo, bool hasValue){
    // since this parameter doesn't have a default 
    // value, it is compulsory
    params_[Environment::FlagPrefix_ + description] = new TypedParam<ParamType>(saveTo, hasValue, true);
  }

  bool contains(const char *description){
    return params_.find(Environment::FlagPrefix_ + string(description)) != params_.end();
  }

  template<typename WorkerClass>
  bool traverse(int argc, char **argv, WorkerClass &worker){
    map<string, Param *>::iterator it;

    // check the parameters for which the user has
    // supplied values
    for(int i = 1; i < argc; ++i){
      string searchString(argv[i]);
      it = params_.find(searchString);

      //CkPrintf("Params: search %s found %d\n", searchString.c_str(), it != params_.end());

      if(it != params_.end()){
        Param *param = it->second;
        if(param->hasValue()){
          //CkPrintf("Params: matchValue %s\n", searchString.c_str());
          worker.matchValue(i, argc, argv, param);
          ++i;
        }
        else{
          //CkPrintf("Params: matchNoValue %s\n", searchString.c_str());
          worker.matchNoValue(i, argc, argv, param);
        }
        param->hasExtractedValue();
      }
      else{
        //CkPrintf("Params: noMatch %s\n", searchString.c_str());
        worker.noMatch(i, argc, argv);
      }
    }

    // set the value of the parameter, either to a
    // user-supplied value, or to the default;
    // complain if user has not supplied values for
    // compulsory parameters
    int nMissing = 0;
    for(it = params_.begin(); it != params_.end(); ++it){
      Param *param = it->second;
      //CkPrintf("Params: check param %s\n", it->first.c_str());
      if(param->isCompulsory() && !param->extractedValue()){
        //CkPrintf("Params: parameter `%s' is compulsory and unspecified: fail\n", it->first.c_str());
        nMissing++;
      }
      else {
        //CkPrintf("Params: parameter `%s' good\n", it->first.c_str());
        it->second->push();
      }
    }

    return (nMissing == 0);
  }

  bool process(int argc, char **argv){
    ParamExtractionWorker worker;
    return traverse(argc, argv, worker);
  }

  bool purge(int argc, char **argv, vector< string > &purgedArgs){
    ParamPurgeWorker worker(&purgedArgs);
    return traverse(argc, argv, worker);
  }
};

#endif // PARAMS_H
