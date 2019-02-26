#ifndef MODELSELECTOR_H
#define MODELSELECTOR_H

#include "prefixspan.h"
#include "learnerspp.h"
#include "database.h"

#include <vector>

using namespace std;

class ModelSelector{
public:
	virtual ~ModelSelector(){
	}
	virtual uint select(const vector<double> &aSolutionPath, const Database &aDatabase, PrefixSpan &aPrefix, LearnerSPP &aLearner, const vector<uint> &aOptions) = 0;
};

#endif
