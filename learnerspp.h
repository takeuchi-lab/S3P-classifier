#ifndef LEARNERSPP_H
#define LEARNERSPP_H

#include "prefixspan.h"

using namespace std;

class LearnerSPP{
public:
	virtual ~LearnerSPP(){
	}
	virtual void learn(PrefixSpan &aPrefixspan, const vector<double> &aLambdas, const vector<uint> &aOptions) = 0;
	virtual double get_lambda_max(PrefixSpan &aPrefix)=0;
	virtual vector<double> predict(const PrefixSpan &aPrefix, const vector<vector<Event> > &aTransaction) const=0;
	virtual vector<vector<double> > get_all_predict(PrefixSpan &aPrefix, const vector<double> &aLambdas, const vector<vector<Event> > &aTransaction, const vector<uint> &aOptions)=0;
};

#endif
