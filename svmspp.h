#ifndef SVMSPP_H
#define SVMSPP_H

#include "prefixspan.h"
#include "database.h"
#include "learnerspp.h"
#include <vector>

using namespace std;

//L1-reguralized L2-SVM
class SVMSPP: public LearnerSPP{
private:

	uint mN;

	double mBias;
	vector<double> mR;
	uint mT;
	uint mMaxIter;
	uint mFreq;
	double mEps;
	double mRatio;
	double clac_sup(const vector<Event> &aSequence, const vector<Event> &aPattern, const uint aSupMode);

public:
	SVMSPP(uint aMaxIter, uint aFreq, double aEps, double aRatio){
		mMaxIter = aMaxIter;
		mFreq = aFreq;
		mEps = aEps;
		mRatio = aRatio;
	}
	;
	//Learning with aLambdas[i - 1] is done by inserting a number from 1 to aLambdas.size() in aOption[0]
	virtual void learn(PrefixSpan &aPrefix, const vector<double> &aLambdas, const vector<uint> &aOptions);
	//Calculate and return the value of λmax
	virtual double get_lambda_max(PrefixSpan &aPrefix);
	//Perform prediction using the model currently being learned
	virtual vector<double> predict(const PrefixSpan &aPrefix, const vector<vector<Event> > &aTransaction) const;
	//Returns the predicted value for all solution paths
	virtual vector<vector<double> > get_all_predict(PrefixSpan &aPrefix, const vector<double> &aLambdas, const vector<vector<Event> > &aTransaction, const vector<uint> &aOptions);
};

#endif
