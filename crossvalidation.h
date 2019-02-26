#ifndef CROSSVALIDATION_H
#define CROSSVALIDATION_H

#include "modelselector.h"
#include <random>

#include "svmspp.h"
#include "lassospp.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

class CrossValidation: public ModelSelector{
private:
	vector<uint> mTrainIdx;
	vector<uint> mTestIdx;
	vector<uint> mNumDiv;

public:
	CrossValidation(){};
	void calc_accuracy(const vector<vector<double> > &aYHats, const vector<double> &aY, vector<double> &aCnt,const uint &aLossType);
	void next_train_test();

	template<class T> uint select(const vector<double> &aSolutionPath, const Database &aDatabase, PrefixSpan &aPrefix, T &aLearner, const vector<uint> &aOptions){
		uint tK = aOptions[0];
		uint tAve = aOptions[1];
		uint tLossType = aOptions[2];

		vector<vector<Event> > tTransaction = aDatabase.get_transaction();
		vector<double> tY = aDatabase.get_y();
		uint tN = tY.size();
		if(tN < tK){
			cout << "error:CV,n<k" << '\n';
			exit(1);
		}

		vector<double> tCnt(aSolutionPath.size() - 1, 0);

		vector<uint> tNumDiv;
		for(uint j = 0; j < tK; ++j){
			tNumDiv.push_back((tN + j) / tK);
		}

		#pragma omp parallel
		{
			vector<double> tCnt_private(aSolutionPath.size() - 1, 0);

			#pragma omp for
			for(uint i = 0; i < tAve; ++i){
				vector<uint> tTrainIdx;
				for(uint j = 0; j < tN; ++j){
					tTrainIdx.push_back(j);
				}

				for(uint j = 0; j < tN; ++j){
					uint k = j + (rand() % (tN - j));
					swap(tTrainIdx[j], tTrainIdx[k]);
				}

				k_fold(aSolutionPath, aPrefix, aLearner, tTransaction, tY, tK, tCnt_private, tTrainIdx, tNumDiv,tLossType);
			}

			#pragma omp critical
			{
				for(uint j = 0; j < tCnt.size(); ++j){
					tCnt[j] += tCnt_private[j];
				}
			}

		}

		vector<double>::iterator tIter;
		cout << '\n';
		switch (tLossType) {
      // L2SVM
			case 1:
			for(uint i = 0; i < tCnt.size(); ++i){
				cout << "λ[" << i + 1 << "]=" << aSolutionPath[i + 1] << " " << tCnt[i] << "/" << aOptions[1] * tN << "=" << tCnt[i] / (double) (aOptions[1] * tN) << '\n';
			}

			tIter = max_element(tCnt.begin(), tCnt.end());
			break;
      //Lasso
			case 2:
			for(uint i = 0; i < tCnt.size(); ++i){
				cout << "λ[" << i + 1 << "]=" << aSolutionPath[i + 1] << " RMSE:" <<sqrt(tCnt[i] / (double) (aOptions[1] * tN)) << '\n';
			}

			tIter = min_element(tCnt.begin(), tCnt.end());
			break;
			default:
			std::cout << "error:CV output" << '\n';
		}




		return distance(tCnt.begin(), tIter) + 1;
	}

	template<class T> void k_fold(const vector<double> &aSolutionPath, PrefixSpan &aPrefix, T &aLearner, const vector<vector<Event>> &aTransaction, const vector<double> &aY, const uint &aK, vector<double> &aCnt, vector<uint> aTrainIdx, vector<uint> aNumDiv,const uint &aLossType){

		vector<uint> tSVMOption = { (uint) aSolutionPath.size() };

		#pragma omp parallel
		{
			PrefixSpan tPrefix_private = aPrefix;

			T tLearner_private = aLearner;

			#pragma omp for
			for(uint i = 0; i < aK; ++i){
				vector<uint> tTestIdx;
				vector<uint> tTrainIdx = aTrainIdx;
				uint tSumIdx = 0;
				for(uint j = 0; j < i; ++j){
					tSumIdx += aNumDiv[j];
				}

				copy(tTrainIdx.begin() + tSumIdx, tTrainIdx.begin() + tSumIdx + aNumDiv[i], back_inserter(tTestIdx));
				tTrainIdx.erase(tTrainIdx.begin() + tSumIdx, tTrainIdx.begin() + tSumIdx + aNumDiv[i]);

				vector<vector<Event> > tTrainTransaction(tTrainIdx.size());
				vector<double> tTrainY(tTrainIdx.size());
				for(uint j = 0; j < tTrainIdx.size(); ++j){
					tTrainTransaction[j] = aTransaction[tTrainIdx[j]];
					tTrainY[j] = aY[tTrainIdx[j]];
				}

				vector<vector<Event>> tTestTransaction(tTestIdx.size());
				vector<double> tTestY(tTestIdx.size());
				for(uint j = 0; j < tTestIdx.size(); ++j){
					tTestTransaction[j] = aTransaction[tTestIdx[j]];
					tTestY[j] = aY[tTestIdx[j]];
				}

				tPrefix_private.init(tTrainTransaction, tTrainY);

				vector<vector<double> > tYHats = tLearner_private.get_all_predict(tPrefix_private, aSolutionPath, tTestTransaction, tSVMOption);

				calc_accuracy(tYHats, tTestY, aCnt,aLossType);
			}
		}

	}

	void k_fold(const vector<double> &aSolutionPath, PrefixSpan &aPrefix, LearnerSPP &aLearner, const vector<vector<Event> > &aTransaction, const vector<double> &aY, const uint &aK, vector<double> &aCnt,const uint &aLossType);

	uint select(const vector<double> &aSolutionPath, const Database &aDatabase, PrefixSpan &aPrefix, LearnerSPP &aLearner, const vector<uint> &aOptions);
};

#endif
