#include "crossvalidation.h"

uint CrossValidation::select(const vector<double> &aSolutionPath, const Database &aDatabase, PrefixSpan &aPrefix, LearnerSPP &aLearner, const vector<uint> &aOptions){
	uint tK = aOptions[0];
	uint tAve = aOptions[1];
	uint tLossType = aOptions[2];

	vector<double> tCnt(aSolutionPath.size() - 1, 0);

	vector<vector<Event> > tTransaction = aDatabase.get_transaction();
	vector<double> tY = aDatabase.get_y();
	uint tN = tY.size();
	if(tN < tK){
		cout << "error:CV,n<k" << '\n';
		exit(1);
	}

	cout << "CV now ";

	for(uint i = 0; i < tAve; ++i){
		if((i + 1) % 5 == 0){
			cout << i + 1 << flush;
		}else{
			cout << "." << flush;
		}

		for(uint j = 0; j < tK; ++j){
			mNumDiv.push_back((tN + j) / tK);
		}

		for(uint j = 0; j < tN; ++j){
			mTrainIdx.push_back(j);
		}

		for(uint j = 0; j < tN; ++j){
			uint k = j + (rand() % (tN - j));
			swap(mTrainIdx[j], mTrainIdx[k]);
		}

		for(uint j = 0; j < mNumDiv.back(); ++j){
			mTestIdx.push_back(mTrainIdx.back());
			mTrainIdx.pop_back();
		}
		mNumDiv.pop_back();

		k_fold(aSolutionPath, aPrefix, aLearner, tTransaction, tY, tK, tCnt, tLossType);

		mNumDiv.clear();
		mTrainIdx.clear();
		mTestIdx.clear();
	}

	vector<double>::iterator tIter;
	cout << '\n';
	switch(tLossType){
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
				cout << "λ[" << i + 1 << "]=" << aSolutionPath[i + 1] << " RMSE:" << sqrt(tCnt[i] / (double) (aOptions[1] * tN)) << '\n';
			}

			tIter = min_element(tCnt.begin(), tCnt.end());
			break;
		default:
			std::cout << "error:CV output" << '\n';
	}

	return distance(tCnt.begin(), tIter) + 1;
}

void CrossValidation::calc_accuracy(const vector<vector<double> > &aYHats, const vector<double> &aY, vector<double> &aCnt, const uint &aLossType){
	if(aCnt.size() != aYHats.size()){
		cout << "error:calc_accuracy" << aCnt.size() << ":" << aYHats.size() << '\n';
		exit(1);
	}

	for(uint j = 0; j < aYHats.size(); ++j){
		for(uint k = 0; k < aY.size(); ++k){
			if(aLossType == 1 && aY[k] * aYHats[j][k] > 0){
				aCnt[j]++;
			}else if(aLossType == 2){
				aCnt[j] += (aY[k] - aYHats[j][k]) * (aY[k] - aYHats[j][k]);
			}
		}
	}
}

void CrossValidation::next_train_test(){
	if(mNumDiv.empty()){
		cout << "error:CV next_train_test" << '\n';
		exit(1);
	}
	vector<uint> tTmp;
	copy(mTestIdx.begin(), mTestIdx.end(), back_inserter(tTmp));
	mTestIdx.clear();

	for(uint i = 0; i < mNumDiv.back(); ++i){
		mTestIdx.push_back(mTrainIdx.back());
		mTrainIdx.pop_back();
	}
	mNumDiv.pop_back();

	copy(mTrainIdx.begin(), mTrainIdx.end(), back_inserter(tTmp));
	mTrainIdx.clear();
	copy(tTmp.begin(), tTmp.end(), back_inserter(mTrainIdx));
}

void CrossValidation::k_fold(const vector<double> &aSolutionPath, PrefixSpan &aPrefix, LearnerSPP &aLearner, const vector<vector<Event>> &aTransaction, const vector<double> &aY, const uint &aK, vector<double> &aCnt, const uint &aLossType){

	vector<uint> tSVMOption = { (uint) aSolutionPath.size() };

	for(uint i = 0; i < aK; ++i){
		vector<vector<Event> > tTrainTransaction(mTrainIdx.size());
		vector<double> tTrainY(mTrainIdx.size());
		for(uint j = 0; j < mTrainIdx.size(); ++j){
			tTrainTransaction[j] = aTransaction[mTrainIdx[j]];
			tTrainY[j] = aY[mTrainIdx[j]];
		}

		vector<vector<Event> > tTestTransaction(mTestIdx.size());
		vector<double> tTestY(mTestIdx.size());
		for(uint j = 0; j < mTestIdx.size(); ++j){
			tTestTransaction[j] = aTransaction[mTestIdx[j]];
			tTestY[j] = aY[mTestIdx[j]];
		}

		aPrefix.init(tTrainTransaction, tTrainY);
		vector<vector<double> > tYHats = aLearner.get_all_predict(aPrefix, aSolutionPath, tTestTransaction, tSVMOption);
		calc_accuracy(tYHats, tTestY, aCnt, aLossType);

		if(i != aK - 1){
			next_train_test();
		}
	}

}
