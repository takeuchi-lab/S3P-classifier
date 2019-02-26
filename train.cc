
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
#include <memory>

#include "database.h"
#include "prefixspan.h"

//learnerspp
#include "svmspp.h"
#include "lassospp.h"

//modelselector
#include "crossvalidation.h"

using namespace std;
using uint = unsigned int;

inline void exit_with_help(){
	cout << "-- Safe Pattern Pruning for sequential pattern mining --\n"
			"Usage: train [-options] input_file\n"
			"options:\n"
			"    -u : problem type (default 1)\n"
			"        1 -- regularization path computation for L1-reguralized L2-SVM\n"
			"        2 -- regularization path computation for Lasso\n"
			"        3 -- regularization path computation for L1-reguralized logistic regression(coming soon)\n"
			"    -t : learning lambda index(when do not cv) (default:most minimum lambda)\n"
			"    -m : minimum supportSum (default 1)\n"
			"    -L : maximum length of pattern (default 10)\n"
			"    -T : the number of regularization parameter (default 100)\n"
			"    -r : lambda min ratio (default 2.0)\n"
			"    -i : max outer iteration in optimization (default 1000000)\n"
			"    -f : frequency of calculate duality gap and convergence check (default 50)\n"
			"    -e : convergence criterion of duality gap (default 1e-6)\n"
			"    -F : name of reslut file (default output/result.csv)\n"
			"    -p : maximum interval of event (default 0|-1:none)\n"
			"    -c : whether to do cross validation (default 0:do not|1:do)\n"
			"    -k : k-fold of k(>1) (when do cross validation)(default:10)\n"
			"    -a : times to do cross validation(default:10)\n"
			"    -C : whether to do CloSpan (default 0:do not|1:do)\n"
			"    -M : whether to do Multiprocess  (default 0:do not|1:do)\n"
			"    -P : whether to do cross validation for Performance evaluation (default 0:do not|1:do)\n"
			"     	: Note that CV is always executed in multi process\n"
			"    -s : the Mode of counting supportSum(default 0)\n"
			"     	: 0 is 0 or 1 per record,1 is the number of pattern\n" << endl;

	exit(1);
}

void print_info(const uint &aLearnerType,const uint &aMaxpat){
	switch (aLearnerType) {
		case 1:
		cout << "Learner:L2SVM" << '\n';
		break;
		case 2:
		cout << "Learner:LASSO" << '\n';
		default:
		break;
	}
	cout << "Max pattern length:" << aMaxpat<<'\n';
}

void cv_for_lambda(const uint aLearnerType, const uint aMinsup, const uint aMaxpat, const int aInterval, const uint aSupMode, const uint aK, const uint aAve, uint& aTargetLambda, const uint aCloSpan, const uint aMP, const uint aPE, const Database aDatabase, const vector<vector<Event>> aTransaction, const vector<double> aY, LearnerSPP* aLearner, vector<double> aLambdas, const uint &aLossType){

	cout << "******Cross validation for select of lambda******" << endl;

	PrefixSpan tPrefix4cv(aMinsup, aMaxpat, aInterval, aSupMode, aCloSpan);
	tPrefix4cv.init(aTransaction, aY);

	vector<uint> tCVOptions = { aK, aAve, aLossType};

	cout << aK << "-fold cross validation" << endl;

	//measure CV
	chrono::system_clock::time_point start, end;
	start = chrono::system_clock::now();

	//lambda select
	if(aMP == 1){
		switch(aLearnerType){
			case 1: {
				CrossValidation tCV4svm;
				SVMSPP* tSvmLearner = dynamic_cast<SVMSPP *>(aLearner);

				if(tSvmLearner){
					aTargetLambda = tCV4svm.select<SVMSPP>(aLambdas, aDatabase, tPrefix4cv, *tSvmLearner, tCVOptions);
				}else{
					cout << "Error: dynamic_cast (pointer " << tSvmLearner << ")" << endl;
					exit(1);
				}
				break;
			}
			case 2: {
				CrossValidation tCV4lasso;
				LASSOSPP* tLassoLearner = dynamic_cast<LASSOSPP *>(aLearner);
				if(tLassoLearner){
					aTargetLambda = tCV4lasso.select<LASSOSPP>(aLambdas, aDatabase, tPrefix4cv, *tLassoLearner, tCVOptions);
				}else{
					cout << "Error: dynamic_cast" << endl;
					exit(1);
				}
				break;
			}
			default: {
				cerr << "Error unknown problem type" << endl;
				exit_with_help();
				break;
			}
		}
	}else{
		CrossValidation tCrossValidation;
		aTargetLambda = tCrossValidation.select(aLambdas, aDatabase, tPrefix4cv, *aLearner, tCVOptions);
	}

	end = chrono::system_clock::now();
	double elapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
	elapsed *= 0.001;

	cout << "*********selected**********" << endl;
	cout << "λ[" << aTargetLambda << "]=" << aLambdas[aTargetLambda] << endl;
	cout << "****************************" << endl;
	cout << "CVtime(sec):" << elapsed << endl;
	cout << "******Cross validation end******" << '\n' << endl;
	aTargetLambda++;

}

int main(int argc, char **argv){

	//default
	uint tLearnerType = 1;
	uint tMinsup = 1;
	uint tMaxpat = 10;
	uint tT = 100;
	uint tMaxiter = 1000000;
	uint tFreq = 50;
	double tEps = 1e-6;
	double tRatio = 2.0;
	int tInterval = 0;
	uint tSupMode = 0;
	uint tCV = 0;
	uint tK = 10;
	uint tAve = 10;
	uint tTargetLambda = tT;
	uint tCloSpan = 0;
	uint tMP = 0;
	uint tPE = 0;
	string tFilename = "result.csv";


	omp_set_nested(1);

	int tI;
	for(tI = 1; tI < argc; tI++){
		if(argv[tI][0] != '-'){
			break;
		}
		if(++tI >= argc){
			exit_with_help();
		}
		switch(argv[tI - 1][1]){
			case 't':
				tTargetLambda = atoi(argv[tI]);
				break;
			case 'm':
				tMinsup = atoi(argv[tI]);
				break;
			case 'L':
				tMaxpat = atoi(argv[tI]);
				break;
			case 'T':
				tT = atoi(argv[tI]);
				break;
			case 'r':
				tRatio = atof(argv[tI]);
				break;
			case 'i':
				tMaxiter = atoi(argv[tI]);
				break;
			case 'f':
				tFreq = atoi(argv[tI]);
				break;
			case 'e':
				tEps = atof(argv[tI]);
				break;
			case 'F':
				tFilename = argv[tI];
				break;
			case 'p':
				tInterval = atoi(argv[tI]);
				break;
			case 's':
				tSupMode = atoi(argv[tI]);
				break;
			case 'c':
				tCV = atoi(argv[tI]);
				break;
			case 'k':
				tK = atoi(argv[tI]);
				break;
			case 'a':
				tAve = atoi(argv[tI]);
				break;
			case 'u':
				tLearnerType = atoi(argv[tI]);
				break;
			case 'C':
				tCloSpan = atoi(argv[tI]);
				break;
			case 'M':
				tMP = atoi(argv[tI]);
				break;
			case 'P':
				tPE = atoi(argv[tI]);
				break;
			default:
				cerr << "Error unknown option: -" << argv[tI - 1][1] << endl;
				exit_with_help();
				break;
		}
	}

	if(tI >= argc){
		cerr << "Error please input filename" << endl;
		exit_with_help();
	}

	//read data
	Database tDatabase;
	tDatabase.read(argv[tI]);
	vector<vector<Event> > tTransaction = tDatabase.get_transaction();
	vector<double> tY = tDatabase.get_y();

	uint tPlusN = 0;
	for(auto it : tY){
		if(it > 0) tPlusN++;
	}

	//output information
	print_info(tLearnerType,tMaxpat);

	if(tPE == 1){
		uint tN = tY.size();
		if(tN < tK){
			cout << "error:CV,n<k" << '\n';
			exit(1);
		}

		cout << "******Cross validation for performance evaluation******" << endl;
		cout << tK << "-fold cross validation" << endl;

		vector<uint> tNumDiv;
		for(uint j = 0; j < tK; ++j){
			tNumDiv.push_back((tN + j) / tK);
		}

		vector<uint> tTrainIdx;
		for(uint j = 0; j < tN; ++j){
			tTrainIdx.push_back(j);
		}

		for(uint j = 0; j < tN; ++j){
			uint k = j + (rand() % (tN - j));
			swap(tTrainIdx[j], tTrainIdx[k]);
		}

		double tCorrect = 0;
		uint tTruePositive = 0;
		uint tFalsePositive = 0;

		//measure CV
		chrono::system_clock::time_point start, end;
		start = chrono::system_clock::now();

		#pragma omp parallel reduction(+:tCorrect, tTruePositive, tFalsePositive)
		{
			//make prefixspan
			PrefixSpan tPrefix(tMinsup, tMaxpat, tInterval, tSupMode, tCloSpan);

			#pragma omp for
			for(uint i = 0; i < tK; ++i){
				vector<uint> tTestIdx;
				vector<uint> tTrainIdx_private = tTrainIdx;
				uint tSumIdx = 0;
				for(uint j = 0; j < i; ++j){
					tSumIdx += tNumDiv[j];
				}

				copy(tTrainIdx_private.begin() + tSumIdx, tTrainIdx_private.begin() + tSumIdx + tNumDiv[i], back_inserter(tTestIdx));
				tTrainIdx_private.erase(tTrainIdx_private.begin() + tSumIdx, tTrainIdx_private.begin() + tSumIdx + tNumDiv[i]);

				vector<vector<Event> > tTrainTransaction(tTrainIdx_private.size());
				vector<double> tTrainY(tTrainIdx_private.size());
				for(uint j = 0; j < tTrainIdx_private.size(); ++j){
					tTrainTransaction[j] = tTransaction[tTrainIdx_private[j]];
					tTrainY[j] = tY[tTrainIdx_private[j]];
				}

				vector<vector<Event>> tTestTransaction(tTestIdx.size());
				vector<double> tTestY(tTestIdx.size());
				for(uint j = 0; j < tTestIdx.size(); ++j){
					tTestTransaction[j] = tTransaction[tTestIdx[j]];
					tTestY[j] = tY[tTestIdx[j]];
				}

				tPrefix.init(tTrainTransaction, tTrainY);

				LearnerSPP* tLearner = nullptr;

				switch(tLearnerType){
					case 1: {
						tLearner = new SVMSPP(tMaxiter, tFreq, tEps, tRatio);
						break;
					}
					case 2: {
						tLearner = new LASSOSPP(tMaxiter, tFreq, tEps, tRatio);
						break;
					}
					default: {
						cerr << "Error unknown problem type" << endl;
						exit_with_help();
						break;
					}
				}

				//calculate lambda max
				double tLamMax = tLearner->get_lambda_max(tPrefix);

				//calculate lambda sequence
				vector<double> tLambdas(tT);
				for(uint t = 0; t < tT; ++t){
					tLambdas[t] = tLamMax * pow(10, -tRatio * t / (tT - 1));
				}

				//cross validation
				if(tCV == 1){
					cv_for_lambda(tLearnerType, tMinsup, tMaxpat, tInterval, tSupMode, tK, tAve, tTargetLambda, tCloSpan, tMP, tPE, tDatabase, tTransaction, tY, tLearner, tLambdas,tLearnerType);
				}

				//learn
				vector<uint> tLearningOption = { tTargetLambda };
				tLearner->learn(tPrefix, tLambdas, tLearningOption);

				vector<double> tYHat = tLearner->predict(tPrefix, tTestTransaction);

				switch (tLearnerType) {
					case 1:
					for(uint j = 0; j < tTestY.size(); ++j){
						if(tTestY[j] * tYHat[j] > 0){
							tCorrect++;
							if(tTestY[j] > 0) tTruePositive++;
						}else if(tYHat[j] > 0){
							tFalsePositive++;
						}
					}
					break;
					case 2:
					for(uint j = 0; j < tTestY.size(); ++j){
						tCorrect += (tTestY[j]-tYHat[j])*(tTestY[j]-tYHat[j]);
					}
					break;
					default:
					break;
				}

				delete tLearner;
			}
		}

		end = chrono::system_clock::now();
		double elapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		elapsed *= 0.001;
		cout << "****************************" << endl;


		double tPrecision = (double)tTruePositive / (tTruePositive + tFalsePositive);
		double tRecall = (double)tTruePositive / tPlusN;
		double tF = 2 * tPrecision * tRecall / (tPrecision + tRecall);

		switch (tLearnerType) {
			case 1:
			cout << "Accuracy is " << tCorrect << " / " << tN << " = " << tCorrect / tN << endl;
			cout << "Precision is " << tPrecision << endl;
			cout << "Recall is " << tRecall << endl;
			cout << "F-measure is " << tF << endl;
			break;
			case 2:
			cout << "RMSE is " << sqrt(tCorrect / tN) << endl;
			default:
			break;
		}

		cout << "****************************" << endl;
		cout << "CVtime(sec):" << elapsed << endl;
		cout << "******Cross validation for performance evaluation end******" << '\n' << endl;

	}else{

		//make prefixspan
		PrefixSpan tPrefix(tMinsup, tMaxpat, tInterval, tSupMode, tCloSpan);
		tPrefix.init(tTransaction, tY);

		LearnerSPP* tLearner = nullptr;

		switch(tLearnerType){
			case 1: {
				tLearner = new SVMSPP(tMaxiter, tFreq, tEps, tRatio);
				break;
			}
			case 2: {
				tLearner = new LASSOSPP(tMaxiter, tFreq, tEps, tRatio);
				break;
			}
			default: {
				cerr << "Error unknown problem type" << endl;
				exit_with_help();
				break;
			}
		}

		//calculate lambda max
		double tLamMax = tLearner->get_lambda_max(tPrefix);

		//calculate lambda sequence
		vector<double> tLambdas(tT);
		for(uint t = 0; t < tT; ++t){
			tLambdas[t] = tLamMax * pow(10, -tRatio * t / (tT - 1));
		}

		//cross validation
		if(tCV == 1){
			cv_for_lambda(tLearnerType, tMinsup, tMaxpat, tInterval, tSupMode, tK, tAve, tTargetLambda, tCloSpan, tMP, tPE, tDatabase, tTransaction, tY, tLearner, tLambdas,tLearnerType);
		}

		//learn
		vector<uint> tLearningOption = { tTargetLambda };
		tLearner->learn(tPrefix, tLambdas, tLearningOption);
		delete tLearner;


		//write outcome
		tPrefix.printTree(tFilename);
	}
	return 0;
}
