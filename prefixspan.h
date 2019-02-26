#ifndef PREFIXSPAN_H
#define PREFIXSPAN_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <algorithm>
#include <unordered_map>

#include "database.h"

using namespace std;

using uint = unsigned int;
using uchar = unsigned char;

class PrefixSpan{
private:


	const uint MAXVAL = 0xffffffff;
	uint mItemSize;
	double mMinsup;
	uint mMaxpat;
	int mInterval;
	uint mSupMode;
	uint mCloSpan;


	int mType;
	vector<double> mR;
	double mAlpha;
	double mRadius;

	vector<vector<Event> > mTransaction;
	vector<Event> mPattern;


	uint mFlagItemsetExist = 0;

	Event mWildEvent;

	string pat2str(void){
		stringstream ss;
		for(uint i = 0; i < mPattern.size(); ++i){
			ss << (i ? " " : "");
			if(mPattern[i].first.size() > 0){
				ss << "(";
				for(uint j = 0; j < mPattern[i].first.size(); ++j){
					ss << mPattern[i].first[j];
					ss << (j + 1 < mPattern[i].first.size() ? "_" : "");
				}
				ss << ")";
				if(mPattern[i].second.size() > 0) ss << ":";
			}
			for(uint j = 0; j < mItemSize; j++){
				if(mPattern[i].second[j] == MAXVAL) ss << "*";
				else ss << mPattern[i].second[j];
				ss << (j + 1 < mItemSize ? ":" : "");
			}
		}
		return ss.str();
	}

	string pat2str(const vector<Event> aPattern){
		stringstream ss;
		for(uint i = 0; i < aPattern.size(); ++i){
			ss << (i ? " " : "");
			if(aPattern[i].first.size() > 0){
				ss << "(";
				for(auto it : aPattern[i].first){
					ss << it;
					if(it != aPattern[i].first.back()) ss << "_";
				}
				ss << ")";
				if(aPattern[i].second.size() > 0) ss << ":";
			}
			for(auto it : aPattern[i].second){
				if(it == MAXVAL) ss << "*";
				else ss << it;
				if(it != aPattern[i].second.back()) ss << ":";
			}
		}
		return ss.str();
	}


	double calcSup(uint aId, vector<Event> aPattern);

	int calcPat(uint aId, vector<Event> aPattern);

	/*
	 * Check whether the input sequence is inclusive.
	 * Returns which of the first or second is shorter.
	 * @return 1 or 2 or 0(Not inclusive or same sequence)
	 */
	uint isSubsequence(const vector<Event> aSeq1, const vector<Event> aSeq2);

	bool isInclude(const Event aEvent1, const Event aEvent2);

public:

	struct Position{
		uint sequence_index;
		uint event_index;
		uint itemset_index;
		bool operator==(const Position& x){
			return (sequence_index == x.sequence_index && event_index == x.event_index && itemset_index == x.itemset_index);
		}
		bool operator!=(const Position& x){
			return (sequence_index != x.sequence_index || event_index != x.event_index || itemset_index != x.itemset_index);
		}
	};

	using projectDB = vector< Position >;


	struct Node{
		list<Node>::iterator parent;
		list<list<Node>::iterator> child;
		vector<Event> pattern;
		string patternStr;
		double supportSum;
		vector<double> support;
		projectDB pdb;
		vector<uint> x;
		double w;
		double val;
		int addLambda = -1;
		bool closed; //CloSpan
	};


	list<Node> mTree;
	uint mN;
	vector<double> mY;
	Node mMaxnode;
	using TreeIter = list<Node>::iterator;
	vector<TreeIter> mActive;

	uint mFlagCScheckEnd = 0;
	unordered_map<uint, list<TreeIter>> mCheckPDB;

	PrefixSpan(double aMinsup, uint aMaxpat, int aInterval,uint aSupMode, uint aCloSpan) {
		mMinsup = aMinsup;
		mMaxpat = aMaxpat;
		mInterval = aInterval;
		mSupMode = aSupMode;
		mCloSpan = aCloSpan;
	};

	void init(const vector<vector<Event> > aTransaction, const vector<double> aY);
	void get_maxnode(const vector<double> &_v, int solver);
	void safe_screening(vector<double> &_v, double _alpha, double _radius, int solver);
	void printTree(string aFilename);
	void add_history(uint aLambda);

	void project(const TreeIter aNodeIter);
	void project_new(projectDB& aPdb, const TreeIter aParent);

	bool calculate(Node& aNode);
	bool calculate_new(Node& aNode);

	int calcSup(const vector<Event> &aTransaction, const vector<Event> &aPattern) const;

	bool isParent(const TreeIter aNodeIter, const TreeIter aChildIter);

	void checkProjectedDB(const TreeIter aCurrentIter);

	uint sum_items_pdb(const projectDB aPdb);

	void childPatternUpdate(const TreeIter aChildIter);

	void searchTree_forMaxnode(const TreeIter aNodeIter);
	void searchTree_forSafeScreening(const TreeIter aNodeIter);
};

#endif
