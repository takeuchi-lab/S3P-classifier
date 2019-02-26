#include "prefixspan.h"

void PrefixSpan::printTree(string aFilename){
	ofstream tFile(aFilename);
	uint tCnt = 0;
	for(auto it : mActive){
		if(it->w != 0){
			tCnt++;
		}
	}

	tFile << "mPattern,supportSum,add lambda,w,list[id:sup]" << '\n';

	for(auto it : mActive){
		if(it->w != 0){
			tFile << it->patternStr << "," << it->supportSum << "," << it->addLambda << "," << it->w << ",[";
			for(uint j = 0; j < it->support.size(); j++)
				tFile << it->x[j] << ":" << it->support[j] << ",";
			tFile << "]" << '\n';
		}
	}

	cout << "Size of data :" << mN << '\n';
	cout << "Size of tree :" << mTree.size() << '\n';
	cout << "Size of Active :" << tCnt << '\n';
}

void PrefixSpan::add_history(uint aLambda){

	for(auto it : mActive){
		if(it->w != 0){
			if(it->addLambda == -1){
				it->addLambda = aLambda;
			}
		}else{
			it->addLambda = -1;
		}
	}
}

//for predict
int PrefixSpan::calcSup(const vector<Event> &aTransaction, const vector<Event> &aPattern) const{
	uint patSize = aPattern.size();
	uint k;
	int num = 0;
	int interval;

	if(aTransaction.size() < patSize) return num;

	for(uint i = 0; i <= aTransaction.size() - patSize; i++){

		interval = 0;

		k = 0;

		for(uint j = i; j < aTransaction.size(); j++){

			if(mInterval >= 0 && interval > mInterval){
				break;
			}
			if(aTransaction[j] == aPattern[k]){
				k++;
				interval = -1;
				if(k == patSize){
					//for 0/1
					if(mSupMode == 0){
						return 1;
					}
					k = 0;
					num++;
					i = j;
					break;
				}
			}
			interval++;
		}
	}
	return num;
}

double PrefixSpan::calcSup(uint aId, vector<Event> aPattern){
	switch(mSupMode){
		case 1:
			return calcPat(aId, aPattern);
			break;
		default:
			return 1;
			break;
	}
}

int PrefixSpan::calcPat(uint aId, vector<Event> aPattern){
	uint patSize = aPattern.size();
	uint k;
	int num = 0;
	int interval;

	if(mTransaction[aId].size() < patSize) return num;

	for(uint i = 0; i <= mTransaction[aId].size() - patSize; i++){

		interval = 0;

		k = 0;
		for(uint j = i; j < mTransaction[aId].size(); j++){

			if(mInterval >= 0 && interval > mInterval){
				break;
			}
			if(mTransaction[aId][j] == aPattern[k]){
				k++;
				interval = -1;
				if(k == patSize){
					k = 0;
					num++;
					i = j;
					break;
				}
			}
			interval++;
		}
	}
	return num;
}

void PrefixSpan::init(const vector<vector<Event> > aTransaction, const vector<double> aY){
	mTransaction = aTransaction;
	mY = aY;
	mItemSize = mTransaction[0][0].second.size();
	mN = mTransaction.size();
	if(mTransaction[0][0].first.size() > 0) mFlagItemsetExist = 1;

	Itemset tItemset;
	vector<uint> tWild(mItemSize, MAXVAL);
	Event tEvent(tItemset, tWild);
	mWildEvent = tEvent;

	if(mTransaction.empty() || mY.empty()){
		cout << "error:Data or label is empty." << '\n';
		exit(1);
	}
	mPattern.clear();
	mTree.clear();
	mActive.clear();
	mR.clear();

	//make root of tree
	map<Event, projectDB> tCounter;
	map<Event, uint> dupCheck;
	//When there is no restriction on the interval between events, duplication of initial nodes is not permitted
	if(mInterval < 0){
		for(uint i = 0; i < mN; ++i){
			for(uint j = 0, size = mTransaction[i].size(); j < size; ++j){
				Event tEvent;
				vector<uint> tItemset;
				tEvent.second = mTransaction[i][j].second;
				uint k = 0;
				if(mTransaction[i][j].first.size() > 0){
					for(; k < mTransaction[i][j].first.size(); ++k){
						tItemset.push_back(mTransaction[i][j].first[k]);
						tEvent.first = tItemset;
						if(dupCheck.find(tEvent) == dupCheck.end()){
							dupCheck[tEvent] = 0;

							Position tPosition;
							tPosition.sequence_index = i;
							tPosition.event_index = j;
							tPosition.itemset_index = k;
							tCounter[tEvent].push_back(tPosition);
						}
						tItemset.clear();
					}
				}else{
					tEvent.first = tItemset;
					if(dupCheck.find(tEvent) == dupCheck.end()){
						dupCheck[tEvent] = 0;
						Position tPosition;
						tPosition.sequence_index = i;
						tPosition.event_index = j;
						tPosition.itemset_index = k;
						tCounter[tEvent].push_back(tPosition);
					}
				}
			}
			dupCheck.clear();
		}
	}else{
		for(uint i = 0; i < mN; ++i){
			for(uint j = 0, size = mTransaction[i].size(); j < size; ++j){
				Event tEvent;
				vector<uint> tItemset;
				tEvent.second = mTransaction[i][j].second;
				uint k = 0;
				if(mTransaction[i][j].first.size() > 0){
					for(; k < mTransaction[i][j].first.size(); ++k){
						tItemset.push_back(mTransaction[i][j].first[k]);
						tEvent.first = tItemset;

						Position tPosition;
						tPosition.sequence_index = i;
						tPosition.event_index = j;
						tPosition.itemset_index = k;
						tCounter[tEvent].push_back(tPosition);
						tItemset.clear();
					}
				}else{
					tEvent.first = tItemset;
					Position tPosition;
					tPosition.sequence_index = i;
					tPosition.event_index = j;
					tPosition.itemset_index = k;
					tCounter[tEvent].push_back(tPosition);
				}
			}
		}

	}
	mPattern.clear();

	Node tRoot;
	mTree.push_back(tRoot);

	for(auto it = tCounter.begin(), end = tCounter.end(); it != end; ++it){
		mPattern.push_back(it->first);
		Node tNode;
		tNode.parent = mTree.begin();
		tNode.pattern = mPattern;
		tNode.patternStr = pat2str();
		tNode.supportSum = 0;
		tNode.pdb = it->second;
		tNode.w = 0;
		tNode.val = 0;
		tNode.closed = true;

		uint oid = MAXVAL;
		for(uint i = 0, size = tNode.pdb.size(); i < size; ++i){
			uint id = tNode.pdb[i].sequence_index;
			if(id != oid){
				double tSup = calcSup(id, mPattern);
				tNode.supportSum += tSup;
				tNode.support.push_back(tSup);
				tNode.x.push_back(id);
			}
			oid = id;
		}

		bool tAddFlag = true;
		if(mMinsup >= 1 && mMinsup > tNode.supportSum) tAddFlag = false;
		else if(mMinsup < 1 && mMinsup > (double) tNode.supportSum / mN) tAddFlag = false;

		if(tAddFlag){
			mTree.push_back(tNode);
			mTree.begin()->child.push_back(prev(mTree.end()));
		}
		mPattern.pop_back();
	}

}

uint PrefixSpan::isSubsequence(const vector<Event> aSeq1, const vector<Event> aSeq2){
	vector<Event> tShort_seq;
	vector<Event> tLong_seq;
	uint tSub_seq;

	if(aSeq1.size() < aSeq2.size()){
		tShort_seq = aSeq1;
		tLong_seq = aSeq2;
		tSub_seq = 1;
	}else if(aSeq1.size() > aSeq2.size()){
		tShort_seq = aSeq2;
		tLong_seq = aSeq1;
		tSub_seq = 2;
	}else if(mFlagItemsetExist){
		uint tSum1 = 0;
		uint tSum2 = 0;
		for(auto it : aSeq1){
			tSum1 += it.first.size();
		}
		for(auto it : aSeq2){
			tSum2 += it.first.size();
		}

		if(tSum1 > tSum2){
			tShort_seq = aSeq2;
			tLong_seq = aSeq1;
			tSub_seq = 2;
		}else if(tSum1 < tSum2){
			tShort_seq = aSeq1;
			tLong_seq = aSeq2;
			tSub_seq = 1;
		}else{
			return 0;
		}

	}else{
		return 0;
	}

	uint tCount = 0;
	uint diff_size = tLong_seq.size() - tShort_seq.size();

	for(uint i = 0; i <= diff_size; ++i){
		for(uint it = 0; it < tShort_seq.size(); ++it){
			if(isInclude(tLong_seq[it + i], tShort_seq[it]) || tShort_seq[it] == mWildEvent || tLong_seq[it + i] == mWildEvent){
				tCount++;
			}else{
				tCount = 0;
				break;
			}

			if(tCount == tShort_seq.size()){
				return tSub_seq;
			}
		}
	}

	return 0;
}

bool PrefixSpan::isInclude(const Event aEvent1, const Event aEvent2){
	if(aEvent1.first.size() < aEvent2.first.size()){
		return false;
	}else if(aEvent1.second != aEvent2.second){
		return false;
	}else{
		uint i = 0;
		bool tExistFlag = false;
		for(auto itr : aEvent2.first){
			tExistFlag = false;
			if(i == aEvent1.first.size()){
				return false;
			}
			while(itr >= aEvent1.first[i]){
				if(itr == aEvent1.first[i]){
					tExistFlag = true;
					i++;
					break;
				}else{
					i++;
				}

				if(i == aEvent1.first.size()){
					return false;
				}
			}

			if(!tExistFlag){
				return false;
			}
		}
		return true;
	}

}

bool PrefixSpan::isParent(const TreeIter aNodeIter, const TreeIter aChildIter){
	auto current = aChildIter->parent;
	while(current->pattern.size() >= aNodeIter->pattern.size()){
		if(current == aNodeIter) return true;

		current = current->parent;
	}
	return false;
}

uint PrefixSpan::sum_items_pdb(const projectDB aPdb){
	uint tSum = 0;
	for(auto itr : aPdb){
		uint id = itr.sequence_index;
		uint j = itr.event_index;
		if(mTransaction[id][j].first.size() > 0){
			tSum += mTransaction[id][j].first.size() - 1 - itr.itemset_index;
		}else{
			tSum++;
		}
		for(++j; j < mTransaction[id].size(); ++j){
			if(mTransaction[id][j].first.size() > 0){
				tSum += mTransaction[id][j].first.size();
			}else{
				tSum++;
			}
		}
	}
	return tSum;

}

/**
 * @return True if node can be cut
 */
bool PrefixSpan::calculate(Node& aNode){
	if(aNode.supportSum < mMinsup) return true;

	aNode.val = 0;
	double p = 0, m = 0;
	for(uint i = 0; i < aNode.x.size(); ++i){
		uint id = aNode.x[i];
		if(mType == 1 || mType == 2){ // svm
			if(mR[id] > 0){
				double val = mAlpha * mR[id] * mY[id] * aNode.support[i];
				aNode.val += val;
				(val > 0) ? p += val : m += val;
			}
		}else if(mType == 3 || mType == 4){ //lasso
			double val = mAlpha * mR[id] * aNode.support[i];
			aNode.val += val;
			(val > 0) ? p += val : m += val;
		}else if(mType == 5 || mType == 6){ //logistic
			double val = (mY[id] > 0) ? 1 / (1 + mR[id]) : 1 / (1 + mR[id]) - 1;
			val *= mAlpha * aNode.support[i];
			aNode.val += val;
			(val > 0) ? p += val : m += val;
		}

	}
	aNode.val = fabs(aNode.val);

	if(mType == 1 || mType == 3 || mType == 5){ // get_mMaxnode
		if(max(p, -m) < mMaxnode.val) return true;
	}else if(mType == 2 || mType == 4 || mType == 6){ // safe_screening
		if(max(p, -m) + mRadius * sqrt(aNode.supportSum) < 1) return true;
	}

	return false;
}

/**
 * @fn
 * Used only within project_new
 * Counting support
 * @return True if node can be a leaf
 */
bool PrefixSpan::calculate_new(Node &aNode){
	uint oid = MAXVAL;
	double p = 0, m = 0;
	for(uint i = 0, size = aNode.pdb.size(); i < size; ++i){
		uint id = aNode.pdb[i].sequence_index;
		if(oid != id){
			double tSup = calcSup(id, aNode.pattern);
			aNode.supportSum += tSup;
			aNode.support.push_back(tSup);
			aNode.x.push_back(id);
			if(mType == 1 || mType == 2){ // svm
				if(mR[id] > 0){
					double val = mAlpha * mR[id] * mY[id] * tSup;
					aNode.val += val;
					(val > 0) ? p += val : m += val;
				}
			}else if(mType == 3 || mType == 4){ // lasso
				double val = mAlpha * mR[id] * tSup;
				aNode.val += val;
				(val > 0) ? p += val : m += val;
			}else if(mType == 5 || mType == 6){ // logistic
				double val = (mY[id] > 0) ? 1 / (1 + mR[id]) : 1 / (1 + mR[id]) - 1;
				val *= mAlpha * tSup;
				aNode.val += val;
				(val > 0) ? p += val : m += val;
			}
		}
		oid = id;
	}

	aNode.val = fabs(aNode.val);

	if(mMinsup > aNode.supportSum) return true;
	else if(mMinsup < 1 && mMinsup > (double) aNode.supportSum / mN) return true;

	if(mType == 1 || mType == 3 || mType == 5){
		if(max(p, -m) < mMaxnode.val) return true;
	}else if(mType == 2 || mType == 4 || mType == 6){
		if(max(p, -m) + mRadius * sqrt(aNode.supportSum) < 1) return true;
	}

	return false;
}

void PrefixSpan::project(const TreeIter aNodeIter){
	projectDB tPdb = aNodeIter->pdb;

	// scan projected database

	//I-Extension
	if(mFlagItemsetExist){
		map<Event, projectDB> tCounter;
		for(uint i = 0, size = tPdb.size(); i < size; ++i){
			uint id = tPdb[i].sequence_index;
			uint j = tPdb[i].event_index;
			if(mTransaction[id][j].first.size() - 1 > tPdb[i].itemset_index){
				uint k = tPdb[i].itemset_index + 1;
				for(; k < mTransaction[id][j].first.size(); ++k){
					Event tEvent;
					vector<uint> tItemset = mPattern.back().first;
					tItemset.push_back(mTransaction[id][j].first[k]);
					tEvent.first = tItemset;
					tEvent.second = mTransaction[id][j].second;
					Position tPos;
					tPos.sequence_index = id;
					tPos.event_index = j;
					tPos.itemset_index = k;
					tCounter[tEvent].push_back(tPos);
				}
			}
		}

		// project: next event
		Event tSaveEvent = mPattern.back();
		for(auto it = tCounter.begin(), end = tCounter.end(); it != end; ++it){
			mPattern.pop_back();
			mPattern.push_back(it->first);
			project_new(it->second, aNodeIter);
		}
		mPattern.pop_back();
		mPattern.push_back(tSaveEvent);
	}

	if(mPattern.size() < mMaxpat){
		//S-Extension
		map<Event, projectDB> tCounter;

		if(mInterval < 0){
			map<Event, uint> dupCheck;
			for(uint i = 0, size = tPdb.size(); i < size; ++i){
				uint id = tPdb[i].sequence_index;
				uint trsize = mTransaction[id].size();
				uint j = tPdb[i].event_index + 1;
				for(; j < trsize; j++){
					Event tEvent;
					vector<uint> tItemset;
					tEvent.second = mTransaction[id][j].second;
					uint k = 0;

					Position tPosition;
					tPosition.sequence_index = id;
					tPosition.event_index = j;

					if(mTransaction[id][j].first.size() > 0){
						for(; k < mTransaction[id][j].first.size(); ++k){
							tItemset.push_back(mTransaction[id][j].first[k]);
							tEvent.first = tItemset;
							if(dupCheck.find(tEvent) == dupCheck.end()){
								dupCheck[tEvent] = 0;
								tPosition.itemset_index = k;
								tCounter[tEvent].push_back(tPosition);
							}
							tItemset.clear();
						}
					}else{
						tEvent.first = tItemset;
						if(dupCheck.find(tEvent) == dupCheck.end()){
							dupCheck[tEvent] = 0;
							tPosition.itemset_index = k;
							tCounter[tEvent].push_back(tPosition);
						}
					}
				}
				dupCheck.clear();
			}
		}else{
			vector<Position> dupCheck;
			for(uint i = 0, size = tPdb.size(); i < size; ++i){
				uint id = tPdb[i].sequence_index;
				uint trsize = mTransaction[id].size();
				uint j = tPdb[i].event_index + 1;
				uint j_tmp = j;

				for(; j < trsize; j++){
					if((int) j - j_tmp > mInterval){
						break;
					}

					Event tEvent;
					vector<uint> tItemset;
					tEvent.second = mTransaction[id][j].second;
					uint k = 0;

					Position tPos;
					tPos.sequence_index = id;
					tPos.event_index = j;

					if(mTransaction[id][j].first.size() > 0){
						for(; k < mTransaction[id][j].first.size(); ++k){
							tItemset.push_back(mTransaction[id][j].first[k]);
							tEvent.first = tItemset;
							tPos.itemset_index = k;
							if(find(dupCheck.begin(), dupCheck.end(), tPos) == dupCheck.end()){
								dupCheck.push_back(tPos);
								tCounter[tEvent].push_back(tPos);
							}
							tItemset.clear();
						}
					}else{
						tPos.itemset_index = k;
						tEvent.first = tItemset;
						if(find(dupCheck.begin(), dupCheck.end(), tPos) == dupCheck.end()){
							dupCheck.push_back(tPos);
							tCounter[tEvent].push_back(tPos);
						}
					}
				}
			}
		}

		// project: next event
		for(auto it = tCounter.begin(), end = tCounter.end(); it != end; ++it){
			mPattern.push_back(it->first);
			project_new(it->second, aNodeIter);
			mPattern.pop_back();
		}
	}

}

void PrefixSpan::project_new(projectDB& aPdb, const TreeIter aParent){
	Node tNode;
	tNode.parent = aParent;
	tNode.pattern = mPattern;
	tNode.patternStr = pat2str();
	tNode.supportSum = 0;
	tNode.pdb = aPdb;
	tNode.w = 0;
	tNode.val = 0;
	tNode.closed = true;

	TreeIter tCurrent = mTree.insert(mTree.end(), tNode);

	bool tFlag = calculate_new(*tCurrent);

	if(mCloSpan == 1){
		if(tCurrent->parent->supportSum == tCurrent->supportSum){
			tCurrent->parent->closed = false;
		}
		mFlagCScheckEnd = 0;

		checkProjectedDB(tCurrent);
		if(mFlagCScheckEnd == 2){
			mTree.pop_back();
			return;
		}
	}

	aParent->child.push_back(tCurrent);

	if(tFlag){
		return;
	}

	if(mType == 1 || mType == 3 || mType == 5){ // get_mMaxnode
		if(tCurrent->val > mMaxnode.val){
			mMaxnode = *tCurrent;
		}
	}else if(mType == 2 || mType == 4 || mType == 6){ // safe_screening
	//SPPC
		double score = tCurrent->val + mRadius * sqrt(tCurrent->supportSum);
		if(score >= 1){
			bool tAdd_flag = true;
			for(uint i = 0; i < mActive.size(); ++i){
				if(mActive[i]->supportSum == tCurrent->supportSum){
					uint tWhichSub = isSubsequence(tCurrent->pattern, mActive[i]->pattern);
					if(tWhichSub == 1){
						tAdd_flag = false;
						break;
					}else if(tWhichSub == 2){
						mActive.erase(mActive.begin() + i);
						i--;
						break;
					}
				}
			}

			if(tAdd_flag){
				mActive.push_back(tCurrent);
			}

		}
	}

	if(mFlagCScheckEnd == 0){
		project(tCurrent);
	}
}

void PrefixSpan::checkProjectedDB(const TreeIter aCurrentIter){
	uint tSumSequence = 0;
	for(auto id : aCurrentIter->x){
		for(uint j = 0; j < mTransaction[id].size(); ++j){
			if(mTransaction[id][j].first.size() > 0){
				tSumSequence += mTransaction[id][j].first.size();
			}else{
				tSumSequence += 1;
			}
		}
	}

	uint tKey = sum_items_pdb(aCurrentIter->pdb) + tSumSequence;

	if(mCheckPDB.find(tKey) != mCheckPDB.end()){
		for(auto itr = mCheckPDB[tKey].begin(); itr != mCheckPDB[tKey].end();){
			if(aCurrentIter->supportSum == (*itr)->supportSum){
				uint tWhichSub = isSubsequence(aCurrentIter->pattern, (*itr)->pattern);
				if(tWhichSub == 1){
					mFlagCScheckEnd = 2;
					return;
				}else if(tWhichSub == 2){
					aCurrentIter->child = (*itr)->child;
					(*itr)->parent->child.erase(find((*itr)->parent->child.begin(), (*itr)->parent->child.end(), *itr));
					mTree.erase(*itr);
					itr = mCheckPDB[tKey].erase(itr);

					for(auto it : aCurrentIter->child){
						it->parent = aCurrentIter;
					}

					for(auto it : aCurrentIter->child){
						childPatternUpdate(it);
					}
					mFlagCScheckEnd = 1;
					return;
				}
			}

			itr++;
		}
	}

	mCheckPDB[tKey].push_back(aCurrentIter);
	return;
}

void PrefixSpan::childPatternUpdate(const TreeIter aChildIter){
	Event tChildLast = aChildIter->pattern.back();

	Event tParentLast = mPattern.back();
	bool tIE = false;
	if(tChildLast.first.size() > tParentLast.first.size()){
		mPattern.pop_back();
		tIE = true;
	}
	mPattern.push_back(tChildLast);
	aChildIter->pattern = mPattern;
	aChildIter->patternStr = pat2str();
	for(auto it : aChildIter->child){
		childPatternUpdate(it);
	}
	mPattern.pop_back();
	if(tIE) mPattern.push_back(tParentLast);
}

void PrefixSpan::searchTree_forMaxnode(const TreeIter aNodeIter){
	mPattern = aNodeIter->pattern;
	bool flag = calculate(*aNodeIter);
	if(!flag){
		if(aNodeIter->closed){
			if(aNodeIter->val > mMaxnode.val) mMaxnode = *aNodeIter;
		}

		if(aNodeIter->child.size() == 0){
			project(aNodeIter);
		}else{
			for(auto it : aNodeIter->child){
				searchTree_forMaxnode(it);
			}
		}
	}
}

void PrefixSpan::searchTree_forSafeScreening(const TreeIter aNodeIter){
	mPattern = aNodeIter->pattern;
	bool flag = calculate(*aNodeIter);
	if(!flag){
		if(aNodeIter->closed){
			double score = aNodeIter->val + mRadius * sqrt(aNodeIter->supportSum);
			if(score >= 1){
				bool tAdd_flag = true;
				for(uint i = 0; i < mActive.size(); ++i){
					if(mActive[i]->supportSum == aNodeIter->supportSum){
						uint tWhichSub = isSubsequence(aNodeIter->pattern, mActive[i]->pattern);
						if(tWhichSub == 1){
							tAdd_flag = false;
							break;
						}else if(tWhichSub == 2){
							mActive.erase(mActive.begin() + i);
							i--;
							break;
						}
					}
				}

				if(tAdd_flag){
					mActive.push_back(aNodeIter);
				}

			}else{
				aNodeIter->w = 0;
			}
		}

		if(aNodeIter->child.size() == 0){
			project(aNodeIter);
		}else{
			for(auto it : aNodeIter->child){
				searchTree_forSafeScreening(it);
			}
		}
	}else{
		aNodeIter->w = 0;
	}

}

void PrefixSpan::get_maxnode(const vector<double> &_v, int solver){
	switch(solver){
		case 1: // svm
			mType = 1;
			break;

		case 2: //lasso
			mType = 3;
			break;
		case 3: //logistic
			mType = 5;
			break;
		default:
			cerr << "Error: unknown solver type" << endl;
			exit(1);
			break;
	}
	mR = _v;
	mAlpha = 1;
	mMaxnode.val = 0;
	mPattern.clear();

	TreeIter root = mTree.begin();

	for(auto it : root->child){
		mPattern = (it->pattern);
		searchTree_forMaxnode(it);
	}

}

void PrefixSpan::safe_screening(vector<double> &_v, double _alpha, double _radius, int solver){
	switch(solver){
		case 1: // svm
			mType = 2;
			break;
		case 2: // lasso
			mType = 4;
			break;
		case 3: //logistic
			mType = 6;
			break;
		default:
			cerr << "Error: unknown solver type" << endl;
			exit(1);
			break;
	}
	mR = _v;
	mAlpha = _alpha;
	mRadius = _radius;
	mPattern.clear();
	mActive.clear();

	TreeIter root = mTree.begin();

	for(auto it : root->child){
		mPattern = (it->pattern);
		searchTree_forSafeScreening(it);
	}
}
