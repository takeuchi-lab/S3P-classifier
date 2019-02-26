#include "database.h"

void Database::read(const char *aFilename){
	ifstream tFile(aFilename);
	if(!tFile){
		cerr << "Error: cannot open file" << endl;
		exit(1);
	}

	uint tItemSize = 0;
	double tLabel;
	string tLine;
	vector<Event> tSequence;

	while(getline(tFile, tLine)){
		tSequence.clear();
		stringstream ss1(tLine);
		ss1 >> tLabel;
		mY.push_back(tLabel);
		string eventstring;


		while(ss1 >> eventstring){
			stringstream ss2(eventstring);
			Itemset tItemset;
			vector<uint> tItem;
			string itemstring;
			int tmp;
			while(getline(ss2, itemstring, ':')){
				if(contain(itemstring, '(')){
					string tString;
					stringstream ss3(itemstring);
					while(getline(ss3, tString, '_')){
						if(contain(tString, '(')){
							if(contain(tString, ')')) break;

							tString.erase(tString.begin());
						} else if(contain(tString, ')')){
							tString.pop_back();
						}
						tItemset.push_back(stoi(tString));
					}
				}else{
				tmp = stoi(itemstring);
				uint tVal = (tmp < 0) ? 0xffffffff : tmp; // wild card
					tItem.push_back(tVal);
				}
			}

			if(!tItemSize){
				tItemSize = tItem.size();
			}else{
				if(tItemSize != tItem.size()){
					cerr << "Format Error: different Event Size at line: " << mTransaction.size() << ", event: " << tSequence.size() << endl;
					exit(-1);
				}
			}

			Event tEvent = make_pair(tItemset, tItem);
			tSequence.push_back(tEvent);
		}
		mTransaction.push_back(tSequence);
	}
}

vector<vector<Event> > Database::get_transaction() const{
	return mTransaction;
}

vector<double> Database::get_y() const{
	return mY;
}
