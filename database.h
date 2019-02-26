#ifndef DATABASE_H
#define DATABASE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

using uint = unsigned int;
using Itemset = vector<uint>;
using Event = pair<Itemset, vector<uint>>;

class Database{
private:
	vector<vector<Event> > mTransaction;
	vector<double> mY;
public:
	void read(const char *aFilename);
	vector<vector<Event> > get_transaction() const;
	vector<double> get_y() const;

	template<class T, class U>

	bool contain(const std::basic_string<T>& s, const U& v){
		return s.find(v) != std::basic_string<T>::npos;
	}
};

#endif
