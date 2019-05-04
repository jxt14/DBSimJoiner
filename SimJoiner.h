#ifndef __EXP2_SIMJOINER_H__
#define __EXP2_SIMJOINER_H__

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <algorithm>
#include <set>
#include <cmath>
#include <map>

const int SUCCESS = 0;
const int FAILURE = 1;


struct trie
{
	int qsize,ql,sl;
	trie* node[129];
	std::vector<int> qgram;
	trie()
	{
		qsize = 0;
		ql = -1;
		sl = -1;
		for(int i = 0; i <= 128; i++)node[i] = NULL;
		qgram.clear();
		qgram.push_back(0);
	}
};

typedef std::pair<int,trie*> pii;

template <typename IDType, typename SimType>
struct JoinResult {
    IDType id1;
    IDType id2;
    SimType s;
};

typedef JoinResult<unsigned, double> JaccardJoinResult;
typedef JoinResult<unsigned, unsigned> EDJoinResult;

const int SUCCESS = 0;
const int FAILURE = 1;

class SimJoiner {
public:
    SimJoiner();
    ~SimJoiner();
	std::vector<std::string> data1,data2;
	unsigned qlimit;
	trie* qroot;
	trie* jacroot;
	std::map<int, std::set<unsigned>> jacset;	

	pii qlists[200011];
	int querycheck[200011],occurtime[200011];
	int qsize;//the qgram size of the query string

	void insert(trie*, const char*, int, int);
	void search(trie*, const char*, int);
	void createIndex(const char *filename, std::vector<std::string>& datas);
	int CalCulateED(const char*, int, const char*, int, int);
    
	double timebuild,timequery,timedp;
    int joinJaccard(const char *filename1, const char *filename2, double threshold, std::vector<JaccardJoinResult> &result);
    int joinED(const char *filename1, const char *filename2, unsigned threshold, std::vector<EDJoinResult> &result);
};

#endif
