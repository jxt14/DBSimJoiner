#include "SimJoiner.h"

using namespace std;

SimJoiner::SimJoiner() {
}

SimJoiner::~SimJoiner() {
}

inline int my_abs(int x)
{
	if (x > 0) return x;
	else return -x;
}

inline int my_min(int a, int b)
{
	if (a < b) return a;
	else return b;
}

inline int my_max(int a, int b)
{
	if (a > b) return a;
	else return b;
}

inline int min_3(int a, int b, int c)
{
	return my_min(my_min(a,b), c);
}

int SimJoiner::CalCulateED(const char* s1, int m, const char* s2, int n, int threshold)
{
    if (my_abs(m - n) > threshold)
        return threshold + 1;
    int f[m+1][n+1];
    for (int i = 0; i <= my_min(threshold, m); i++)f[i][0] = i;
    for (int i = 0; i <= my_min(threshold, n); i++)f[0][i] = i;
    for (int i = 1; i <= m; i++)
    {
        int limmin = my_max(i - threshold, 1);
        int limmax = my_min(i + threshold, n);
        for (int j = limmin; j <= limmax; j++)
        {
			int f1,f2;
			if (my_abs(i - 1 - j) > threshold) f1 = threshold + 1;
			else f1 = f[i-1][j];
			if (my_abs(i - j + 1) > threshold) f2 = threshold + 1;
			else f2 = f[i][j-1];
            f[i][j] = min_3(f1 + 1, f2 + 1, f[i-1][j-1] + (s1[i-1] != s2[j-1]));
        }
    }
    return f[m][n];

}

void SimJoiner::createIndex(const char *filename, vector<string>& datas)
{
	char c;
	bool p;
	int len;
	ifstream fp;
	p = false;

    datas.clear();

	fp.open(filename);
	
	string temp;
	while (getline(fp, temp)) {
		len = temp.length();
		datas.push_back(temp);
	}
	fp.close();

}

void SimJoiner::insert(trie* rt, const char* s, int id, int tlim)
{
	for (int i = 0; i < tlim; i++){
		if (rt -> node[(int)s[i]] == NULL){
			rt->node[(int)s[i]] = new trie();
		}
		rt = rt -> node[(int)s[i]];
	}
	if (rt->ql != id){
		rt->qgram.push_back(id);
		rt->qsize++;
		rt->ql = id;

	}
}

void SimJoiner::search(trie* rt, const char* s, int tlim)
{
	for (int i = 0; i < tlim; i++){
		if (rt->node[(int)s[i]] == NULL) return;
		rt = rt -> node[(int)s[i]];
	}
	int tid;

}

int SimJoiner::joinJaccard(const char *filename1, const char *filename2, double threshold, vector<JaccardJoinResult> &result) {
    result.clear();
    return SUCCESS;
}

int SimJoiner::joinED(const char *filename1, const char *filename2, unsigned threshold, vector<EDJoinResult> &result) {
    result.clear();
    createIndex(filename1, data1);
    createIndex(filename2, data2);


    return SUCCESS;
}
