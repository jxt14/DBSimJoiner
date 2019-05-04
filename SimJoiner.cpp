#include "SimJoiner.h"

using namespace std;

SimJoiner::SimJoiner() {
    qlimit = 3;
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
 //   cout << "input " << filename << endl;

	fp.open(filename);
	
	string temp;
	while (getline(fp, temp)) {
		len = temp.length();
		datas.push_back(temp);
 //       cout << temp << endl;
	}
	fp.close();

//    cout << "input end" << endl;
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
	    rt->ql = id;
        rt->qgram.push_back(id);
	    rt->qsize++;
    }
}

trie* SimJoiner::search(trie* rt, const char* s, int tlim)
{
    for (int i = 0; i < tlim; i++){
        if (rt -> node[(int)s[i]] == NULL)return NULL;
        rt = rt->node[(int)s[i]];
    }
    return rt;
}

bool qcmp(qgramsort a, qgramsort b)
{
    return a.s < b.s;
}

void SimJoiner::clean(trie* rt)
{
    for (int i = 0; i <= 128; i++){
        if (rt -> node[i] != NULL)clean(rt->node[i]);
    }
    delete rt;
}

void SimJoiner::BuildED()
{
    string temp;
    int len;
    qgramsort a[311];
    for (int i = 0; i < sz2; i++) {
        len = data2[i].length();
        for (int j = 0; j < len-qlimit+1; j++){
            a[j+1].id = j;
            a[j+1].s = data2[i].substr(j, qlimit);
        }
        if (len-qlimit+1 > prethresh) {
            sort(a+1, a+1+len-qlimit+1, qcmp);
            for (int j = 1; j <= prethresh; j++) {
                temp = data2[i].substr(a[j].id, qlimit);
                insert(qroot, temp.c_str(), i, qlimit);
            }
        }
        else {
            for (int j = 1; j <= len-qlimit+1; j++) {
                temp = data2[i].substr(a[j].id, qlimit);
                insert(qroot, temp.c_str(), i, qlimit);
            }     
        }
    }
}

int SimJoiner::joinJaccard(const char *filename1, const char *filename2, double threshold, vector<JaccardJoinResult> &result) {
    result.clear();
    return SUCCESS;
}

int SimJoiner::joinED(const char *filename1, const char *filename2, unsigned threshold, vector<EDJoinResult> &result) {
    
    trie* sans;
    string tempt;
    vector<int> cand;

    int len,lt,t,val;
    qgramsort a[311];
    EDJoinResult res;

    result.clear();
    qroot = new trie();

    prethresh = threshold * qlimit + 1;
    createIndex(filename2, data2);
    sz2 = data2.size();
    BuildED();
    createIndex(filename1, data1);
    sz1 = data1.size();

    for (int i = 0; i < sz2; i++)querycheck[i] = 0;
    querytime = 0;

    for (int i = 0; i < sz1; i++) {
        cand.clear();
        querytime++;
        len = data1[i].length();
        for (int j = 0; j < len - qlimit + 1; j++){
            a[j+1].id = j;
            a[j+1].s = data1[i].substr(j, qlimit);
        }
        if (len-qlimit+1 >= prethresh)sort(a+1, a+1+len-qlimit+1, qcmp);
        lt = my_min(prethresh, len-qlimit+1);
        for (int j = 1; j <= lt; j++) {
            tempt = data1[i].substr(a[j].id, qlimit);
            sans = search(qroot, tempt.c_str(), qlimit);
            if (sans != NULL) {
                for(int k = 0; k < sans->qsize; k++){
                    t = sans->qgram[k];
                    if (querycheck[t] != querytime){
                        querycheck[t] = querytime;
                        cand.push_back(t);
                        occurtime[t] = 0;
                    }
                    occurtime[t]++;
                }
            }
        }
        sort(cand.begin(), cand.end());
        for (int j = 0; j < cand.size(); j++){
            t = cand[j];
            qthresh = my_max(data1[i].length(), data2[t].length()) - qlimit + 1 - threshold * qlimit;
            int suf1,suf2;
            suf1 = my_max(data1[i].length() - qlimit + 1 - prethresh, 0);
            suf2 = my_max(data2[t].length() - qlimit + 1 - prethresh, 0);
            if (occurtime[t] + my_min(suf1, suf2) >= qthresh) {
                val = CalCulateED(data1[i].c_str(), data1[i].length(), data2[t].c_str(), data2[t].length(), threshold);
                if (val <= threshold){
                    res.id1 = i;
                    res.id2 = t;
                    res.s = val;
                    result.push_back(res);                    
                }
            }
        }
    }
/*
    for (int i = 0; i < result.size(); i++){
        printf("%d %d %d\n", result[i].id1, result[i].id2, result[i].s);
    }
*/
    clean(qroot);
    return SUCCESS;
}
