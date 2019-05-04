#include "SimJoiner.h"

using namespace std;

SimJoiner::SimJoiner() {
    qlimit = 4;
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

void SimJoiner::createIndex(const char *filename, vector<string>& datas, map<string,int>& qgrams)
{
	char c;
	bool p;
	int len,tot;
	ifstream fp;
	p = false;
    map<string, int> qcheck;

    qcheck.clear();
    datas.clear();
    qgrams.clear();
 //   cout << "input " << filename << endl;

	fp.open(filename);
	
	string temp,tq;
    tot = 0;
	while (getline(fp, temp)) {
        tot++;
		len = temp.length();
		datas.push_back(temp);
        for (int i = 0; i < len-qlimit+1; i++){
            tq = temp.substr(i, qlimit);
            if (qcheck[tq] != tot){
                qcheck[tq] = tot;
                qgrams[tq]++;
            }
        }
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

void SimJoiner::clean(trie* rt)
{
    for (int i = 0; i <= 128; i++){
        if (rt -> node[i] != NULL)clean(rt->node[i]);
    }
    delete rt;
}

void SimJoiner::BuildED(unsigned threshold)
{
    string temp;
    int len;
    pi a[311];
    les2.clear();
    for (int i = 0; i < sz2; i++) {
        len = data2[i].length();
        if (len - qlimit + 1 - threshold * qlimit <= 0) {
            les2.push_back(i);
        }
        else{
            for (int j = 0; j < len-qlimit+1; j++){
                a[j+1].second = j;
                a[j+1].first = qgram2[data2[i].substr(j, qlimit)];
            }
            sort(a+1, a+1+len-qlimit+1);
            for (int j = 1; j <= prethresh; j++) {
                temp = data2[i].substr(a[j].second, qlimit);
                insert(qroot, temp.c_str(), i, qlimit);
            }
        }
    }
//    cout << "Build end " << endl;
}

int SimJoiner::joinJaccard(const char *filename1, const char *filename2, double threshold, vector<JaccardJoinResult> &result) {
    result.clear();
    return SUCCESS;
}

int SimJoiner::checkED(int id1, int id2, unsigned threshold)
{
    int val = CalCulateED(data1[id1].c_str(), data1[id1].length(), data2[id2].c_str(), data2[id2].length(), threshold);
    return val;
}

int SimJoiner::joinED(const char *filename1, const char *filename2, unsigned threshold, vector<EDJoinResult> &result) {
    
    trie* sans;
    string tempt;
    vector<int> cand;

    int len,lt,t,val;
    pi a[311];
    EDJoinResult res;

    result.clear();
    qroot = new trie();

    prethresh = threshold * qlimit + 1;
    createIndex(filename2, data2, qgram2);
    int tt = 2;
 //  cout << "input end " << endl;
    sz2 = data2.size();
    BuildED(threshold);
    createIndex(filename1, data1, qgram1);
    sz1 = data1.size();

    for (int i = 0; i < sz2; i++)querycheck[i] = 0;
    querytime = 0;

    for (int i = 0; i < sz1; i++) {
        cand.clear();
        querytime++;
        len = data1[i].length();
        if (len - qlimit + 1 - threshold * qlimit <= 0) {
            for (int j = 0; j < sz2; j++)cand.push_back(j);
        }
        else {
            for (int j = 0; j < len - qlimit + 1; j++){
                a[j+1].second = j;
                a[j+1].first = qgram1[data1[i].substr(j, qlimit)];
            }
            sort(a+1, a+1+len-qlimit+1);
            for (int j = 1; j <= prethresh; j++) {
                tempt = data1[i].substr(a[j].second, qlimit);
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
            for (int j = 0; j < les2.size(); j++)cand.push_back(les2[j]);
        }
        sort(cand.begin(), cand.end());
 //       printf("cand %d ", i);
        for (int j = 0; j < cand.size(); j++){
            t = cand[j];
            qthresh = my_max(data1[i].length(), data2[t].length()) - qlimit + 1 - threshold * qlimit;
            int suf1,suf2;
            suf1 = my_max(data1[i].length() - qlimit + 1 - prethresh, 0);
            suf2 = my_max(data2[t].length() - qlimit + 1 - prethresh, 0);
            if (occurtime[t] + my_min(suf1, suf2) >= qthresh) {
                val = checkED(i, t, threshold);
                if (val <= threshold){
                    res.id1 = i;
                    res.id2 = t;
                    res.s = val;
                    result.push_back(res);
                }
            }
        }
 //       printf("\n");
    }
/*
    for (int i = 0; i < result.size(); i++){
        printf("%d %d %d\n", result[i].id1, result[i].id2, result[i].s);
    }
*/
 //   clean(qroot);
    return SUCCESS;
}
