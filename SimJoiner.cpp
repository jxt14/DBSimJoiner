#include "SimJoiner.h"

using namespace std;

SimJoiner::SimJoiner() {
    qlimit = 4;
    qgramcount.clear();
}

SimJoiner::~SimJoiner() {
}

bool ccmp(EDJoinResult a, EDJoinResult b)
{
    if(a.id1 != b.id1)return a.id1 < b.id1;
    else return a.id2 < b.id2;
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
	int len,tot;
	ifstream fp;
	p = false;

    datas.clear();
 //   cout << "input " << filename << endl;

	fp.open(filename);
	
	string temp,tq;
    unsigned hq;
    tot = 0;
	while (getline(fp, temp)) {
        tot++;
		len = temp.length();
		datas.push_back(temp);
 //       cout << temp << endl;
	}
	fp.close();

//    cout << "input end" << endl;
}

void SimJoiner::insert(trie* rt, const char* s, int id, int tlim, int kd)
{
	for (int i = 0; i < tlim; i++){
		if (rt -> node[(int)s[i]] == NULL){
			rt->node[(int)s[i]] = new trie();
		}
		rt = rt -> node[(int)s[i]];
	}
    if (kd == 0) {
        if (rt->id == -1) {
            rt->id = ++qtot;
            if (qgramcount.size() < qtot){
                qgramcount.push_back(0);
            }
            else qgramcount[qtot - 1] = 0;
        }
        if (rt->sl != id) {
            rt->sl = id;
            qgramcount[rt->id - 1]++;
        }
    }
    else {
        if (rt->ql != id){
	        rt->ql = id;
            rt->qgram.push_back(id);
	        rt->qsize++;
        }
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

bool SimJoiner::checkmid(const char* smid, const char* s2)
{
    trie *p1,*p2;
    int l1,l2,id1,id2;
    p1 = search(qroot, smid, qlimit);
    if (p1 == NULL || p1->id == -1){
        l1 = 0;
        id1 = -1;
    }
    else{
        l1 = qgramcount[p1->id - 1];
        id1 = p1->id;
    }
    p2 = search(qroot, s2, qlimit);
    if (p2 == NULL || p2->id == -1){
        l2 = 0;
        id2 = -1;
    }
    else{
        l2 = qgramcount[p2->id - 1];
        id2 = p2->id;
    }
    if (l1 != l2) return (l1>l2);
    else return (id1 >= id2);
}

void SimJoiner::clean(trie* rt)
{
    for (int i = 0; i <= 128; i++){
        if (rt -> node[i] != NULL)clean(rt->node[i]);
    }
    delete rt;
}

void SimJoiner::BuildED(int threshold)
{
    string temp;
    int len;
    trie* qsearch;
    pi a[311];
    les.clear();
    lesable.clear();
    for (int i = 0; i < sz2; i++) {
        len = data2[i].length();
        if (len <= threshold * qlimit + qlimit - 1 + threshold){
            lesable.push_back(i);
        }
        if (len - qlimit + 1 - threshold * qlimit <= 0) {
            les.push_back(i);
        }
        else{
            for (int j = 0; j < len-qlimit+1; j++){
                a[j+1].second = j;
                temp = data2[i].substr(j, qlimit);
                qsearch = search(qroot, temp.c_str(), qlimit);
                if (qsearch == NULL || qsearch->id == -1){
                    a[j+1].first.first = 0;
                    a[j+1].first.second = -1;
                }
                else{
                    a[j+1].first.first = qgramcount[qsearch->id - 1];
                    a[j+1].first.second = qsearch->id;
                }
            }
            sort(a+1, a+1+len-qlimit+1);
            for (int j = 1; j <= len-qlimit+1; j++)qsrt[i][j] = a[j].second;
            for (int j = 1; j <= prethresh; j++) {
                temp = data2[i].substr(a[j].second, qlimit);
                insert(qroot, temp.c_str(), i, qlimit, 1);
            }
        }
    }
//    cout << "Build end " << endl;
}

bool SimJoiner::sufcheck(int id1, int* qs1, int id2, int thresh)
{
    int L1,R1,L2,R2,mc,mid,ans;
    string sc,st;
    if (thresh <= 0)return true;
    L1 = prethresh + 1;
    R1 = data1[id1].length() - qlimit + 1;
    mc = (L1+R1)/2;
    sc = data1[id1].substr(qs1[mc], qlimit);
    while (mc > L1) {
        if(sc != data1[id1].substr(qs1[mc-1], qlimit))break;
        mc--;
    }
    L2 = prethresh + 1;
    R2 = data2[id2].length() - qlimit + 1;
    while (L2<=R2) {
        mid = (L2+R2)/2;
        st = data2[id2].substr(qsrt[id2][mid], qlimit);
        if(checkmid(st.c_str(), sc.c_str())) R2 = mid-1;
        else L2 = mid+1;
    }
    ans = L2;
    L2 = prethresh + 1;
    R2 = data2[id2].length() - qlimit + 1;
    st = data2[id2].substr(qsrt[id2][ans], qlimit);
    int numcheck;
    if (st == sc){
        numcheck = 1 + my_max(mc - L1, ans - L2) + my_min(R1 - mc, R2 - ans);
    }
    else {
        numcheck = my_max(mc - L1, ans - L2) + my_min(R1 - mc, R2 - ans + 1);
    }
    if (numcheck >= thresh)return true;
    else return false;
}

int SimJoiner::joinJaccard(const char *filename1, const char *filename2, double threshold, vector<JaccardJoinResult> &result) {
    result.clear();
    return SUCCESS;
}

int SimJoiner::joinED(const char *filename1, const char *filename2, unsigned threshold, vector<EDJoinResult> &result) {
    
    trie* sans;
    trie* qsearch;
    string tempt;
    unsigned hq;
    int cand[200011];
    int candtot,suf1,suf2;
    int* tempq;
    qtot = 0;

    int Len,t,val;
    pi a[311];
    EDJoinResult res;

    result.clear();
    qroot = new trie();

    prethresh = threshold * qlimit + 1;
    createIndex(filename2, data2);
    sz2 = data2.size();
    for (int i = 0; i < sz2; i++){
        tempq = new int[data2[i].length()+1];
        if(qsrt.size() <= i)qsrt.push_back(tempq);
        else qsrt[i] = tempq;
    }
    for (int i = 0; i < sz2; i++){
        Len = data2[i].length();
        for (int j = 0; j < Len - qlimit + 1; j++){
            tempt = data2[i].substr(j, qlimit);
            insert(qroot, tempt.c_str(), i, qlimit, 0);
        }
    }

    BuildED(threshold);
    createIndex(filename1, data1);

    sz1 = data1.size();

    for (int i = 0; i < sz2; i++)querycheck[i] = 0;
    querytime = 0;

    for (int i = 0; i < sz1; i++) {
        candtot = 0;
        querytime++;
        Len = data1[i].length();
        if (Len - qlimit + 1 - threshold * qlimit <= 0) {
            for (int j = 0; j < lesable.size(); j++){
                t = lesable[j];
                val = CalCulateED(data1[i].c_str(), data1[i].length(), data2[t].c_str(), data2[t].length(), threshold);
                if (val <= threshold){
                    res.id1 = i;
                    res.id2 = t;
                    res.s = val;
                    result.push_back(res);
                }
            }
        }
        else {
            for (int j = 0; j < Len - qlimit + 1; j++){
                a[j+1].second = j;
                tempt = data1[i].substr(j, qlimit);
                qsearch = search(qroot, tempt.c_str(), qlimit);
                if (qsearch == NULL || qsearch->id == -1){
                    a[j+1].first.first = 0;
                    a[j+1].first.second = -1;
                }
                else{
                    a[j+1].first.first = qgramcount[qsearch->id - 1];
                    a[j+1].first.second = qsearch->id;
                }
            }
            sort(a+1, a+1+Len-qlimit+1);
            tempq = new int[Len-qlimit+1+1];
            for (int j = 1; j <= Len-qlimit+1; j++)tempq[j] = a[j].second;
            for (int j = 1; j <= prethresh; j++) {
                tempt = data1[i].substr(a[j].second, qlimit);
                sans = search(qroot, tempt.c_str(), qlimit);
                if (sans != NULL) {
                    for(int k = 0; k < sans->qsize; k++){
                        t = sans->qgram[k];
                        if (querycheck[t] != querytime){
                            querycheck[t] = querytime;
                            cand[++candtot] = t;
                            occurtime[t] = 0;
                        }
                        occurtime[t]++;
                    }
                }
            }
            if (Len <= threshold * qlimit + qlimit - 1 + threshold){
                for (int j = 0; j < les.size(); j++){
                    t = les[j];
                    val = CalCulateED(data1[i].c_str(), data1[i].length(), data2[t].c_str(), data2[t].length(), threshold);
                    if (val <= threshold){
                        res.id1 = i;
                        res.id2 = t;
                        res.s = val;
                        result.push_back(res);
                    }
                }
            }
            for (int j = 1; j <= candtot; j++){
                t = cand[j];
                qthresh = my_max(data1[i].length(), data2[t].length()) - qlimit + 1 - threshold * qlimit;
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
            delete tempq;
        }
    }

    sort(result.begin(), result.end(), ccmp);
/*
    for (int i = 0; i < result.size(); i++){
        printf("%d %d %d\n", result[i].id1, result[i].id2, result[i].s);
    }
    printf("\n");
*/
    clean(qroot);
    for (int i = 0; i < sz2; i++)delete qsrt[i];

    return SUCCESS;
}
