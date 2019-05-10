#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<cstring>
#include<string>
#include<set>
#include<utility>
#include<vector>
#include<algorithm>

using namespace std;

int main()
{
    vector<pair<int,string>> sv;
    sv.clear();
    sv.push_back(make_pair(1,"gaa"));
    sv.push_back(make_pair(1,"aaa"));
    sv.push_back(make_pair(2,"aa"));
    sort(sv.begin(), sv.end());
    for (int i = 0; i < sv.size(); i++){
        cout << sv[i].first << " " << sv[i].second << endl;
    }
}