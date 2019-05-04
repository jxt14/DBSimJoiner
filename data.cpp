#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<cstring>
#include<string>
#include<iostream>
#include<fstream>

using namespace std;

int main()
{
    int lim;
    scanf("%d", &lim);

    freopen("sp.txt", "w", stdout);
    string s;
    char c;
    int k;
    srand(unsigned(time(NULL)));
    for(int i = 1; i <= lim; i++){
        for (int j = 1; j <= 256; j++){
            k = rand()%26;
            c = 'a' + k;
            if (j == 1) s = c;
            else s = s + c; 
        }
        cout << s << endl;
    }
    ofstream fp;
    fp.open("readin.txt");
    for(int i = 1; i <= lim; i++){
        for (int j = 1; j <= 256; j++){
            k = rand()%6;
            c = 'a' + k;
            if (j == 1) s = c;
            else s = s + c; 
        }
        fp << s << endl;
    }
    fp.close();
}