// Example program
#include <bits/stdc++.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

vector<int> findChaines(string line) {
    // Return posistions of &
    vector<int> pos;

    for (int i = 0; i < line.length(); i++) {
        char car = line[i];
        if (car == '&') {
            pos.push_back(i);
        }
    }

    return pos;
}

void sep_chains(string line) {
    /*
    Breaks paired base which are not on the same chain (sparated by '&')
    */

    vector<tuple<int, int>> chains;

    int curChaine = 0, pos, chain;
    for (int i = 0; i < line.length(); i++) {
        char car = line[i];
        if (car == '(') {
            chains.push_back({i, curChaine});
        } else if (car == '&') {
            curChaine++;
        } else if (car == ')') {
            tie(pos, chain) = chains.back();
            chains.pop_back();
            if (chain != curChaine) {
                line[i] = '.';
                line[pos] = '.';
            }
        }
    }

    size_t p = 0;
    vector<string> res;
    while ((p = line.find('&')) != string::npos) {
        res.push_back(line.substr(0, p));
        line.erase(0, p + 1);
    }
    res.push_back(line);

    for (int i = 0; i < res.size(); i++) {
        cout << res[i] << endl;
    }
}

int main() {
    string str = ".((((((((((((....((((((((.....((((((...(.....)...))))..))....)))))).)).((.......((.(((((...))))).)).......))..)))))))))))).&((((((......((((((((((.....(((..(((((((((((((......((((((.....[[[)))...].((((((...(((]].....)))......))))))...)))....(((....)))(((((((......))))))).(((((((((.........))))))).)).))....((((((.....(...[..)......))))).).........((....))...((((......((.]...)).....))))...((((...((((((((....))))))))....{.(((((((((.....((((........(((((((((.(((...(.((((......)))))[[(...).(((.....]])))...)))...))))))))).)))).)))))..)))))....((((((((..[......))))))))....}.(((((].[[[[)))))..)))......)).))))))).))...............((....))....))......(..((((....))))............)).....)))))))))).....(.((((..((((..[...))))..)))).<.(.(........).).((((((..(.(((((((.(((((..(((...(((.....)))...)))...((....))..((((.....))))......))))).)))))))..((.(...((((((.....((((((((((.((((...((((((......))))))...)<)))...(((.(((((.......";
    sep_chains(str);
}
