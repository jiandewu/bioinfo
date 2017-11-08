#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include "alphanum_comp.h"
using namespace std;

int main(int argc, const char * argv[])
{
    // heterogeneous comparison:
    struct Comp0
    {
        bool operator() (const vector<string>& s1, const vector<string>& s2 ) const {
            int i = doj::alphanum_comp(s1[1], s2[1]);
            if (i != 0) return i < 0 ? true : false;
            return stoi(s1[2]) < stoi(s2[2]);
        }
    };
    struct Comp1
    {
        bool operator() ( const vector<string>& s, string str ) const { return doj::alphanum_comp(s[1], str) < 0 ? true : false; }
        bool operator() ( string str, const vector<string>& s ) const { return doj::alphanum_comp(str, s[1]) < 0 ? true : false; }
    };
    struct Comp2
    {
        bool operator() ( const vector<string>& s, string str ) const { return stoi(s[2]) < stoi(str); }
        bool operator() ( string str, const vector<string>& s ) const { return stoi(str) < stoi(s[2]); }
    };

    string inputfile = argv[1];
    string outputfile = inputfile + "_gene2++.txt";

    ifstream reffile("Gene_location.txt");
    string line;
    vector<vector<string> > refgene;
    while (getline(reffile, line)) {
        istringstream iss(line);
        vector<string> tokens;
        copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));
        refgene.push_back(tokens);
    }
    
    sort(refgene.begin(),refgene.end(),Comp0());

    ifstream sourcefile(inputfile);
    ofstream out(outputfile);

    string output;
    while (getline(sourcefile, line)) {
        istringstream iss(line);
        vector<string> tokens;
        copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));
        
        if (tokens[3] == "B") {
            auto p = equal_range(refgene.begin(), refgene.end(), tokens[0], Comp1());
            auto up = upper_bound(p.first, p.second, tokens[1], Comp2());

            for (auto low = p.first; low != up; ++low) {
                if (stoi((*low)[3]) >= stoi(tokens[1])) {
                   output += ((*low)[4] + '\t' + line + '\n');
                }
            }

            // // may faster methord
            // vector<size_t> vecinx(up - p.first);
            // iota(vecinx.begin(), vecinx.end(), p.first - refgene.begin());
            // sort(vecinx.begin(), vecinx.end(), [&refgene](size_t i1, size_t i2) {return stoi(refgene[i1][3]) < stoi(refgene[i2][3]);});
            // auto upp = upper_bound(vecinx.begin(), vecinx.end(), tokens[1], [&refgene](string i2, size_t i1) {return stoi(i2) < stoi(refgene[i1][3]);});
            // for (auto low = upp; low != vecinx.end(); ++low) {
            //     output += (refgene[*low][4] + '\t' + line + '\n');
            // }
        }
    }

    out << output;
    out.close();
    return 0;
}