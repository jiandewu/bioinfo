#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include "alphanum_comp.h"
using namespace std;

int main(int argc, const char *argv[])
{
    struct Gene_ref
    {
        string gene_id;
        string chromosome;
        size_t start_bp;
        size_t end_bp;
        string gene_name;
        Gene_ref(string s1, string s2, size_t s3, size_t s4, string s5) : gene_id(s1), chromosome(s2), start_bp(s3), end_bp(s4), gene_name(s5) {}
    };

    // heterogeneous comparison:
    struct Comp0
    {
        bool operator()(const Gene_ref &s1, const Gene_ref &s2) const
        {
            int i = doj::alphanum_comp(s1.chromosome, s2.chromosome);
            if (i != 0)
                return i < 0 ? true : false;
            return s1.start_bp < s2.start_bp;
        }
    };
    struct Comp1
    {
        bool operator()(const Gene_ref & s, string str) const { return doj::alphanum_comp(s.chromosome, str) < 0 ? true : false; }
        bool operator()(string str, const Gene_ref & s) const { return doj::alphanum_comp(str, s.chromosome) < 0 ? true : false; }
    };
    struct Comp2
    {
        bool operator()(const Gene_ref & s1, size_t s2) const { return s1.start_bp < s2; }
        bool operator()(size_t s2, const Gene_ref & s1) const { return s2 < s1.start_bp; }
    };

    string inputfile = argv[1];
    string outputfile = inputfile + "_gene++.txt";

    ifstream reffile("Gene_location.txt");
    string line;
    getline(reffile, line);

    vector<Gene_ref> vec_gene;
    while (getline(reffile, line))
    {
        istringstream iss(line);
        vector<string> tokens;
        copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));
        // Chromosome name format: chr1
        Gene_ref gene(tokens[0], "chr" + tokens[1], stoi(tokens[2]), stoi(tokens[3]), tokens[4]);
        // Chromosome name format: 1
        // Gene_ref gene(tokens[0], tokens[1], stoi(tokens[2]), stoi(tokens[3]), tokens[4]);
        vec_gene.push_back(gene);
    }

    sort(vec_gene.begin(), vec_gene.end(), Comp0());

    ifstream sourcefile(inputfile);
    ofstream out(outputfile);

    string output;
    while (getline(sourcefile, line))
    {
        istringstream iss(line);
        vector<string> tokens;
        copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

        if (tokens[3] == "B")
        {
            size_t position = stoi(tokens[1]);
            auto p = equal_range(vec_gene.begin(), vec_gene.end(), tokens[0], Comp1());
            auto up = upper_bound(p.first, p.second, position, Comp2());

            for (auto low = p.first; low != up; ++low) {
                if (low->end_bp >= position) {
                   output += (low->gene_name + '\t' + line + '\n');
                }
            }

            // // may faster methord
            // vector<size_t> vecinx(up - p.first);
            // iota(vecinx.begin(), vecinx.end(), p.first - vec_gene.begin());
            // sort(vecinx.begin(), vecinx.end(), [&vec_gene](size_t i1, size_t i2) {return vec_gene[i1].end_bp < vec_gene[i2].end_bp;});
            // auto upp = upper_bound(vecinx.begin(), vecinx.end(), position, [&vec_gene](size_t i2, size_t i1) {return i2 < vec_gene[i1].end_bp;});
            // for (auto low = upp; low != vecinx.end(); ++low) {
            //     output += (vec_gene[*low].gene_name + '\t' + line + '\n');
            // }
        }
    }

    out << output;
    out.close();
    return 0;
}