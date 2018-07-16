#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include "alphanum_comp.h"
#include <getopt.h>
using namespace std;

/* Flag set by ‘--verbose’. */
static int verbose_flag;

void help(int exitcode)
{
    cout << "Usage: genematch OPTION... [FILE]...\n";
    cout << "Print matched chromosome position from each FILE with reference gene symble file to standard output.\n";
    cout << "\n  Options:\n";
    cout << "    -f,\t--reffile=FILE\trequired, provide reference gene symble file\n";
    cout << "    \t\t\tgene symble file format (column seprator is TAB):\n";
    cout << "    \t\t\tGene stable ID\tChromosome\tGene start\tGene end\tGene name\n";
    cout << "    \t\t\tENSG00000283891\t15\t\t55372940\t55373034\tMIR628\n";
    cout << "    -c,\t--chr\t\tchromesome string format: 'chr1'\n";
    cout << "    -m,\t--matched\tprint matched to file\n";
    cout << "    -u,\t--unmatched\tprint unmatched to file\n";
    cout << "    \t--help\t\tdisplay this help and exit\n";
    cout << "    \t--version\toutput version information and exit\n";
    cout << "\nWith no FILE, read standard input.\n";
    cout << "FILE format: first 2 columns should be (column seprator is TAB) chr & position\n";
    cout << "chr17\t51391264\nOR...\n";
    cout << "17\t51391264\n\n";
    cout << "Examples:\n";
    cout << "gene_match -cmuf Gene_location.txt mutation.txt\n";
    cout << "gene_match -cf Gene_location.txt mutation.txt > out_matched.txt\n\n";
    cout << "Report bugs to: jwu2@lsuhsc.edu\n";
    cout << "pkg home page: <https://www.medschool.lsuhsc.edu/bioinformatics/>\n";
    cout << "General help using GNU software: <http://www.gnu.org/gethelp/>\n\n";
    exit(exitcode);
}

std::string trim(const std::string &str,
                 const std::string &whitespace = " \t")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

std::string reduce(const std::string &str,
                   const std::string &fill = " ",
                   const std::string &whitespace = " \t")
{
    // trim first
    auto result = trim(str, whitespace);

    // replace sub ranges
    auto beginSpace = result.find_first_of(whitespace);
    while (beginSpace != std::string::npos)
    {
        const auto endSpace = result.find_first_not_of(whitespace, beginSpace);
        const auto range = endSpace - beginSpace;

        result.replace(beginSpace, range, fill);

        const auto newStart = beginSpace + fill.length();
        beginSpace = result.find_first_of(whitespace, newStart);
    }

    return result;
}

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
    bool operator()(const Gene_ref &s, string str) const { return doj::alphanum_comp(s.chromosome, str) < 0 ? true : false; }
    bool operator()(string str, const Gene_ref &s) const { return doj::alphanum_comp(str, s.chromosome) < 0 ? true : false; }
};
struct Comp2
{
    bool operator()(const Gene_ref &s1, size_t s2) const { return s1.start_bp < s2; }
    bool operator()(size_t s2, const Gene_ref &s1) const { return s2 < s1.start_bp; }
};

void match(istream &sourcefile, const vector<Gene_ref> &vec_gene, string &output, string &outnomatch)
{
    string line;
    // size_t linedebug = 1;
    while (getline(sourcefile, line))
    {
        // line = trim(line);
        if (line.at(0) == '#')
            continue;

        istringstream iss(line);
        vector<string> tokens;
        copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

        bool outline = false;
        // if (tokens[3] == "B")
        if (true)
        {
            try
            {
                size_t position = stoi(tokens[1]);
                auto p = equal_range(vec_gene.begin(), vec_gene.end(), tokens[0], Comp1());
                auto up = upper_bound(p.first, p.second, position, Comp2());

                for (auto low = p.first; low != up; ++low)
                {
                    if (low->end_bp >= position)
                    {
                        output += (low->gene_name + '\t' + line + '\n');
                        outline = true;
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
            catch (...)
            {
                // String could not be read properly as an int.
            }
        }
        if (!outline)
        {
            outnomatch += line + '\n';
        }
    }

    // string outputfile = inputfile + "_gene++.txt";
    // ofstream out(outputfile);
    // out << output;
    // out.close();

    // outputfile = inputfile + "_gene++_nomatch.txt";
    // ofstream out2(outputfile);
    // out2 << outnomatch;
    // out2.close();
}

void outfile(string matchedstr, string nonmatchedstr, bool outmatched_flag, bool outunmatched_flag, string inputname)
{
    if (outmatched_flag)
    {
        string outputfile = inputname + "_gene++.txt";
        ofstream out(outputfile);
        out << matchedstr;
        out.close();
    }
    else
    {
        cout << matchedstr;
    }
    if (outunmatched_flag)
    {
        string outputfile = inputname + "_gene++_unmatch.txt";
        ofstream out(outputfile);
        out << nonmatchedstr;
        out.close();
    }
}

int main(int argc, char *argv[])
{
    string reffilename;
    vector<string> inputfiles;
    bool chrformat_flag = false;
    bool outmatched_flag = false;
    bool outunmatched_flag = false;

    static struct option long_options[] =
        {
            /* These options set a flag. */
            {"verbose", no_argument, &verbose_flag, 1},
            {"brief", no_argument, &verbose_flag, 0},
            /* These options don’t set a flag.
             We distinguish them by their indices. */
            {"version", no_argument, 0, 1000},
            {"help", no_argument, 0, 1001},
            {"matched", no_argument, 0, 'm'},
            {"unmatched", no_argument, 0, 'u'},
            {"chrformat", no_argument, 0, 'c'},
            {"refile", required_argument, 0, 'f'},
            {0, 0, 0, 0}};
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "mucf:", long_options, &option_index)) != -1)
    {
        switch (c)
        {
        case 1000:
            cout << "Gene Match 0.32\n";
            cout << "Copyright (C) 2017 LSUHSC Bioinformatics and Genomics (BIG) Program.\n";
            cout << "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n";
            cout << "This is free software: you are free to change and redistribute it.\n";
            cout << "There is NO WARRANTY, to the extent permitted by law.\n";
            exit(0);
            break;
        case 1001:
            help(0);
            break;
        case 'f':
            reffilename = trim(optarg);
            break;
        case 'c':
            chrformat_flag = true;
            break;
        case 'm':
            outmatched_flag = true;
            break;
        case 'u':
            outunmatched_flag = true;
            break;
        default:
            // cerr << "Please provide right options.\n\n";
            help(1);
        }
    }
    /* Print any remaining command line arguments (not options). */
    while (optind < argc)
        inputfiles.push_back(trim(argv[optind++]));

    // where is the ref file ................
    if (reffilename == "")
    {
        cerr << "Please provide reference gene symbol file.\n";
        help(1);
    }

    ifstream reffile(reffilename);
    if (!reffile.is_open())
    {
        cerr << "Can't open " << reffilename << "!\n";
        exit(1);
    }

    string line;
    getline(reffile, line);

    vector<Gene_ref> vec_gene;
    while (getline(reffile, line))
    {
        istringstream iss(line);
        vector<string> tokens;
        copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));
        // For chromosome name format: 1
        Gene_ref gene(tokens[0], tokens[1], stoi(tokens[2]), stoi(tokens[3]), tokens[4]);
        // For chromosome name format: chr1
        if (chrformat_flag)
            gene.chromosome = "chr" + gene.chromosome;
        vec_gene.push_back(gene);
    }

    sort(vec_gene.begin(), vec_gene.end(), Comp0());

    string matchedstr;
    string nonmatchedstr;
    if (inputfiles.empty())
    {
        match(cin, vec_gene, matchedstr, nonmatchedstr);
        outfile(matchedstr, nonmatchedstr, outmatched_flag, outunmatched_flag, "gene_match");
    }
    else
    {
        for (auto inputfile : inputfiles)
        {
            matchedstr = "";
            nonmatchedstr = "";
            ifstream sourcefile(inputfile);
            if (!sourcefile.is_open())
            {
                cerr << "Can't open " << inputfile << "!\n";
                continue;
            }

            match(sourcefile, vec_gene, matchedstr, nonmatchedstr);
            outfile(matchedstr, nonmatchedstr, outmatched_flag, outunmatched_flag, inputfile);
        }
    }

    return 0;
}