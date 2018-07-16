#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <getopt.h>
#include <ctime>
#include <cstdlib>
#include <numeric>
#include <iomanip>

using namespace std;

/* Flag set by ‘--verbose’. */
static int verbose_flag;
// random generator function:
int myrandom(int i) { return std::rand() % i; }

void help(int exitcode)
{
    cout << "Usage: genematch OPTION... [FILE]...\n";
    cout << "Print matched chromosome position from each FILE with reference gene symble file to standard output.\n";
    cout << "\n  Options:\n";
    cout << "    -c,\t--chr\t\tchromesome string format like 'chr1'\n";
    cout << "    -f,\t--refile=FILE\trequired, provide reference gene symble file\n";
    cout << "    -m,\t--matched\tprint matched to file\n";
    cout << "    -u,\t--unmatched\tprint unmatched to file\n";
    cout << "    \t--help\t\tdisplay this help and exit\n";
    cout << "    \t--version\toutput version information and exit\n";
    cout << "\nWith no FILE, read standard input.\n\n";
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

int main(int argc, char *argv[])
{
    std::srand(unsigned(std::time(0)));
    string datafile;
    long treated = 0;
    long perm_num = 100;
    static struct option long_options[] =
        {
            /* These options set a flag. */
            {"verbose", no_argument, &verbose_flag, 1},
            {"brief", no_argument, &verbose_flag, 0},
            /* These options don’t set a flag.
             We distinguish them by their indices. */
            {"version", no_argument, 0, 1000},
            {"help", no_argument, 0, 1001},
            {"datafile", required_argument, 0, 'f'},
            {"treated", required_argument, 0, 't'},
            {"perm_num", required_argument, 0, 'p'},
            {0, 0, 0, 0}};
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "p:t:f:", long_options, &option_index)) != -1)
    {
        switch (c)
        {
        case 1000:
            cout << "Permutation Test 0.3\n";
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
            datafile = trim(optarg);
            break;
        case 't':
            treated = stol(optarg);
            break;
        case 'p':
            perm_num = stol(optarg);
            break;
        default:
            // cerr << "Please provide right options.\n\n";
            help(1);
        }
    }

    // where is the ref file ................
    if (datafile == "")
    {
        cerr << "Please provide data file.\n";
        help(1);
    }
    if (treated == 0)
    {
        cerr << "Please provide number of treated.\n";
        help(1);
    }

    ifstream dataf(datafile);
    if (!dataf.is_open())
    {
        cerr << "Can't open " << datafile << "!\n";
        exit(1);
    }

    string line;
    string outerstr;
    while (getline(dataf, line))
    {
        istringstream iss(line);
        vector<double> tokens;
        //copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(stod(tokens));
        double n;
        //string tmp;
        //iss >> tmp;
        while (iss >> n)
            tokens.push_back(n);
        long tokensize = tokens.size();
        double treated_accu = std::accumulate(tokens.begin(), tokens.begin() + treated, 0.0);
        long greater = 0;

        for (long i = 0; i < perm_num; ++i)
        {
            // using myrandom:
            std::random_shuffle(tokens.begin(), tokens.end(), myrandom);
            if (std::accumulate(tokens.begin(), tokens.begin() + treated, 0) > treated_accu)
                greater++;
        }
        if (greater == 0) greater++;
        double pval = ((double)greater) / (perm_num + 1);
        if (pval > 0.5) pval = 1 - pval;
        ostringstream oss;
        oss << std::setprecision (15) << pval << std::endl;
        //string tmstr = to_string(pval) + "\n";
        //cout << tmstr;
        outerstr += oss.str();
    }
    cout << outerstr;
}