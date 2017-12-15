
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include "alphanum_comp.h"
#include <unordered_map>
#include "json.hpp"

// for convenience
using json = nlohmann::json;
using namespace std;

int main(int argc, const char *argv[])
{
    // read a JSON file
    std::ifstream i("./files.2017-12-12T19_25_53.669792.json");
    json j;
    i >> j;

    unordered_map<std::string, std::string> id_file;
    // iterate the array
    for (json::iterator it = j.begin(); it != j.end(); ++it)
    {
        json str = (*it)["cases"];
        string ss = (*it)["file_name"];
        id_file[str[0]["case_id"]] = ss.substr(0, ss.length() - 3);
    }

    typedef unordered_map<std::string, std::string> gene_expression;
    typedef unordered_map<std::string, std::string> clinical;
    unordered_map<std::string, gene_expression> id_gene;
    unordered_map<std::string, clinical> id_clin;
    vector<string> patient_uuids;
    vector<string> clin_prop;

    ifstream cutfile("./gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/BRCA.merged_only_clinical_clin_format_cut.txt");
    string line;

    // first line must be patient_uuid
    getline(cutfile, line);
    istringstream iss(line);
    vector<string> tokens;
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));
    string prop = tokens[0];
    clin_prop.push_back(prop);
    clin_prop.push_back("filename");
    for (auto i = tokens.begin() + 1; i < tokens.end(); ++i)
    {
        patient_uuids.push_back(*i);
        clinical clin;
        clin[prop] = *i;
        id_clin[*i] = clin;
    }

    while (getline(cutfile, line))
    {
        istringstream iss(line);
        vector<string> tokens;
        copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

        prop = tokens[0];
        clin_prop.push_back(prop);
        for (size_t i = 0; i < patient_uuids.size(); ++i)
        {
            id_clin[patient_uuids[i]][prop] = tokens[i + 1];
        }
    }

    ifstream normalfile("./normal-list.txt");
    vector<string> normal_list;
    while (getline(normalfile, line))
    {
        normal_list.push_back(line);
    }

    prop = "normal";
    clin_prop.push_back(prop);
    unordered_map<std::string, double> gene_name;
    auto ii = begin(patient_uuids);
    while (ii != end(patient_uuids))
    {
        string i = *ii;
        string exp_file = id_file[i];

        ifstream expfile("./Tumor-txt/" + exp_file);
        if (expfile.is_open())
        {
            id_clin[i][prop] = "0";
        }
        else
        {
            expfile.close();
            expfile.open("./Normal-txt/" + exp_file);
            if (expfile.is_open())
            {
                id_clin[i][prop] = "1";
                normal_list.erase(remove(normal_list.begin(), normal_list.end(), exp_file), normal_list.end());
            }
            else
            {
                cout << exp_file << "\tdoes not exist.\n";
                ii = patient_uuids.erase(ii);
            }
        }
        if (expfile.is_open())
        {
            ++ii;
            id_clin[i]["filename"] = exp_file;

            gene_expression gene_exp;
            double sum_exp = 0.0;
            while (getline(expfile, line))
            {
                istringstream iss(line);
                vector<string> tokens;
                copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

                gene_exp[tokens[0]] = tokens[1];
                sum_exp += stof(tokens[1]);
                gene_name[tokens[0]] = 1;
            }
            gene_exp["sum_exp"] = to_string(sum_exp);
            id_gene[i] = gene_exp;
        }
    }

    string fake_uuid = "fake-uuid-";
    for (size_t i = 0; i < normal_list.size(); ++i)
    {
        string uuid = fake_uuid + to_string(i);
        patient_uuids.push_back(uuid);
        string exp_file = normal_list[i];
        ifstream expfile("./Normal-txt/" + exp_file);
        if (expfile.is_open())
        {
            clinical clin;
            id_clin[uuid] = clin;
            for (auto &ii : clin_prop)
            {
                id_clin[uuid][ii] = "NA";
            }
            id_clin[uuid]["patient.bcr_patient_uuid"] = uuid;
            id_clin[uuid]["filename"] = exp_file;
            id_clin[uuid]["normal"] = "2";

            gene_expression gene_exp;
            double sum_exp = 0.0;
            while (getline(expfile, line))
            {
                istringstream iss(line);
                vector<string> tokens;
                copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

                gene_exp[tokens[0]] = tokens[1];
                sum_exp += stof(tokens[1]);
                gene_name[tokens[0]] = 1;
            }
            gene_exp["sum_exp"] = to_string(sum_exp);
            id_gene[uuid] = gene_exp;
        }
    }

    prop = "trip_negative";
    clin_prop.push_back(prop);
    for (auto &i : patient_uuids)
    {
        int trip = 0;
        if (id_clin[i]["patient.breast_carcinoma_estrogen_receptor_status"] == "negative")
            trip++;
        if (id_clin[i]["patient.breast_carcinoma_progesterone_receptor_status"] == "negative")
            trip++;
        if (id_clin[i]["patient.lab_proc_her2_neu_immunohistochemistry_receptor_status"] == "negative")
            trip++;
        id_clin[i][prop] = to_string(trip);
    }

    sort(patient_uuids.begin(), patient_uuids.end(), [&](string a, string b) {
        if (id_clin[a]["normal"] == id_clin[a]["normal"])
            return id_clin[a]["trip_negative"] < id_clin[b]["trip_negative"];
        return id_clin[a]["normal"] < id_clin[a]["normal"]; });
    // sort(patient_uuids.begin(), patient_uuids.end(), [&](string a, string b) { return id_clin[a]["trip_negative"] < id_clin[b]["trip_negative"]; });

    for (auto &i : gene_name)
    {
        size_t zeros = 0;
        for (auto j : patient_uuids)
        {
            auto t = id_gene[j].find(i.first);
            if (t == id_gene[j].end())
            {
                id_gene[j][i.first] = "0.0";
                zeros++;
            }
            else
            {
                if (id_gene[j][i.first] == "0.0")
                    zeros++;
            }
        }
        i.second = double(zeros) / patient_uuids.size();
    }

    string out_ep, out_trip;

    for (auto &i : clin_prop)
    {
        out_ep += i + '\t';
        out_trip += i + '\t';
        for (auto j : patient_uuids)
        {
            if (id_clin[j]["normal"] == "2")
            {
                out_ep += id_clin[j][i] + '\t';
                out_trip += id_clin[j][i] + '\t';
            }
            else if (id_clin[j]["normal"] == "0")
            {
                if (id_clin[j]["trip_negative"] == "3")
                    out_trip += id_clin[j][i] + '\t';
                else
                    out_ep += id_clin[j][i] + '\t';
            }
        }
        out_ep += '\n';
        out_trip += '\n';
    }

    for (auto &i : gene_name)
    {
        if (i.second < 0.9)
        {
            out_ep += i.first + '\t';
            out_trip += i.first + '\t';
            for (auto j : patient_uuids)
            {
                double norm = stof(id_gene[j][i.first]) / stof(id_gene[j]["sum_exp"]) * 1000000.0;
                string nor = to_string(norm) + '\t';
                // output += id_gene[j][i.first] + '\t';
                if (id_clin[j]["normal"] == "2")
                {
                    out_ep += nor;
                    out_trip += nor;
                }
                else if (id_clin[j]["normal"] == "0")
                {
                    if (id_clin[j]["trip_negative"] == "3")
                        out_trip += nor;
                    else
                        out_ep += nor;
                }
            }
            out_ep += '\n';
            out_trip += '\n';
        }
    }

    ofstream out_epf("./clin_ep_TPM.txt");
    ofstream out_tripf("./clin_trip_TPM.txt");
    out_epf << out_ep;
    out_epf.close();
    out_tripf << out_trip;
    out_tripf.close();
    return 0;
}