/*
 * Copyright 2012-2013 TopCoder, Inc.

 *
 * This code was developed under U.S. government contract NNH10CD71C. 

 *
 * Licensed under the Apache License, Version 2.0 (the "License");

 * You may not use this file except in compliance with the License.

 * You may obtain a copy of the License at:
 *     http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software

 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

 * See the License for the specific language governing permissions and

 * limitations under the License.
 */

#include <vector>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "place1.h"
#include "place2.h"
#include "place3.h"
#include "place5.h"

using namespace std;

string inttostr(int x) {
    ostringstream oss;
    oss << x;
    return oss.str();
}

void error(string s) {
    cout << "ERROR: " << s << endl;
    exit(0);
}

vector<string> read_data(string file_name) {
    ifstream fin(file_name.c_str(), ifstream::in);

    if (fin.fail()) {
        error("Unable to read data from \"" + file_name + "\".");
    }

    string s;
    getline(fin, s);
    istringstream iss(s);
    int n;
    iss >> n;

    vector<string> res;
    for (int i = 0; i < n; i++) {
        getline(fin, s);
        res.push_back(s);
    }

    return res;
}

const int MIN_DAY = 5440;
const int MIN_REC_DAY = 11284;
const int MAX_REC_DAY = 17644;
const int MAX_DAY = 17858;

int main(int argc, char *argv[]) {
    if (argc != 7) {
        error("This program is to be executed as follows:\nmain <input folder> <output folder> <solution> <first training day> <first testing day> <last testing day>");
    }

    string data_folder = argv[1];
    string out_folder = argv[2];

    vector<string> events_data = read_data(data_folder + "/events.txt");
    vector<string> regions_data = read_data(data_folder + "/regions.txt");

    vector<vector<string>> events(MAX_DAY + 1);

    for (int i = 0; i < events_data.size(); i++) {
        istringstream iss(events_data[i]);
        int day, country, region;
        string lat, lon;
        iss >> day >> lat >> lon >> country >> region;
        ostringstream oss;
        oss << lat << " " << lon << " " << country << " " << region;
        events[day].push_back(oss.str());
    }

    int start_train, start_test, end_test;
    istringstream(argv[4]) >> start_train;
    istringstream(argv[5]) >> start_test;
    istringstream(argv[6]) >> end_test;

    if (start_train < MIN_DAY) {
        error("The value of <first training day> parameter must be " + inttostr(MIN_DAY) + " or above.");
    }

    if (start_train < MIN_REC_DAY) {
        cout << "WARNING: The value of <first training day> parameter is recommended to be " << MIN_REC_DAY << " or above." << endl;
    }

    if (start_test < start_train) {
        error("The value of <first testing day> parameter must be greater than or equal to the value of <first training day> parameter.");
    }

    if (end_test < start_test) {
        error("The value of <last testing day> parameter must be greater than or equal to the value of <first testing day> parameter.");
    }

    if (end_test > MAX_REC_DAY) {
        cout << "WARNING: The value of <last testing day> parameter is recommended to be " << MAX_REC_DAY << " or below." << endl;
    } 

    if (end_test > MAX_DAY) {
        error("The value of <last testing day> parameter must be " + inttostr(MAX_DAY) + " or below.");
    }

    string solution = argv[3];
    BaseMassAtrocityPredictor* obj;
    if (solution == "1") {
        obj = new place1::MassAtrocityPredictor();
    } else if (solution == "2") {
        obj = new place2::MassAtrocityPredictor();
    } else if (solution == "3") {
        obj = new place3::MassAtrocityPredictor();
    } else if (solution == "5") {
        obj = new place5::MassAtrocityPredictor();
    } else {
        error("The value of <solution> parameter must be one of \"1\", \"2\", \"3\", or \"5\" (quotes for clarity).");
        return 0;
    }

    obj->receiveData(2, start_train, regions_data);

    for (int cur_day = start_train; cur_day <= end_test; cur_day++) {
        cout << "Day = " << cur_day << endl;

        vector<string> data = read_data(data_folder + "/data_" + inttostr(cur_day) + ".txt");

        obj->receiveData(1, cur_day, data);
        obj->receiveData(0, cur_day, events[cur_day]);

        if (cur_day >= start_test) {
            vector<double> res = obj->predictAtrocities(cur_day);
            string res_file_name = out_folder + "/res_" + inttostr(cur_day) + ".txt";
            ofstream ofs(res_file_name, ofstream::out);
            if (ofs.fail()) {
                error("Unable to save data to \"" + res_file_name + "\".");
            }
            for (int i = 0; i < res.size(); i++) {
                ofs << res[i] << endl;
            }
            ofs.close();
        }
    }

    return 0;
}