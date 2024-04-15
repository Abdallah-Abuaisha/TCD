//
// Created by Abdallah Abu Aisha on 13/9/2022.
//

#ifndef PTVFTCSA2_STATION_H
#define PTVFTCSA2_STATION_H

#include <string>
#include <vector>
#include <cstdint>

using namespace std;

class Station {
public:
    vector<int> depEvents; // connection indices
    vector<uint16_t> stops;
    vector<int> neighbours;
    int nbhd; // neighbourhood

    Station();

    ~Station();

    Station(vector<int> depEvents_val, vector<uint16_t> stops_val, vector<int> neighbours_val, int nbhd_val);
};


#endif //PTVFTCSA2_STATION_H
