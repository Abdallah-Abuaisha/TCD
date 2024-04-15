//
// Created by Abdallah Abu Aisha on 13/9/2022.
//

#include "Station.h"

Station::Station() {

}

Station::~Station() {

}

Station::Station(vector<int> depEvents_val, vector<uint16_t> stops_val, vector<int> neighbours_val, int nbhd_val) {
    depEvents = depEvents_val;
    stops = stops_val;
    neighbours = neighbours_val;
    nbhd = nbhd_val;
}
