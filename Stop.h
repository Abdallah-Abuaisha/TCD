//
// Created by Abdallah Abu Aisha on 19/10/2022.
//

#ifndef VBBFTCSAVTT_PARENT_STATION_H
#define VBBFTCSAVTT_PARENT_STATION_H

#include <cstdint>
#include <vector>
#include "Transfer.h"

using namespace std;

class Stop {
public:
    uint16_t parentStation;
    vector<Transfer> transfers;

    vector<Transfer> minTransfers; // journey replan only!


    Stop();

    Stop(uint16_t parent_val, vector<Transfer> transfers_val, vector<Transfer> minTransfers_val);

    ~Stop();
};


#endif //VBBFTCSAVTT_PARENT_STATION_H
