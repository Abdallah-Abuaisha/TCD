//
// Created by Abdallah Abu Aisha on 13/9/2022.
//

#ifndef PTVFTCSA2_TRIP_H
#define PTVFTCSA2_TRIP_H

#include <vector>
#include "StopEvent.h"

using namespace std;

class Trip {
public:
    vector<StopEvent> stopEvents;
    uint16_t route;

    Trip();

    ~Trip();

    Trip(vector<StopEvent> stopEvents_val, uint16_t route_val);
};


#endif //PTVFTCSA2_TRIP_H
