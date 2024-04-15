//
// Created by Abdallah Abu Aisha on 13/9/2022.
//

#include "Trip.h"

Trip::Trip() {

}

Trip::~Trip() {

}

Trip::Trip(vector<StopEvent> stopEvents_val, uint16_t route_val) {
    stopEvents = stopEvents_val;
    route = route_val;
}