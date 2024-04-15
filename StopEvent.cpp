//
// Created by Abdallah Abu Aisha on 6/11/2022.
//

#include "StopEvent.h"

StopEvent::StopEvent() {

}

StopEvent::StopEvent(uint16_t stop_val, int arr_val, int dep_val) {
    stop = stop_val;
    arrival = arr_val;
    departure = dep_val;
}

StopEvent::~StopEvent() {

}