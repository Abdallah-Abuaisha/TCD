//
// Created by Abdallah Abu Aisha on 6/11/2022.
//

#ifndef VBBFTCSAVTT3_STOPEVENT_H
#define VBBFTCSAVTT3_STOPEVENT_H

#include <cstdint>

class StopEvent {
public:
    uint16_t stop;
    int arrival;
    int departure;

    StopEvent();

    StopEvent(uint16_t stop_val, int arr_val, int dep_val);

    ~StopEvent();
};


#endif //VBBFTCSAVTT3_STOPEVENT_H
