//
// Created by Abdallah Abu Aisha on 13/9/2022.
//

#ifndef TCD_LEG_H
#define TCD_LEG_H

#include <cstdint>

class Leg {
public:
    int from;
    int to;
    int from_sp;
    int to_sp;
    int departure;
    int arrival;
    int trip;


    Leg();

    ~Leg();

    Leg(int from_val, int to_val, int from_sp_val, int to_sp_val, int departure_val, int arrival_val, int trip_val);
};


#endif //TCD_LEG_H
