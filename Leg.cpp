//
// Created by Abdallah Abu Aisha on 13/9/2022.
//

#include "Leg.h"

Leg::Leg() {

}

Leg::~Leg() {

}

Leg::Leg(int from_val, int to_val, int from_sp_val, int to_sp_val, int departure_val, int arrival_val, int trip_val){
    from = from_val;
    to = to_val;
    from_sp = from_sp_val;
    to_sp = to_sp_val;
    departure = departure_val;
    arrival = arrival_val;
    trip = trip_val;
}