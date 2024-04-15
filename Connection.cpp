//
// Created by Abdallah Abu Aisha on 13/9/2022.
//

#include "Connection.h"

Connection::Connection() {

}

Connection::~Connection() {

}

Connection::Connection(uint16_t from_val, uint16_t to_val, int departure_val, int arrival_val, int trip_val){
    from = from_val;
    to = to_val;
    departure = departure_val;
    arrival = arrival_val;
    trip = trip_val;
}