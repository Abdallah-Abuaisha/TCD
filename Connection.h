//
// Created by Abdallah Abu Aisha on 13/9/2022.
//

#ifndef PTVFTCSA2_CONNECTION_H
#define PTVFTCSA2_CONNECTION_H

#include <cstdint>

using namespace std;

class Connection {
public:
    uint16_t from;
    uint16_t to;
    int departure;
    int arrival;
    int trip;

    Connection();

    Connection(uint16_t from_val, uint16_t to_val, int departure_val, int arrival_val, int trip_val);

    ~Connection();
};


#endif //PTVFTCSA2_CONNECTION_H
