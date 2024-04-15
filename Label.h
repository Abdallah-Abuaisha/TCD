//
// Created by Abdallah Abu Aisha on 13/9/2022.
//

#ifndef PTVFTCSA2_STATION_LABEL_H
#define PTVFTCSA2_STATION_LABEL_H

#include <vector>
#include <cstdint>

using namespace std;

class Label {
public:
    int connection;
    uint16_t FT;

    Label();

    ~Label();

    Label(int connection_val, uint16_t FT_val);
};


#endif //PTVFTCSA2_STATION_LABEL_H
