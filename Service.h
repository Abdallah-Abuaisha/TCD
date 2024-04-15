//
// Created by Abdallah Abu Aisha on 8/11/2022.
//

#ifndef VBBFTCSAVTT3_SERVICE_H
#define VBBFTCSAVTT3_SERVICE_H

#include <cstdint>
#include <vector>
#include <string>

using namespace std;

class Service {
public:
    string name;
    vector<uint16_t> stops;

    Service();

    Service(string name_val, vector<uint16_t> stops_val);

    ~Service();
};


#endif //VBBFTCSAVTT3_SERVICE_H
