//
// Created by Abdallah Abu Aisha on 18/10/2022.
//

#ifndef VBBFTCSAVTT_TRANSFER_H
#define VBBFTCSAVTT_TRANSFER_H

#include <cstdint>

using namespace std;

class Transfer {
public:
    uint16_t destination;
    uint16_t duration;

    Transfer();

    Transfer(uint16_t destination_val, uint16_t duration_val);

    ~Transfer();
};


#endif //VBBFTCSAVTT_TRANSFER_H
