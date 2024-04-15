//
// Created by Abdallah Abu Aisha on 19/10/2022.
//

#include "Stop.h"

Stop::Stop() {

}

Stop::Stop(uint16_t parent_val, vector<Transfer> transfers_val, vector<Transfer> minTransfers_val){
    parentStation = parent_val;
    transfers = transfers_val;
    minTransfers = minTransfers_val; // ONLY FOR JOURNEY REPLAN
}

Stop::~Stop() {

}