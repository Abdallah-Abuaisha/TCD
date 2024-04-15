//
// Created by Abdallah Abu Aisha on 8/11/2022.
//

#include "Temporary_Label.h"

Temporary_Label::Temporary_Label() {

}

Temporary_Label::Temporary_Label(int con_val, int dep_val, uint16_t FT_val, int EA_val, uint16_t stop_val) {
    con = con_val;
    dep = dep_val;
    FT = FT_val;
    EA = EA_val;
    stop = stop_val;
}

Temporary_Label::~Temporary_Label() {

}