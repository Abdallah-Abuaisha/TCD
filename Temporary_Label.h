    //
    // Created by Abdallah Abu Aisha on 8/11/2022.
    //

    #ifndef VBBFTCSAVTT3_TEMPORARY_LABEL_H
    #define VBBFTCSAVTT3_TEMPORARY_LABEL_H

    #include <cstdint>
    #include <vector>

    using namespace std;

    class Temporary_Label {
    public:
        int con;
        int dep;
        uint16_t FT;
        int EA;
        uint16_t stop;

        Temporary_Label();

        Temporary_Label(int con_val, int dep_val, uint16_t FT_val, int EA_val, uint16_t stop_val);

        ~Temporary_Label();
    };


    #endif //VBBFTCSAVTT3_TEMPORARY_LABEL_H
