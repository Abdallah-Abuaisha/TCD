// 07/04/24
// Clean baseline algorithm for TCD to be uploaded to GitHub


#include "main.h"


// *********************************************************************************************************************

int main() {
    // define and build public transport network
    vector<Station> allStations;
    vector<Stop> allStops;
    vector<Connection> allConnections;
    vector<vector<int>> allTrips, allNeighbourhoods;
    vector<int> connectionTimeIndex;

    buildNetwork("stops", "stopEvents", "transitiveTransfers", allStations, allStops, allConnections, allTrips,
                 allNeighbourhoods, connectionTimeIndex);

    vector<vector<vector<Label>>> firstTransferTable(allNeighbourhoods.size(), vector<vector<Label>>(allStations.size()));
    buildFirstTransferTable(firstTransferTable, allConnections, allStations, allStops, allNeighbourhoods, connectionTimeIndex);

    solveQueries("queries", firstTransferTable, allConnections, allStations, allStops, allTrips, connectionTimeIndex);







//     storing labels on disk
//    store_labels(all_labels, "labels_min");

//    cout << "Total Preprocessing space(MB): " << compute_size(all_labels) << endl;

// *********************************************************************************************************************
//     load precomputed labels from disk
//    auto start0 = chrono::steady_clock::now();
//
//    vector<vector<vector<Label>>> all_labels = load_labels("labels", all_stations.size(), all_nbhds.size());
//
//    auto end0 = chrono::steady_clock::now();
//    auto elapsed0 = chrono::duration_cast<chrono::seconds>(end0 - start0).count();
//    cout << "elapsed time (sec): " << elapsed0 << endl;
//
//    // Manage and reallocate unused memory
//    all_labels.shrink_to_fit();
//    for (int O{0}; O < all_labels.size(); O++) {
//        all_labels[O].shrink_to_fit();
//        for (int D{0}; D < all_labels.size(); D++) {
//            all_labels[O][D].shrink_to_fit();
//        }
//    }
//
//    cout << "Total Preprocessing space(MB): " << compute_size(all_labels) << endl;



//     run each query at the different 8 times

//     5 runs in a row

//    vector<vector<int>> query_times(5);
//    vector<vector<int>> vehicle_legs(5); // to record number of legs for each query at each hour
//    vector<vector<int>> lab_touch(5); // to record number of evaluated labels

//    vector<vector<int>> query_orig(5);
//    vector<vector<int>> query_dest(5);

    vector<int> EA;


//    for(int i{0}; i < lab_touch.size(); i++){
//        for (int q{0}; q < queries[0].size(); q++) {

//            auto start = chrono::high_resolution_clock::now();
//            vector<Segment> path = find_path_wait(queries[0][q], queries[1][q], queries[2][q], all_labels, all_trips, con_ind,
//                             all_connections, all_stops, all_stations, all_routes, labels_touched);
//            pureCSAmod(queries[0][q], queries[1][q], queries[2][q], all_connections, all_stations, all_stops, con_ind);
//            auto end = chrono::high_resolution_clock::now();
//            auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - start).count();

//            query_orig[i].push_back(path[0].from_sp);
//            query_dest[i].push_back(path.back().to_sp);


//            query_times[i].push_back(elapsed);
//            vehicle_legs[i].push_back(count_segments(path));
//            lab_touch[i].push_back(labels_touched);

//            EA.push_back(path.back().arrival);

//            cout << "query " << q << " is completed!" << endl;

//        }
//    }

//    // Exporting results to a text file
//
//    ofstream final_output3{"../EAminWait.txt"};
//    for (int q{0}; q < queries[0].size(); q++) {
//        final_output3 << EA[q] << endl;
//    }
//    final_output3.close();


//    ofstream final_output1{"../orig_stops.txt"};
//    for (int q{0}; q < query_orig[0].size(); q++) {
//        final_output1 << query_orig[0][q] << " " << query_orig[1][q] << " " << query_orig[2][q] << " " <<
//        query_orig[3][q] << " " << query_orig[4][q] << " " << endl;
//    }
//    final_output1.close();


// *********************************************************************************************************************


    return 0;
}

// *********************************************************************************************************************
// *********************************************************************************************************************
// *********************************************************************************************************************

// functions definitions

void buildNetwork (string stopsFile, string stopEventsFile, string transitiveTransfersFile, vector<Station> &stations,
                   vector<Stop> &stops, vector<Connection> &connections, vector<vector<int>> &trips,
                   vector<vector<int>> &neighbourhoods, vector<int> &connectionTimeIndex){

    defineStopsAndStations(stopsFile, stations, stops);
    vector<vector<int>> stopEvents = defineStopEvents(stopEventsFile);
    checkStopEvents(stopEvents, stops);
    defineConnections(stopEvents, connections);
    int numberOfTrips = connections.back().trip + 1;
    sort(connections.begin(), connections.end(), sort_con);
    defineTrips(trips, numberOfTrips, connections);
    defineTransfers(transitiveTransfersFile, stops);
    defineNeighbours(stations, stops);
    defineNeighbourhoods(neighbourhoods, stations);
    defineStationDepartures(connections, stations, stops);
    defineConnectionTimeIndex(connections, connectionTimeIndex);

    cout << "Network has been built!" << endl;
}



// read stops file, define stops and stations, and assign stops and stations
void defineStopsAndStations(string stopsFile, vector<Station> &stations, vector<Stop> &stops){
    vector<vector<int>> stopsAndStations(2);

    ifstream in_file;
    string line;
    int stopID, stationID;
    in_file.open("../" + stopsFile + ".txt");
    if (!in_file) {
        cerr << "Problem Opening File" << endl;
        // return 1;
    }
    while (in_file >> stopID >> stationID) {
        stopsAndStations[0].push_back(stopID);
        stopsAndStations[1].push_back(stationID);
    }
    in_file.close();

    Station s0;
    s0.stops.push_back(stopsAndStations[0][0]);
    stations.push_back(s0);
    Stop p0;
    p0.parentStation = stopsAndStations[1][0];
    stops.push_back(p0);
    for(int i{1}; i < stopsAndStations[0].size(); i++){
        Stop p;
        p.parentStation = stopsAndStations[1][i];
        stops.push_back(p);
        if (stopsAndStations[1][i] == stopsAndStations[1][i-1]){
            stations.back().stops.push_back(stopsAndStations[0][i]);
        } else {
            Station s;
            s.stops.push_back(stopsAndStations[0][i]);
            stations.push_back(s);
        }
    }
}

// read stop events file and define stop events
vector<vector<int>> defineStopEvents(string stopEventsFile){
    vector<vector<int>> stopEvents(4);

    ifstream in_file;
    int stopID, sequence, arrivalTime, departureTime;
    in_file.open("../" + stopEventsFile + ".txt");
    if (!in_file) {
        cerr << "Problem Opening File" << endl;
        //return 1;
    }
    while (in_file >> stopID >> sequence >> arrivalTime >> departureTime) {
        stopEvents[0].push_back(stopID);
        stopEvents[1].push_back(sequence);
        stopEvents[2].push_back(arrivalTime);
        stopEvents[3].push_back(departureTime);
    }
    in_file.close();

    return stopEvents;
}

// check stop events and modify inconsistent events
void checkStopEvents(vector<vector<int>> &stopEvents, vector<Stop> &stops){
    // delete all bad connections having same departure and arrival stop (or parent station)!
    for (int i{0}; i < stopEvents[0].size() - 1; i++) {
        if ((stops[stopEvents[0][i]].parentStation == stops[stopEvents[0][i + 1]].parentStation ||
                stopEvents[0][i] == stopEvents[0][i + 1]) && stopEvents[1][i] < stopEvents[1][i + 1]) {
            stopEvents[2][i + 1] = stopEvents[2][i];
            stopEvents[0].erase(stopEvents[0].begin() + i);
            stopEvents[1].erase(stopEvents[1].begin() + i);
            stopEvents[2].erase(stopEvents[2].begin() + i);
            stopEvents[3].erase(stopEvents[3].begin() + i);
            i--;
        }
    }

    // fix the issue where two or more stops on the same trip have the same exact departure and arrival times
    for (int i{1}; i < stopEvents[0].size(); i++) {
        if (stopEvents[1][i] > stopEvents[1][i - 1] && stopEvents[2][i] <= stopEvents[2][i - 1]) {
            int diff = stopEvents[2][i - 1] - stopEvents[2][i];
            stopEvents[2][i] += diff + 1;
            stopEvents[3][i] += diff + 1;
        }
    }
}

// define connections from stop events
void defineConnections(vector<vector<int>> &stopEvents, vector<Connection> &connections){
    int trip_number {0};
    for (int i{0}; i < stopEvents[0].size() - 1; i++){
        if (stopEvents[1][i] < stopEvents[1][i + 1]){
            Connection c;
            c.from = stopEvents[0][i];
            c.to = stopEvents[0][i+1];
            c.departure = stopEvents[3][i];
            c.arrival = stopEvents[2][i+1];
            c.trip = trip_number;
            connections.push_back(c);
        } else {
            trip_number++;
        }
    }
}

// define trips as sequence of connections
void defineTrips(vector<vector<int>> &trips, const int &numberOfTrips, const vector<Connection> &connections){
    trips.resize(numberOfTrips);
    for (int c{0}; c < connections.size(); c++){
        trips[connections[c].trip].push_back(c);
    }
}

// define transitive transfers and add them to respective stops
void defineTransfers(string transitiveTransfersFile, vector<Stop> &stops){
    vector<vector<int>> transfers(3);
    ifstream in_file;
    int from, to, duration;

    in_file.open("../" + transitiveTransfersFile + ".txt");
    if (!in_file) {
        cerr << "Problem Opening File" << endl;
        //return 1;
    }
    while(in_file >> from >> to >> duration){
        transfers[0].push_back(from);
        transfers[1].push_back(to);
        transfers[2].push_back(duration);
    }
    in_file.close();

    for (int i{0}; i < transfers[0].size(); i++) {
        Transfer tran;
        tran.destination = transfers[1][i];
        tran.duration = transfers[2][i];
        stops[transfers[0][i]].transfers.push_back(tran);
    }

    // sort all transfers in all stops by their destination (OPTIONAL)
    for (int s{0}; s < stops.size(); s++) {
        sort(stops[s].transfers.begin(), stops[s].transfers.end(), sort_transfer);
    }
}

// define neighbouring stations of each station, if any, based on existence of footpaths between their stops
void defineNeighbours(vector<Station> &stations, const vector<Stop> &stops){
    for (int s{0}; s < stations.size(); s++) {
        for (int t{0}; t < stops[stations[s].stops[0]].transfers.size(); t++) {
            if (stops[stops[stations[s].stops[0]].transfers[t].destination].parentStation != s &&
                count(stations[s].neighbours.begin(), stations[s].neighbours.end(),
                      stops[stops[stations[s].stops[0]].transfers[t].destination].parentStation) == 0) {
                stations[s].neighbours.push_back(
                        stops[stops[stations[s].stops[0]].transfers[t].destination].parentStation);
            }
        }
    }
}

// define neighbourhood stations based on footpaths/transfers between neighbouring stations
void defineNeighbourhoods(vector<vector<int>> &neighbourhoods, vector<Station> &stations){
    vector<bool> checked(stations.size(), false);
    for (int s{0}; s < stations.size(); s++) {
        if (stations[s].neighbours.size() == 0) {
            neighbourhoods.push_back({s});
            checked[s] = true;
        } else {
            if (!checked[s]) {
                neighbourhoods.push_back(stations[s].neighbours);
                neighbourhoods.back().push_back(s);
                for (int n{0}; n < stations[s].neighbours.size(); n++) {
                    checked[stations[s].neighbours[n]] = true;
                }
            }
        }
    }

    // add the respective neighbourhood for each station
    for (int N{0}; N < neighbourhoods.size(); N++) {
        for (int s{0}; s < neighbourhoods[N].size(); s++) {
            stations[neighbourhoods[N][s]].nbhd = N;
        }
    }
}

// define departures for each station by passing the connection index
void defineStationDepartures(const vector<Connection> &connections, vector<Station> &stations, const vector<Stop> &stops){
    for (int c{0}; c < connections.size(); c++) {
        stations[stops.at(connections[c].from).parentStation].depEvents.push_back(c);
    }

    // order connections by departure time but in reverse order (latest to earliest)! as required by TCD preprocessing
    for (int s{0}; s < stations.size(); s++) {
        reverse(stations[s].depEvents.begin(), stations[s].depEvents.end());
    }
}

// create connection-time index to store the index of the first connection departing at each second/minute in the day
void defineConnectionTimeIndex(const vector<Connection> &connections, vector<int> &connectionTimeIndex){
    connectionTimeIndex.resize(connections.back().departure + 1, -1);
    for (int i{0}; i < connections.size(); i++) {
        if (i == 0 || connections[i].departure != connections[i - 1].departure){
            connectionTimeIndex[connections[i].departure] = i;
        }
    }

    for (int i{static_cast<int>(connectionTimeIndex.size() - 2)}; i >= 0; i--) {
        if (connectionTimeIndex[i] == -1) {
            connectionTimeIndex[i] = connectionTimeIndex[i + 1];
        }
    }
}

void buildFirstTransferTable(vector<vector<vector<Label>>> &firstTransferTable, const vector<Connection> &connections,
                             const vector<Station> &stations, const vector<Stop> &stops,
                             const vector<vector<int>> &neighbourhoods, const vector<int> &connectionTimeIndex){
    for (int i{0}; i < neighbourhoods.size(); i++) {
        createLabels(i, connections, stations, stops, firstTransferTable, connectionTimeIndex, neighbourhoods);
        cout << i + 1 << " neighbourhoods completed out of " << neighbourhoods.size() << " " <<
             int(double(i + 1) / neighbourhoods.size() * 100) << "%" << endl;
    }

    // Manage and reallocate unused memory
    firstTransferTable.shrink_to_fit();
    for (int O{0}; O < firstTransferTable.size(); O++) {
        firstTransferTable[O].shrink_to_fit();
        for (int D{0}; D < firstTransferTable.size(); D++) {
            firstTransferTable[O][D].shrink_to_fit();
        }
    }

    cout << "Oracle has been built!" << endl;
}

// define queries then solve them using the first transfer table
void solveQueries(string queriesFile, const vector<vector<vector<Label>>> &firstTransferTable,
                  const vector<Connection> &connections, const vector<Station> &stations, const vector<Stop> &stops,
                  vector<vector<int>> &trips, vector<int> &connectionTimeIndex){
    vector<vector<int>> queries(3);
    queries = defineQueries(queriesFile);

    vector<int> allSolutions;
    int solution{-1};
    for (int q{0}; q < queries[0].size(); q++){
        solution = findPath(queries[0][q], queries[1][q], queries[2][q], firstTransferTable, connections, stations, stops, trips,
                  connectionTimeIndex);
        allSolutions.push_back(solution);
    }

    exportSolutions(allSolutions);

    cout << "Queries have been solved and saved!" << endl;
}

// read queries from queries file
vector<vector<int>> defineQueries(string queriesFile){
    vector<vector<int>> queries(3);

    ifstream in_file;
    int origin, destination, departureTime;

    in_file.open("../" + queriesFile + ".txt");
    if (!in_file) {
        cerr << "Problem Opening File" << endl;
        //return 1;
    }
    while(in_file >> origin >> destination >> departureTime){
        queries[0].push_back(origin);
        queries[1].push_back(destination);
        queries[2].push_back(departureTime);
    }
    in_file.close();

    return queries;
}

// extract optimal path/journey from first transfer table that solves a given query (origin station, destination station,
// earliest departure time)
int findPath(const int &origin, const int &destination, const int &departureTime,
              const vector<vector<vector<Label>>> &firstTransferTable, const vector<Connection> &connections,
              const vector<Station> &stations, const vector<Stop> &stops, const vector<vector<int>> &trips,
              const vector<int> &connectionTimeIndex){
    // initialisation
    vector<Leg> path;
    Leg walkToDestination(INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX);
    int curCon, labelIndex, arrival, stop1, stop2, FTStop, FTStation;
    int curTime{departureTime}, curFTStation{origin}, curNbhd{stations[origin].nbhd}; // cur = current
    bool pruneFlag{false};

    // main loop
    while (curFTStation != destination && !pruneFlag) {

        // check if walking to target is possible and update walk_target if needed
        if (count(stations[curFTStation].neighbours.begin(), stations[curFTStation].neighbours.end(), destination)) {
            if (curFTStation == origin) {
                walkToDestination = walking_source_target(curFTStation, destination, stations, stops, curTime);
            } else {
                Leg walkOtoT = walking_other_target(destination, path.back().to_sp, stations, stops, curTime);
                if (walkOtoT.arrival < walkToDestination.arrival) {
                    walkToDestination = walkOtoT;
                }
            }
        }
        // find the candidate label
        labelIndex = find_label(firstTransferTable[curNbhd][destination], connectionTimeIndex[curTime]); // changed to consider L[N][t]
        bool is_reachable{false};
        while (!is_reachable) {
            if (labelIndex == -1) { // && curFTStation == source DELETED!
                if (walkToDestination.arrival !=
                    INT_MAX) { // if directly walking from source to target is always the best/only way to get to target
                    path.push_back(walkToDestination);
                } else { // if no way to get from source to target by any mean including walking
                    Leg error(INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX);
                    path.push_back(error);
                }
//                return path;
                return path.back().arrival;
            }
            curCon = firstTransferTable[curNbhd][destination][labelIndex].connection; // changed to consider L[N][t]
            FTStop = firstTransferTable[curNbhd][destination][labelIndex].FT; // changed to consider L[N][t]
            FTStation = stops[FTStop].parentStation;
            // find the respective transfer time and stop to start from
            int transTime{0};
            stop2 = connections[curCon].from;
            if (curFTStation != origin) {
                stop1 = path.back().to_sp;
                transTime = findTransferTime(stop1, stop2, stops);
            } else if (curFTStation == origin && stops[stop2].parentStation != curFTStation) {
                vector<int> temp = min_transfer_time(curFTStation, stop2, stations, stops);
                transTime = temp[0];
                stop1 = temp[1];
            }
            // check if the current connection is reachable
            if (curTime + transTime <= connections[curCon].departure) {
                is_reachable = true;
                if (stops[stop2].parentStation != curFTStation) { // if connection starts from a neighbour station
                    Leg walk(curFTStation, stops[stop2].parentStation, stop1, stop2, curTime, curTime + transTime,
                                 -1);
                    path.push_back(walk);
                }
                // create PT leg
                arrival = findArrivalTime(trips[connections[curCon].trip], FTStop,
                                        connections[curCon].departure, connections);
                Leg seg(stops[stop2].parentStation, FTStation, stop2, FTStop, connections[curCon].departure,
                            arrival, connections[curCon].trip);
                path.push_back(seg);
                // prepare for next leg
                curTime = seg.arrival;
                curFTStation = FTStation;
                curNbhd = stations[FTStation].nbhd;
                // check if walking to target is better than PT (last leg)
                if (walkToDestination.arrival < seg.arrival) {
                    pruneFlag = true;
                }
            } else { // if the current connection is NOT rechable
                labelIndex++;
                if (labelIndex == firstTransferTable[curNbhd][destination].size()) {
                    pruneFlag = true;
                    is_reachable = true;
                }
            }
        }
    }
    // check whether the last leg should be walking to target
    if (pruneFlag) {
        int non_optimal_index = findNonOptimalInd(path, walkToDestination);
        if (non_optimal_index != INT_MAX) {
            path.erase(path.begin() + (non_optimal_index + 1), path.end());
            path.push_back(walkToDestination);
        } else {
            path.clear();
            path.push_back(walkToDestination);
        }
    }
//    return path;
    return path.back().arrival;
}


void exportSolutions(const vector<int> &solutions){
    ofstream out_file{"../solutions.txt"};

    for (int q{0}; q < solutions.size(); q++) {
        out_file << solutions[q] << endl;
    }
    out_file.close();
};









bool sort_con(const Connection &con1, const Connection &con2) {
    if (con1.departure < con2.departure)
        return true;
    else if (con1.departure == con2.departure && con1.arrival < con2.arrival)
        return true;
    else
        return false;
}



bool sort_transfer(const Transfer &trans1, const Transfer &trans2) {
    if (trans1.destination < trans2.destination)
        return true;
    else
        return false;
}


bool sort_trip(const Trip &trip1, const Trip &trip2) {
    if (trip1.stopEvents[0].departure < trip2.stopEvents[0].departure)
        return true;
    else
        return false;
}

// combine all departure events departing from all stations within a given neighbourhood
vector<int> combineDepartures(const vector<Station> &stations, const vector<vector<int>> &nbhds, const int &src) {
    vector<int> all_departures = stations[nbhds[src][0]].depEvents;
    for (int s{1}; s < nbhds[src].size(); s++) {
        for (int de{0}; de < stations[nbhds[src][s]].depEvents.size(); de++) {
            all_departures.push_back(stations[nbhds[src][s]].depEvents[de]);
        }
    }
    sort(all_departures.begin(), all_departures.end(), greater{});
    return all_departures;
}

// find the last connection to scan in the connections array based on maximum travel time (e.g., 4 hours)
int findLastConnection(const vector<int> &con_ind, const int &dep) {
    int cn{INT_MAX};
    if (dep + 14400 <= con_ind.size()) { // 4 hr - 240 min - 14400 sec
        cn = con_ind[dep + 14400];
    } else {
        cn = con_ind.back();
    }
    return cn;
}

// Preprocessing and building oracle (label sets in OD matrix)
void createLabels(const int &neighbourhood, const vector<Connection> &connections, const vector<Station> &stations,
                   const vector<Stop> &stops, vector<vector<vector<Label>>> &firstTransferTable,
                   const vector<int> &connectionTimeIndex, const vector<vector<int>> &neighbourhoods) {
    // initialisation
    vector<vector<Temporary_Label>> tempo(stations.size());
    vector<int> all_departures = combineDepartures(stations, neighbourhoods, neighbourhood);

    for (int d{0}; d < all_departures.size(); d++) {
        int source = stops[connections[all_departures[d]].from].parentStation;

        int c0 = all_departures[d]; // first connection to start from
        int source_stop = connections[c0].from;
        int dep = connections[c0].departure;
        int trip = connections[c0].trip;
        int cn = findLastConnection(connectionTimeIndex, dep); // last connection (240min)

        vector<int> EA(stops.size(), INT_MAX);
        vector<int> FT_stop(stops.size(), INT_MAX);

        // scan connections
        for (int i{c0}; i < cn; i++) {
            if (EA[connections[i].from] <= connections[i].departure || i == c0) {
                if (connections[i].arrival < EA[connections[i].to]) {
                    // determine first transfer station, if any
                    if (connections[i].trip == trip) {
                        FT_stop[connections[i].to] = connections[i].to;
                    } else {
                        FT_stop[connections[i].to] = FT_stop[connections[i].from];
                    }
                    // relax footpaths
                    for (int j{0}; j < stops[connections[i].to].transfers.size(); j++) {
                        if (connections[i].arrival + stops[connections[i].to].transfers[j].duration <
                            EA[stops[connections[i].to].transfers[j].destination]) {
                            EA[stops[connections[i].to].transfers[j].destination] =
                                    connections[i].arrival + stops[connections[i].to].transfers[j].duration;
                            FT_stop[stops[connections[i].to].transfers[j].destination] = FT_stop[connections[i].to];
                        }
                    }
                }
            }
        }

        // applying dominance checks
        for (int sta{0}; sta < stations.size(); sta++) {
            if (sta != source && sta != stops[source_stop].parentStation) {
                int min_EA = EA[stations[sta].stops[0]];
                int ind = stations[sta].stops[0];
                for (int s{1}; s < stations[sta].stops.size(); s++) {
                    if (EA[stations[sta].stops[s]] < min_EA) {
                        min_EA = EA[stations[sta].stops[s]];
                        ind = stations[sta].stops[s];
                    }
                }

                if (min_EA != INT_MAX && (tempo[sta].empty() || min_EA < tempo[sta].back().EA)) {
                    int FT = FT_stop[ind];
                    Temporary_Label tl(c0, dep, FT, min_EA, source_stop);
                    tempo[sta].push_back(tl);
                } else if (min_EA < INT_MAX && min_EA >= tempo[sta].back().EA && tempo[sta].size() > 0) {
                    if (source_stop != tempo[sta].back().stop) {
                        int index = tempo[sta].size() - 1;
                        bool flag{true};
                        while (index != -1 && min_EA >= tempo[sta][index].EA) {
                            int max_diff = findTransferTime(source_stop, tempo[sta][index].stop, stops);
                            if (tempo[sta][index].dep - dep >= max_diff) {
                                flag = false;
                                break;
                            }
                            index--;
                        }
                        index++;
                        if (flag) {
                            int FT = FT_stop[ind];
                            Temporary_Label tl(c0, dep, FT, min_EA, source_stop);
                            tempo[sta].insert(tempo[sta].begin() + index, tl);
                        }
                    }
                }
            }
        }
        EA.clear();
        EA.shrink_to_fit();
        FT_stop.clear();
        FT_stop.shrink_to_fit();
    }

    // creating final labels
    for (int s{0}; s < stations.size(); s++) {
        if (!tempo[s].empty()) {
            reverse(tempo[s].begin(), tempo[s].end());
            for (int i{0}; i < tempo[s].size(); i++) {
                Label l(tempo[s][i].con, tempo[s][i].FT);
                firstTransferTable[neighbourhood][s].push_back(l);
            }
        }
    }
    tempo.clear();
    tempo.shrink_to_fit();
    all_departures.clear();
    all_departures.shrink_to_fit();
}


int count_segments(const vector<Leg> &path) {
    int number{0};
    for (auto seg: path) {
        if (seg.trip != -1) {
            number++;
        }
    }
    return number;
}

Leg walking_source_target(const int &source, const int &target, const vector<Station> &stations,
                              const vector<Stop> &stops, const int &startTime) {
    vector<int> min(3, INT_MAX);
    for (int sp{0}; sp < stations[source].stops.size(); sp++) {
        for (int tr{0}; tr < stops[stations[source].stops[sp]].transfers.size(); tr++) {
            if (stops[stops[stations[source].stops[sp]].transfers[tr].destination].parentStation == target &&
                stops[stations[source].stops[sp]].transfers[tr].duration < min[0]) {
                min[0] = stops[stations[source].stops[sp]].transfers[tr].duration;
                min[1] = stations[source].stops[sp];
                min[2] = stops[stations[source].stops[sp]].transfers[tr].destination;
            }
        }
    }
    Leg walkStoT(source, target, min[1], min[2], startTime, startTime + min[0], -1);
    return walkStoT;
}

Leg walking_other_target(const int &target, const int &from_stop, const vector<Station> &stations,
                             const vector<Stop> &stops, const int &curTime) {
    vector<int> min_sp(2, INT_MAX);
    for (int sp{0}; sp < stations[target].stops.size(); sp++) {
        for (int tr{0}; tr < stops[stations[target].stops[sp]].transfers.size(); tr++) {
            if (stops[stations[target].stops[sp]].transfers[tr].destination == from_stop) {
                if (stops[stations[target].stops[sp]].transfers[tr].duration < min_sp[0]) {
                    min_sp[0] = stops[stations[target].stops[sp]].transfers[tr].duration;
                    min_sp[1] = stations[target].stops[sp];
                }
                break;
            }
        }
    }
    Leg walkOtoT(stops[from_stop].parentStation, target, from_stop, min_sp[1], curTime, curTime + min_sp[0], -1);
    return walkOtoT;
}

vector<int>
min_transfer_time(const int &source, const int &target, const vector<Station> &stations, const vector<Stop> &stops) {
    vector<int> min_sp(2, INT_MAX);
    for (int sp{0}; sp < stations[source].stops.size(); sp++) {
        for (int tr{0}; tr < stops[stations[source].stops[sp]].transfers.size(); tr++) {
            if (stops[stations[source].stops[sp]].transfers[tr].destination == target) {
                if (stops[stations[source].stops[sp]].transfers[tr].duration < min_sp[0]) {
                    min_sp[0] = stops[stations[source].stops[sp]].transfers[tr].duration;
                    min_sp[1] = stations[source].stops[sp];
                }
                break;
            }
        }
    }
    return min_sp;
}

// find the duration/cost of the transfer between two stops
int findTransferTime(const int &stop1, const int &stop2, const vector<Stop> &stops) {
    for (int i{0}; i < stops[stop1].transfers.size(); i++) {
        if (stops[stop1].transfers[i].destination == stop2) {
            return stops[stop1].transfers[i].duration;
        }
    }
    return INT_MAX;
}

int find_label(const vector<Label> &Labels, const int &con_index) {
    for (int i{0}; i < Labels.size(); i++) {
        if (Labels[i].connection >= con_index) {
            return i;
        }
    }
    return -1;
}

// modified version when trips are defined as sequences of connections (not stop events)
int findArrivalTime(const vector<int> &trip, const int &FTstop, const int &dep, const vector <Connection> &cons) {
    for (int i{0}; i < trip.size(); i++) {
        if (cons[trip[i]].to == FTstop && cons[trip[i]].arrival >= dep) { // modified to solve the "back-in-time issue"
            return cons[trip[i]].arrival;
        }
    }
    return -1;
}

int findNonOptimalInd(const vector<Leg> &path, const Leg &walkT) {
    for (int i{0}; i < path.size(); i++) {
        if (path[i].to_sp == walkT.from_sp) {
            return i;
        }
    }
    return INT_MAX;
}

void store_labels(const vector<vector<vector<Label>>> &labels, string file_name) {
    ofstream output{"../" + file_name + ".txt"};
//    if (!output) {
//        cerr << "Error creating file" << endl;
//    }
    for (int O{0}; O < labels.size(); O++) {
        for (int D{0}; D < labels[O].size(); D++) {
            for (int L{0}; L < labels[O][D].size(); L++) {
                output << O << " " << D << " " << labels[O][D][L].connection << " " << labels[O][D][L].FT << '\n';
            }
        }
    }
//    output.close();
}

vector<vector<vector<Label>>> load_labels(string file_name, size_t num_stations, size_t num_nbhds) {
    ifstream input{"../" + file_name + ".txt"};
//    if (!input) {
//        cerr << "Error creating file" << endl;
//    }

    bool flag{true};
    vector<vector<vector<Label>>> labels(num_nbhds, vector<vector<Label>>(num_stations));
    int con;
    int O, D, FT;
    string line;
    while (!input.eof()) {
        getline(input, line);
        input >> O;
        input >> D;
        input >> con;
        input >> FT;
//        if(FT != INT_MAX)
        labels[O][D].push_back(Label(con, FT));
        if (O >= num_nbhds / 2 && flag) {
            cout << "50% of labels are loaded!" << endl;
            flag = false;
        }
    }
    return labels;
}

double compute_size(const vector<vector<vector<Label>>> &labels) {
    double size{0};
    for (int s{0}; s < labels.size(); s++) {
        for (int t{0}; t < labels[s].size(); t++) {
            size += labels[s][t].size() * 2 * 3; // 1 short integer and 1 int with 2 and 4 bytes
        }
    }
    return (size / 1000000);
}