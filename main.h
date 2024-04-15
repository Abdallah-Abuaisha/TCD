//
// Created by Abdallah Abu Aisha on 13/9/2022.
//

#ifndef TCD_HEADER_H
#define TCD_HEADER_H

#include <iostream>
#include <vector>
#include <climits>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <cstdlib>
#include "Station.h"
#include "Connection.h"
#include "Trip.h"
#include "Label.h"
#include "StopEvent.h"
#include "Leg.h"
#include "Transfer.h"
#include "Stop.h"
#include "Service.h"
#include "Temporary_Label.h"

using namespace std;

// functions

bool sort_con(const Connection &con1, const Connection &con2);

bool sort_transfer(const Transfer &trans1, const Transfer &trans2);

bool sort_trip(const Trip &trip1, const Trip &trip2);

int critical_transfer(const vector<Transfer> &vec1, const vector<Transfer> &vec2, const vector<Stop> &stops,
                      const int &source);

vector<int> combineDepartures(const vector<Station> &stations, const vector<vector<int>> &nebs, const int &src);

int findLastConnection(const vector<int> &con_ind, const int &dep);

void createLabels(const int &neighbourhood, const vector<Connection> &connections, const vector<Station> &stations,
                   const vector<Stop> &stops, vector<vector<vector<Label>>> &firstTransferTable,
                   const vector<int> &connectionTimeIndex, const vector<vector<int>> &neighbourhoods);

Leg walking_source_target(const int &source, const int &target, const vector<Station> &stations,
                              const vector<Stop> &stops, const int &startTime);

Leg walking_other_target(const int &target, const int &from_stop, const vector<Station> &stations,
                             const vector<Stop> &stops, const int &curTime);

int findTransferTime(const int &stop1, const int &stop2, const vector<Stop> &stops);

vector<int> min_transfer_time(const int &source, const int &target, const vector<Station> &stations,
                              const vector<Stop> &stops);

int find_label(const vector<Label> &Labels, const int &con_index);

int findArrivalTime(const vector<int> &trip, const int &FTstop, const int &dep, const vector <Connection> &cons);

int findNonOptimalInd(const vector<Leg> &path, const Leg &walkT);

// save labels
void store_labels(const vector<vector<vector<Label>>> &labels, string file_name);

// load labels
vector<vector<vector<Label>>> load_labels(string file_name, size_t num_stations, size_t num_nbhds);

double compute_size(const vector<vector<vector<Label>>> &labels);


// Final GitHub
void buildNetwork (string stopsFile, string stopEventsFile, string transitiveTransfersFile, vector<Station> &stations,
                   vector<Stop> &stops, vector<Connection> &connections, vector<vector<int>> &trips,
                   vector<vector<int>> &neighbourhoods, vector<int> &connectionTimeIndex);
void defineStopsAndStations(string stopsFile, vector<Station> &stations, vector<Stop> &stops);
vector<vector<int>> defineStopEvents(string stopEventsFile);
void checkStopEvents(vector<vector<int>> &stopEvents, vector<Stop> &stops);
void defineConnections(vector<vector<int>> &stopEvents, vector<Connection> &connections);
void defineTrips(vector<vector<int>> &trips, const int &numberOfTrips, const vector<Connection> &connections);
void defineTransfers(string transitiveTransfersFile, vector<Stop> &stops);
void defineNeighbours(vector<Station> &stations, const vector<Stop> &stops);
void defineNeighbourhoods(vector<vector<int>> &neighbourhoods, vector<Station> &stations);
void defineStationDepartures(const vector<Connection> &connections, vector<Station> &stations, const vector<Stop> &stops);
void defineConnectionTimeIndex(const vector<Connection> &connections, vector<int> &connectionTimeIndex);
void buildFirstTransferTable(vector<vector<vector<Label>>> &firstTransferTable, const vector<Connection> &connections,
                             const vector<Station> &stations, const vector<Stop> &stops,
                             const vector<vector<int>> &neighbourhoods, const vector<int> &connectionTimeIndex);
void solveQueries(string queriesFile, const vector<vector<vector<Label>>> &firstTransferTable,
                  const vector<Connection> &connections, const vector<Station> &stations, const vector<Stop> &stops,
                  vector<vector<int>> &trips, vector<int> &connectionTimeIndex);
vector<vector<int>> defineQueries(string queriesFile);
int findPath(const int &origin, const int &destination, const int &departureTime,
              const vector<vector<vector<Label>>> &firstTransferTable, const vector<Connection> &connections,
              const vector<Station> &stations, const vector<Stop> &stops, const vector<vector<int>> &trips,
              const vector<int> &connectionTimeIndex);
void exportSolutions(const vector<int> &solutions);






#endif //TCD_HEADER_H
