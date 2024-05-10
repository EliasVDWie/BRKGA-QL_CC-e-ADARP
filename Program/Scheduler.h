#pragma once
#include <math.h>
#include <algorithm>
#include "Data.h"

/************************************************************************************
 Method: Scheduler()
 Description: Takes a routing solution as argument and returns the solution with time variables (charging decisions will be included in the method later)
*************************************************************************************/
TSol Scheduler(TSol s, std::vector <TNode>& node, std::vector <std::vector <double> > &dist, std::vector <TVhcl> vehicle, int nbUsers,
    int nbVehicles, int periodLength, std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> > &cons);

/************************************************************************************
 Method: findClosestCharger
 Description: Finds closest charger for which the intial time window of the zeroload is completely available
*************************************************************************************/
TVecSol findClosestCharger(TVecSol node1, TVecSol node2, std::vector <TCstat> cStations, std::vector <std::vector <double> >& dist, std::vector <TNode>& node, std::vector<std::vector<std::vector<unsigned long long>>>& availability, int& stationIndex);

/************************************************************************************
 Method: checkAvailability
 Description: Checks whether chargign station with index 'csIndex' is available for the full period between start and end
 Returns: true if available, false if at max capacity during at least 1 minute of the interval
*************************************************************************************/
bool checkAvailability(float start, float end, int csIndex, std::vector <TCstat> cStations, std::vector<std::vector<std::vector<unsigned long long>>>& availability);

/************************************************************************************
 Method: updateAvailability
 Description: Updates a charging session at charging station with index csIndex from 'start' to 'end'
*************************************************************************************/
void updateAvailability(float start, float end, int csIndex, std::vector <TCstat> cStations, std::vector<std::vector<std::vector<unsigned long long>>>& availability);

/************************************************************************************
Method: ScheduleLatePUEarlyDO
Description: Sets T variables for all routes in a solution, as late as possible for pickups, based on price for charging stations, as early as possible for dropoffs
*************************************************************************************/
void ScheduleLatePUEarlyDO(TSol& s, int nbVehicles, int nbUsers, std::vector <TNode>& node, std::vector <std::vector <double> >& dist, std::vector <double> elecPrices, int periodLength);

/************************************************************************************
Method: ScheduleRKPUEarlyDO
Description: Sets T variables for all routes in a solution, RK-based for pickups, price-based for charging stations, as early as possible for dropoffs.
Note: Pickups directly after a charging stationed are scheduled as late as possible to ensure consistency with LT of charging station (calculated earlier in algorithm)
*************************************************************************************/
void ScheduleRkPUEarlyDO(TSol& s, int nbVehicles, int nbUsers, std::vector <TNode>& node, std::vector <std::vector <double> >& dist, std::vector <double> elecPrices, int periodLength);


struct ZL { // Struct used for representing a zero-load point
    int node;
    double price;
    double window; //window of opportunity during which charging could possibly occur (not taking into account travel times to charging station etc. initially) // LT next pick-up - ET zero-load drop-off - drop-off service time
    double rk;
};


