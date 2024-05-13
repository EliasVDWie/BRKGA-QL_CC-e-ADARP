//#pragma once
#ifndef _DECODER_H
#define _DECODER_H

#define INFINITO 999999999

#include <math.h>
#include <algorithm>
#include "Data.h"
#include "Scheduler.h"


/************************************************************************************
 Method: Decoder()
 Description: Convert a random key solution in a real problem solution
*************************************************************************************/
TSol Decoder(TSol s, int n, std::vector <TNode> node,
    std::vector <std::vector <double> > dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> > cons);


/************************************************************************************
 Method: Dec1
 Description: standard decoder 
*************************************************************************************/
TSol Dec1(TSol s, int n, std::vector <TNode> node,
    std::vector <std::vector <double> > dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> > cons);
/************************************************************************************
 Method: Dec2
 Description: standard decoder 
*************************************************************************************/
TSol Dec2(TSol s, int n, std::vector <TNode> node,
    std::vector <std::vector <double> > dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> > cons);

/************************************************************************************
 Method: objFct
 Description: calculates the objective function value
*************************************************************************************/
double objFct(TSol s, std::vector <TNode> &node, std::vector <std::vector <double> > &dist, int nbVehicles, int nbUsers, int periodLength, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices);

/************************************************************************************
 Method: AssignFinalDepots
 Description: returns vector with the assigned depot for vehicle 0 at index 0, vehicle 1 at index 1, ...
*************************************************************************************/
std::vector <int> AssignFinalDepots(std::vector<TVecRk> rkVec, int nbUsers, int nbVehicles, std::vector<int> fDepotIDs, std::vector <std::vector <double> > dist);

#endif