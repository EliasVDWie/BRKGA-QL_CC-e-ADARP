//#pragma once
#ifndef _LOCALSEARCH_H
#define _LOCALSEARCH_H

#include "Data.h"
#include "Decoder.h"

/************************************************************************************
 Method: LocalSearch
 Description: RVND
*************************************************************************************/
TSol LocalSearch(TSol s, int n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> >& cons);

/************************************************************************************
 Method: LS1
 Description: Consecutive node swap
*************************************************************************************/
TSol LS1(TSol s, std::vector <TNode>& node, std::vector <std::vector <double> >& dist, std::vector <TVhcl> vehicle, int nbUsers,
    int nbVehicles, int periodLength, std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> >& cons);

/************************************************************************************
 Method: LS2
 Description: 2-Opt
*************************************************************************************/
TSol LS2(TSol s, std::vector <TNode>& node, std::vector <std::vector <double> >& dist, std::vector <TVhcl> vehicle, int nbUsers,
    int nbVehicles, int periodLength, std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> >& cons);

/************************************************************************************
 Method: LS3
 Description: Relocate
*************************************************************************************/
TSol LS3(TSol s, int n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> >& cons);

/************************************************************************************
 Method: LS4
 Description: Exchange
*************************************************************************************/
TSol LS4(TSol s, int n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> >& cons);

/************************************************************************************
 Method: Insertion
 Description: Attempting to insert one or more uninserted requests
*************************************************************************************/
TSol Insertion(TSol s, int n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> >& cons);

/************************************************************************************
 Method: UpdateLTBackwards
 Description: Updates LT "backwards" from start as long as needed & returns false if infeasibility is encountered.
*************************************************************************************/
bool UpdateLTBackwards(std::vector<TVecSol>& route, int start, std::vector <std::vector <double> > &dist, std::vector <TNode>& node);

/************************************************************************************
 Method: UpdateETForwards
 Description: Updates ET "forwards" from start as long as needed & returns false if infeasibility is encountered.
*************************************************************************************/
bool UpdateETForwards(std::vector<TVecSol>& route, int start, std::vector <std::vector <double> >& dist, std::vector <TNode>& node);

/************************************************************************************
 Method: BestFeasibleInsertion
 Description: Inserts user pick-up and drop-off into cheapest (TRT) feasible position of a route (pickup-first). Returns false if no feasible insertion found.
 ! Modifies route parameter to best feasible route if one is found and does not change the route if no best feasible insertion is found.
*************************************************************************************/
bool BestFeasibleInsertion(std::vector<TVecSol>& route, int user, int nbUsers, std::vector <std::vector <double> > &dist, std::vector <TNode> &node);

/************************************************************************************
 Method: RemoveRequest
 Description: Removes request with pickup at position 'position' in route 'route' + updates ET and LT
 Returns: ID of the deleted request
*************************************************************************************/
int RemoveRequest(std::vector<TVecSol>& route, int position, int nbUsers, std::vector <std::vector <double> >& dist, std::vector <TNode>& node);

void CalculateLoadCapacity(TSol& s, std::vector <TNode>& node);

double rand(double min, double max);
int irand(int min, int max);

#endif