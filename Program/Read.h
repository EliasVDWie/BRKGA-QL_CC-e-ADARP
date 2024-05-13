//#pragma once
#ifndef _READ_H
#define _READ_H

#include "Data.h"

#include <vector>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>

void ReadData(char nameTable[], int& n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl>& vehicle,
    int& nbVehicles, int& nbUsers, int& periodLength, int& H, std::vector <int>& oDepotIDs, std::vector <int>& fDepotIDs,
    std::vector <TCstat>& cStations, std::vector <int>& maxRideTimes, double(&weights)[3], std::vector <double>& elecPrices, std::vector <std::vector <double> >& cons);

void ReadDataCordeau(char nameTable[], int& n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl>& vehicle,
    int& nbVehicles, int& nbUsers, int& periodLength, int& H, std::vector <int>& oDepotIDs, std::vector <int>& fDepotIDs,
    std::vector <TCstat>& cStations, std::vector <int>& maxRideTimes, double(&weights)[3], std::vector <double>& elecPrices, std::vector <std::vector <double> >& cons);

void ReadDataUber(char nameTable[], int& n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl>& vehicle,
    int& nbVehicles, int& nbUsers, int& periodLength, int& H, std::vector <int>& oDepotIDs, std::vector <int>& fDepotIDs,
    std::vector <TCstat>& cStations, std::vector <int>& maxRideTimes, double(&weights)[3], std::vector <double>& elecPrices, std::vector <std::vector <double> >& cons);

void FreeMemoryProblem(std::vector <TNode>& node, std::vector <std::vector <double> >& dist, std::vector <TVhcl>& vehicle, std::vector <int>& oDepotIDs, std::vector <int>& fDepotIDs,
    std::vector <TCstat>& cStations, std::vector <int>& maxRideTimes, std::vector <double>& elecPrices, std::vector <std::vector <double> >& cons);

#endif