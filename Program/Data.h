//#pragma once
#ifndef _DATA_H
#define _DATA_H

//using namespace std;
#include <vector>
#include <algorithm>    


//------ DEFINITION OF TYPES OF PROBLEM SPECIFIC --------

struct TNode								// struct with node informations
{
	int id;
	double x; // Lat
	double y; // Long
    float d; // Service time
    int l; // Load change
    double arr; // Earliest time of starting service
    double dep; // Latest time of starting service
};

struct TVhcl                                // struct with vehicle informations
{
    int id; // Vehicle ID
    int C; // (Load) Capacity
    double B0; // Initial battery inventory
    double Q; // (Effective) Battery capacity
    double r; // Minimum end battery ratio level
};

struct TCstat                               // struct with charging station specific informations
{
    int id;
    double alpha; // Recharging rate
    int cap; // Number of chargers at this station
};


//------ DEFINITION OF TYPES OF BRKGA-QL --------

/***********************************************************************************
 Struct: TVecRk
 Description: struct to represent the random keys
************************************************************************************/
struct TVecRk
{
    int user;                                // ID of the pick-up from the request
    double rk;                               // random-key of chromosome
};

/***********************************************************************************
 Struct: TVecSol
 Description: struct to represent the solution
************************************************************************************/
struct TVecSol
{
    int stop;                                  // ID of the node from the request
    double T;                                  // Start-time of service at stop
    double w=0;                                  // charging duration at charging station
    double ET;                                 // Tightened earliest service start-time
    double LT;                                 // Tightened latest service start-time
    int C;                                     // Load capacity after service at stop
    double B;                                  // State-of-charge after service at stop
    bool c_station = false;               //indicates if a node is a charging station or not
};


/***********************************************************************************
 Struct: TSol
 Description: struct to represent a solution problem
************************************************************************************/
struct TSol
{
    std::vector <TVecRk> vec;               // id of user/pick-up and random key
    std::vector <std::vector <TVecSol>> sol;// matrix to store solution routes + schedule
    int battery_infeasibles = 0;            // Auxiliary element to indicate how many battery infeasibilities were incurred. So, obj fct can be adjusted to this number
    bool scheduled = false;                 // Auxiliary element for the schedulers
    double fo;                              // objetive function value
    int label;                              // defines a community solution with a number
    int similar;                            // indicates if a solution is similar to other (0 no, 1 yes)
    int flag;                               // indicates if a local search has already been performed on this solution (0 no, 1 yes)
    int promising;                          // indicates if a solution is promising to apply local search
    double TRT;
    double ERT;
    double CC;
};


/***********************************************************************************
 Struct: TQ
 Description: struct to represent a quality matrix
************************************************************************************/
struct TQ
{
    int S;
    double pVar;
    double q;
    int k;
};


#endif