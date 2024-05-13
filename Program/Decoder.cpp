#include "Decoder.h"

// Sort TSol by random-keys
bool sortByRk(const TVecRk& lhs, const TVecRk& rhs) { return lhs.rk < rhs.rk; }

// Sort TSol by user
bool sortByUser(const TVecRk& lhs, const TVecRk& rhs) { return lhs.user < rhs.user; }

TSol Decoder(TSol s, int n, std::vector <TNode> node,
    std::vector <std::vector <double> > dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> > cons)
{

    int numDecoders = 2;

    // Create initial solution (for each vehicle a route between the origin depot and the assigned final depot
    TVecSol tempVecSol;
    std::vector <int> fDepAss = AssignFinalDepots(s.vec, nbUsers, nbVehicles, fDepotIDs, dist);
    for (int i = 0; i < nbVehicles; i++) {
        // Assign origin depot & set correct starting capacities
        tempVecSol.stop = oDepotIDs[i]-1; // Correct origin depot, -1 because indexes start at 0
        //tempVecSol.T = 0.0; // Start-time = 0
        tempVecSol.ET = node[oDepotIDs[i] - 1].arr;
        tempVecSol.LT = node[oDepotIDs[i] - 1].dep;
        tempVecSol.C = vehicle[i].C; // Correct initial capacity
        tempVecSol.B = vehicle[i].B0; // Correct starting state-of-charge
        s.sol[i].push_back(tempVecSol);
        // Assign final depot
        tempVecSol.stop = fDepAss[i];
        tempVecSol.ET = node[fDepAss[i]].arr;
        tempVecSol.LT = node[fDepAss[i]].dep;
        tempVecSol.C = vehicle[i].C; // No change in capacity
        tempVecSol.B = -1.0; // State-of-charge neglected for now
        s.sol[i].push_back(tempVecSol);
    }

    int dec = ceil(s.vec[n].rk*numDecoders + 0.000000000001);
    //printf("\n%d (%.2lf)", dec, s.vec[n].rk);

    s.fo = -1; // Set objective value to -1 for testing feasibility later

    switch (dec)
    {
        case 1: // Ascending insert-all pick-up first
            s = Dec1(s, n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons);
            break;

        case 2: // Ascending insert-all drop-off first
            s = Dec2(s, n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons);
            break;

        default:
            int i = 0;
            break;
    }

    // Sort vec back by user, else vec is not correctly used by the scheduler and the parametric crossover
    sort(s.vec.begin(), s.vec.begin() + nbUsers, sortByUser);

    // if objective value is still -1, it means it has NOT been set to a high value by the decoder to indicate infeasibility
    // So, in this case, scheduling should be applied and objective value should be calculated in the normal way
    if (s.fo == -1) {
        s = Scheduler(s, node, dist, vehicle, nbUsers, nbVehicles, periodLength, cStations, maxRideTimes, weights, elecPrices, cons);
        if (s.scheduled) {
            s.fo = objFct(s, node, dist, nbVehicles, nbUsers, periodLength, maxRideTimes, weights, elecPrices);
        }
    }
    
    return s;
}

TSol Dec1(TSol s, int n, std::vector <TNode> node,
    std::vector <std::vector <double> > dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> > cons)
{
    // sort random-key vector (of the rk's relating to pickups)
    sort(s.vec.begin(), s.vec.begin() + nbUsers, sortByRk);

    // In ascending order, go over users and assign them to vehicles, immediately applying cheapest feasible insertion
    int veh = 0;
    int bestInsert = 0;
    int j = 0;
    float costBest = INFINITO;
    float costInsertion = 0;
    TVecSol precStop;
    TVecSol tempStop{};
    std::vector <TVecSol> tempRoute;
    std::vector <TVecSol> bestRoute;
    std::vector <TVecSol> backUp;
    bool feasible = true;

    for (int i = 0; i < nbUsers; i++) {
        veh = floor(s.vec[i].rk * nbVehicles); // Determine which vehicle to insert into
        bestRoute = s.sol[veh];
        backUp = s.sol[veh];

        for (int nodeID : {s.vec[i].user, nbUsers + s.vec[i].user}) { // Once for pick-up node, once for drop-off node
            costBest = INFINITO;
            // Determine cheapest feasible insertion
            // If drop-off node, start searching where you just inserted the pick-up. Else (if pick-up), start search at first position.
            if (nodeID >= nbUsers) {
                j = bestInsert + 1; // set starting search position to position after where you just inserted the related pickup
            }
            else { j = 1; }
            while (j < s.sol[veh].size()) { // For all possible insertion positions
                feasible = true;
                if (node[nodeID].arr < node[s.sol[veh][j].stop].dep - dist[nodeID][s.sol[veh][j].stop] - node[nodeID].d) { // If soonest time of node you want to insert is earlier than latest time of node that would be visited next if inserted (- travel time & service time). In the other case, insertion will be infeasible anyway.
                    precStop = s.sol[veh][j - 1];
                    // Check if cheapest
                    costInsertion = dist[precStop.stop][nodeID] + dist[nodeID][s.sol[veh][j].stop] - dist[precStop.stop][s.sol[veh][j].stop];
                    // Only check feasibility if cheapest.
                    if (costInsertion < costBest) {
                        // Create a temporary route                  
                        tempRoute = s.sol[veh];
                        tempStop.stop = nodeID;
                        tempRoute.insert(tempRoute.begin() + j, tempStop); // Insert the stop into the temporary route   
                        // Update capacity and check feasibility
                        for (int k = j; k < tempRoute.size(); k++) {
                            tempRoute[k].C = tempRoute[k - 1].C - node[tempRoute[k].stop].l;
                            if (tempRoute[k].C < 0) {
                                feasible = false;
                                break;
                            }
                        }
                        // If not load-feasible, go to next insertion
                        if (!feasible) {
                            j++;
                            continue;
                        }
                        // Update earliest start-time
                        tempRoute[j].ET = std::max(tempRoute[j - 1].ET + node[tempRoute[j - 1].stop].d + dist[tempRoute[j - 1].stop][tempRoute[j].stop], node[tempRoute[j].stop].arr);
                        // Update latest start-time
                        tempRoute[j].LT = std::min(tempRoute[j + 1].LT - node[tempRoute[j].stop].d - dist[tempRoute[j].stop][tempRoute[j + 1].stop], node[tempRoute[j].stop].dep);
                        // If tightened window indicates infeasible, set to infeasibility
                        if (tempRoute[j].LT < tempRoute[j].ET) {
                            feasible = false;
                        }
                        else {
                            // If tightened window does not indicate infeasibility, update ET and LT of other stops and check feasibility again
                            // Update LT "backwards" as long as needed
                            for (int k = j - 1; k >= 0; k--) {
                                if (tempRoute[k + 1].LT - dist[tempRoute[k].stop][tempRoute[k + 1].stop] - node[tempRoute[k].stop].d < tempRoute[k].LT) { // If LT needs to be updated
                                    tempRoute[k].LT = tempRoute[k + 1].LT - dist[tempRoute[k].stop][tempRoute[k + 1].stop] - node[tempRoute[k].stop].d; // Update LT
                                    // Check new tightened window
                                    if (tempRoute[k].LT < tempRoute[k].ET) {
                                        feasible = false;
                                        break;
                                    }
                                }
                                else break; // Else, you can stop backwards updating
                            }
                            // IF STILL FEASIBLE, update ET "forwards" as long as needed
                            if (feasible) {
                                for (int k = j + 1; k < tempRoute.size(); k++) {
                                    if (tempRoute[k - 1].ET + dist[tempRoute[k - 1].stop][tempRoute[k].stop] + node[tempRoute[k - 1].stop].d > tempRoute[k].ET) { // If ET needs to be updated
                                        tempRoute[k].ET = tempRoute[k - 1].ET + dist[tempRoute[k - 1].stop][tempRoute[k].stop] + node[tempRoute[k - 1].stop].d; // Update ET
                                        // Check new tightened window
                                        if (tempRoute[k].LT < tempRoute[k].ET) {
                                            feasible = false;
                                            break;
                                        }
                                    }
                                    else break; // Else, you can stop forwards updating
                                }
                            }
                        }
                        // If temporary route is feasible, copy to best route and set best cost to the cost of this insertion                   
                        if (feasible) {
                            bestRoute = tempRoute;
                            costBest = costInsertion;
                            bestInsert = j;
                        }
                    }
                }
                j++;
            }
            // Check feasibility of this assignment of requests to vehicles
            if (s.sol[veh].size() == bestRoute.size()) {
                // If best route found is as long as current route, no feasible insertion was found.
                // So, objective value is set to high value and we stop decoding further.
                // Because we want to reward attempts that "almost" were feasible, the high value depends on this
                s.fo = (nbUsers - i) * 1000000;
                s.sol[veh] = backUp; // Over-write inserted pick-up if dropoff isn't feasibly insertable.
                return s;

                /*ALTERNATIVE
                if (s.fo == -1) s.fo += 1;
                s.fo += 1000000; // Add 1 million because 1 infeasible insertion found
                bestInsert = 0; // Set to 0 so dropoff can be inserted even if pickup cant
                continue; // Move-on to next node (so, even tho pickup couldn't be inserted, still try to insert drop-off. This way, if both pick-up and drop-off are uninsertable, a larger penalty is set)*/
            }
            // If  feasible inserted, copy best route to solution
            s.sol[veh] = bestRoute;
        }
    }

    return s;
}

TSol Dec2(TSol s, int n, std::vector <TNode> node,
    std::vector <std::vector <double> > dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> > cons)
{
    // sort random-key vector (of the rk's relating to pickups)
    sort(s.vec.begin(), s.vec.begin() + nbUsers, sortByRk);

    // In ascending order, go over users and assign them to vehicles, immediately applying cheapest feasible insertion
    int veh = 0;
    int bestInsert = 0;
    int end = 0;
    float costBest = INFINITO;
    float costInsertion = 0;
    TVecSol precStop;
    TVecSol tempStop{};
    std::vector <TVecSol> tempRoute;
    std::vector <TVecSol> bestRoute;
    std::vector <TVecSol> backUp;
    bool feasible = true;

    for (int i = 0; i < nbUsers; i++) {
        veh = floor(s.vec[i].rk * nbVehicles); // Determine which vehicle to insert into
        bestRoute = s.sol[veh];
        backUp = s.sol[veh];

        for (int nodeID : {nbUsers + s.vec[i].user, s.vec[i].user}) { // Once for drop-off node, once for pick-up node
            costBest = INFINITO;
            // Determine cheapest feasible insertion
            // If pick-up node, end searching where you just inserted the drop-off. Else (if drop-off), end search at last position.
            if (nodeID < nbUsers) {
                end = bestInsert + 1; // set ending search position to position where you just inserted the related drop-off (+1 for range purposes)
            }
            else { end  = s.sol[veh].size(); }
            for (int j = 1; j < end; j++) { // For all possible insertion positions
                feasible = true;
                if (node[nodeID].arr < node[s.sol[veh][j].stop].dep - dist[nodeID][s.sol[veh][j].stop] - node[nodeID].d) { // If soonest time of node you want to insert is earlier than latest time of node that would be visited next if inserted (- travel time & service time). In the other case, insertion will be infeasible anyway.
                    precStop = s.sol[veh][j - 1];
                    // Check if cheapest
                    costInsertion = dist[precStop.stop][nodeID] + dist[nodeID][s.sol[veh][j].stop] - dist[precStop.stop][s.sol[veh][j].stop];
                    // Only check feasibility if cheapest.
                    if (costInsertion < costBest) {
                        // Check load capacity feasibility
                        tempStop.C = precStop.C - node[nodeID].l;
                        if (tempStop.C < 0) {
                            j++; continue;
                        } // If load capacity violated, skip this insertion and check next
                        // Create a temporary route                  
                        tempRoute = s.sol[veh];
                        tempStop.stop = nodeID;
                        tempRoute.insert(tempRoute.begin() + j, tempStop); // Insert the stop into the temporary route 
                        // Update capacity and check feasibility
                        for (int k = j; k < tempRoute.size(); k++) {
                            tempRoute[k].C = tempRoute[k - 1].C - node[tempRoute[k].stop].l;
                            if (tempRoute[k].C < 0) {
                                feasible = false;
                                break;
                            }
                        }
                        // If not load-feasible, go to next insertion
                        if (!feasible) {
                            continue;
                        }
                        // Update earliest start-time
                        tempRoute[j].ET = std::max(tempRoute[j - 1].ET + node[tempRoute[j - 1].stop].d + dist[tempRoute[j - 1].stop][tempRoute[j].stop], node[tempRoute[j].stop].arr);
                        // Update latest start-time
                        tempRoute[j].LT = std::min(tempRoute[j + 1].LT - node[tempRoute[j].stop].d - dist[tempRoute[j].stop][tempRoute[j + 1].stop], node[tempRoute[j].stop].dep);
                        // If tightened window indicates infeasible, set to infeasibility
                        if (tempRoute[j].LT < tempRoute[j].ET) {
                            feasible = false;
                        }
                        else {
                            // If tightened window does not indicate infeasibility, update ET and LT of other stops and check feasibility again
                            // Update LT "backwards" as long as needed
                            for (int k = j - 1; k >= 0; k--) {
                                if (tempRoute[k + 1].LT - dist[tempRoute[k].stop][tempRoute[k + 1].stop] - node[tempRoute[k].stop].d < tempRoute[k].LT) { // If LT needs to be updated
                                    tempRoute[k].LT = tempRoute[k + 1].LT - dist[tempRoute[k].stop][tempRoute[k + 1].stop] - node[tempRoute[k].stop].d; // Update LT
                                    // Check new tightened window
                                    if (tempRoute[k].LT < tempRoute[k].ET) {
                                        feasible = false;
                                        break;
                                    }
                                }
                                else break; // Else, you can stop backwards updating
                            }
                            // IF STILL FEASIBLE, update ET "forwards" as long as needed
                            if (feasible) {
                                for (int k = j + 1; k < tempRoute.size(); k++) {
                                    if (tempRoute[k - 1].ET + dist[tempRoute[k - 1].stop][tempRoute[k].stop] + node[tempRoute[k - 1].stop].d > tempRoute[k].ET) { // If ET needs to be updated
                                        tempRoute[k].ET = tempRoute[k - 1].ET + dist[tempRoute[k - 1].stop][tempRoute[k].stop] + node[tempRoute[k - 1].stop].d; // Update ET
                                        // Check new tightened window
                                        if (tempRoute[k].LT < tempRoute[k].ET) {
                                            feasible = false;
                                            break;
                                        }
                                    }
                                    else break; // Else, you can stop forwards updating
                                }
                            }
                        }
                        // If temporary route is feasible, update load capacities and copy to best route and set best cost to the cost of this insertion  
                        // Update load capacities from insertion onwards
                        for (int k = j + 1; k < tempRoute.size(); k++) {
                            tempRoute[k].C = tempRoute[k - 1].C - node[tempRoute[k].stop].l;
                        }
                        if (feasible) {
                            bestRoute = tempRoute;
                            costBest = costInsertion;
                            bestInsert = j;
                        }
                    }
                }
            }
            // Check feasibility of this assignment of requests to vehicles
            if (s.sol[veh].size() == bestRoute.size()) {
                // If best route found is as long as current route, no feasible insertion was found.
                // So, objective value is set to high value and we stop decoding further.
                // Because we want to reward attempts that "almost" were feasible, the high value depends on this
                s.fo = (nbUsers - i) * 1000000;
                s.sol[veh] = backUp; // Over-write inserted dropoff if pickup can't be inserted.
                return s;

                /*ALTERNATIVE
                if (s.fo == -1) s.fo += 1;
                s.fo += 1000000; // Add 1 million because 1 infeasible insertion found
                bestInsert = 0; // Set to 0 so pickup can be inserted even if drop-off cant
                continue; // Move-on to next node (so, even tho drop-off couldn't be inserted, still try to insert pick-up. This way, if both pick-up and drop-off are uninsertable, a larger penalty is set)*/
            }

            // If feasible inserted, copy best route to solution
            s.sol[veh] = bestRoute;
        }
    }

    return s;
}

double objFct(TSol s, std::vector <TNode>& node, std::vector <std::vector <double> > &dist, int nbVehicles, int nbUsers, int periodLength, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices) {

    // If battery infeasible, no point in calculating objective function
    if (s.battery_infeasibles > 0) {
        return s.battery_infeasibles * 100000;
    }

    double TRT = 0; // Total Ride Time
    double ERT = 0; // Excess Ride Time
    double CC = 0;  // Charging cost
    double penalty = 0; // Penalty for URT violation
    std::vector <double> URTs; // User Ride Times
    int currStop = -1;
    int nextStop = -1;
    int startPeriod;
    int endPeriod;

    for (int u = 0; u < nbUsers; u++) { // Initialize User Ride Times by substracting service time already
        URTs.push_back(-1 * node[u].d);
    }

    for (int v = 0; v < nbVehicles; v++) { // Loop over vehicles
        // First: origin depot
        currStop = s.sol[v][0].stop;
        nextStop = s.sol[v][1].stop;
        TRT += dist[currStop][nextStop]; // Add travel time between current stop and next stop
        // Loop over stops in route of vehicle, except final depot and origin depot
        for (int i = 1; i < s.sol[v].size() - 1; i++) {
            currStop = nextStop;
            nextStop = s.sol[v][i + 1].stop;
            TRT += dist[currStop][nextStop]; // Add travel time between current stop and next stop
            if (currStop < nbUsers) { // If current stop is a pick-up
                URTs[currStop] -= s.sol[v][i].T; // Subtract start-time of service at pick-up
            }
            else if (currStop < nbUsers * 2) { // If current stop is a drop-off
                URTs[currStop-nbUsers] += s.sol[v][i].T; // Add start-time of service at pick-up
            }
            else { // Else, current stop is a charging station
                // Find period in which charging starts and ends
                startPeriod = floor(s.sol[v][i].T / periodLength);
                endPeriod = floor((s.sol[v][i].T + s.sol[v][i].w) / periodLength);
                if (endPeriod == elecPrices.size()) { endPeriod = endPeriod - 1; }  // This is to prevent subscript errors when charging ends exactly at the time horizon -> floor function above acts as if charging ends after H, while it ends just at H (possible because charging stations can coincide with destination depots)
                if (startPeriod == elecPrices.size()) { startPeriod = startPeriod - 1; } // Similar but for case where 'charging' starts at H (nonsensical case but possible because charging stations can coincide
                // If whole session in 1 period
                if (startPeriod == endPeriod) {
                    CC += s.sol[v][i].w * elecPrices[startPeriod];
                }
                else {
                    // Calculate cost during first and last period
                    CC += ((startPeriod + 1) * periodLength - s.sol[v][i].T) * elecPrices[startPeriod];
                    CC += (s.sol[v][i].T + s.sol[v][i].w - endPeriod * periodLength) * elecPrices[endPeriod];
                    // Calculate cost for inbetween periods
                    for (int j = startPeriod + 1; j < endPeriod; j++) {
                        CC += periodLength * elecPrices[j];
                    }
                }
                
            }
        }
    }

    for (int i = 0; i < nbUsers; i++) {
        // Check if Max. URT is violated, then return INFINITO. Else, calculate ERT.
        if (URTs[i] > maxRideTimes[i]) {
            penalty += 1;
        }
        ERT += URTs[i] - dist[i][i + nbUsers]; // Actual Ride Time - Ideal Ride Time
    }

    /////////////////////////// FOR DEBUG PURPOSES //////////////////////////////////////////////////////////
   /* bool dropOffFound;
    for (int v = 0; v < nbVehicles; v++) { // For each vehicle
        for (int i = 1; i < s.sol[v].size() - 1; i++) { // For each stop in the route
            dropOffFound = false;
            if (s.sol[v][i].stop < nbUsers) { // If pickup
                for (int j = i + 1; j < s.sol[v].size() - 1; j++) { // For all following nodes
                    if (s.sol[v][j].stop == s.sol[v][i].stop + nbUsers) { // If corresponding dropoff
                        dropOffFound = true;
                    }
                }
                if (!dropOffFound) {
                    int debug = 0; // DROPOFF WAS NOT FOUND IN THE ROUTE IT SHOULD BE IN
                    return weights[0] * TRT + weights[1] * ERT + weights[2] * CC + penalty * 10000;
                }
            }
            
        }
    }*/
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    return weights[0] * TRT + weights[1] * ERT + weights[2] * CC + penalty * 10000;
}

std::vector <int> AssignFinalDepots(std::vector<TVecRk> rkVec, int nbUsers, int nbVehicles, std::vector<int> fDepotIDs, std::vector <std::vector <double> > dist) {
    std::vector <int> fDepotAssigned(nbVehicles, -1);
    int matchingDepot = -1;
    sort(rkVec.begin() + 2 * nbUsers, rkVec.begin() + 2 * nbUsers + nbVehicles, sortByRk); // Sort vehicles by random key

    // Go over each vehicle and assign the matching final depot, if not yet assigned, else closest to matching depot
    for (int i = 2*nbUsers; i < 2*nbUsers + nbVehicles; i++) {
        matchingDepot = fDepotIDs[floor(rkVec[i].rk * fDepotIDs.size())] -1; // -1 because first node has index 0 instead of 1 and so on
        if (std::count(fDepotAssigned.begin(), fDepotAssigned.end(), matchingDepot)) { // if depot has already been assigned
            // Find closest final depot that has not been assigned yet
            double smallestDist = INFINITO;
            int newMatchDepot = matchingDepot;
            for (int fDepotID : fDepotIDs) {
                if (dist[matchingDepot][fDepotID - 1] < smallestDist) { // IF smallest dist
                    if (!std::count(fDepotAssigned.begin(), fDepotAssigned.end(), fDepotID - 1)){ // IF not yet assigned
                        newMatchDepot = fDepotID - 1;
                        smallestDist = dist[matchingDepot][fDepotID - 1];
                    }
                }
            }
            // Set matchingDepot to this closest found depot
            matchingDepot = newMatchDepot;
        }
        // Assign the matching depot to the vehicle
        fDepotAssigned[rkVec[i].user] = matchingDepot;
    }

    sort(rkVec.begin() + 2 * nbUsers, rkVec.begin() + 2 * nbUsers + nbVehicles, sortByUser); // Sort vehicles by user again (to not hinder cross-over)

    return fDepotAssigned;
}
