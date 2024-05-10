#include "LocalSearch.h"

// Sort TSol by random-keys
bool sortByRkLS(const TVecRk& lhs, const TVecRk& rhs) { return lhs.rk < rhs.rk; }

// Sort TSol by user
bool sortByUserLS(const TVecRk& lhs, const TVecRk& rhs) { return lhs.user < rhs.user; }

TSol LocalSearch(TSol s, int n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> >& cons)
{
    // If obj value > 999999, not all requests could feasibly be inserted
    // Try to insert un-inserted requests in other routes than the rk-decided route
    if (s.fo > 999999) {
        s = Insertion(s, n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons);
        if (s.fo > 999999) { return s; } // If still not feasible, just return
        // Else, if it IS feasible now, continue with LS
    }

    // Delete all charging stations that were inserted during original scheduling + reset T and w variables (just to be sure)
    for (int i = 0; i < nbVehicles; i++) {
        for (int j = 0; j < s.sol[i].size(); j++) {
            if (s.sol[i][j].stop >= nbUsers * 2) { // If not a pick-up or drop-off
                if (s.sol[i][j].c_station == true) { // If charging station
                    s.sol[i].erase(s.sol[i].begin() + j); // Remove stop
                    j--; // Set j one back, since now a new element stands at spot j
                }
                else { // If depot
                    s.sol[i][j].T = 0;
                    s.sol[i][j].w = 0;
                }
            }
            else {
                s.sol[i][j].T = 0;
                s.sol[i][j].w = 0;
            }
        }
    }

    // Also reset ET and LT variables, since these were possibly changed during the first scheduling
    for (int i = 0; i < nbVehicles; i++) {
        s.sol[i][0].ET = node[s.sol[i][0].stop].arr; // Origin depot
        s.sol[i][0].LT = node[s.sol[i][0].stop].dep; // Origin depot
        s.sol[i][s.sol[i].size()-1].ET = node[s.sol[i][s.sol[i].size()-1].stop].arr; // Final depot
        s.sol[i][s.sol[i].size()-1].LT = node[s.sol[i][s.sol[i].size()-1].stop].dep; // Final depot
        if (!UpdateETForwards(s.sol[i], 1, dist, node) || !UpdateLTBackwards(s.sol[i], s.sol[i].size() - 2, dist, node)) {
            // This point should not be reached. If it is, it means a ET-LT overlap has indicated an infeasibility, which should not be possible
            int error = 666;
        }
    }

    // we use a Random Variable Neighborhood Descent (RVND) as local search
    int k = 1;

    // predefined number of neighborhood moves
    std::vector <int> NSL;
    std::vector <int> NSLAux;
    for (int i=1; i<=4; i++)
    {
        NSL.push_back(i);
        NSLAux.push_back(i);
    }

    //printf("\nHeuristics: ");
    while (!NSL.empty())
	{

        // current objective function
        double foCurrent = s.fo;

        // randomly choose a neighborhood
        int pos = irand(0,NSL.size()-1);
        k = NSL[pos];

        switch (k)
        {
        case 1: 
            s = LS1(s, node, dist, vehicle, nbUsers, nbVehicles, periodLength, cStations, maxRideTimes, weights, elecPrices, cons); // Consecutive node swap
            break;

        case 2:
            s = LS2(s, node, dist, vehicle, nbUsers, nbVehicles, periodLength, cStations, maxRideTimes, weights, elecPrices, cons); // Random feasible 2-opt
            break;

        case 3:
            s = LS3(s, n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons); // Exhaustive relocate;
            break;

        case 4:
            s = LS4(s, n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons); // Random feasible exchange
        
        default:
            break;
        }

        // return to first neighborhood if better the current solution
        if (s.fo < foCurrent)
        {
            // refresh NSL
            NSL.clear();
            NSL = NSLAux;
        }
        // next neighborhood, otherwise
        else
        {
            // Remove N(n) from NSL
            NSL.erase(NSL.begin()+pos);
        }
	} //end while

    // Final version has NO schedule and charging stops, in order to keep iteratively applying LS possible.
    // However, fo has been set to fo AFTER scheduling
    // Before returning, schedule has to be set again & charging station need to be implemented
    s = Scheduler(s, node, dist, vehicle, nbUsers, nbVehicles, periodLength, cStations, maxRideTimes, weights, elecPrices, cons);
 	return s;
}

TSol LS1(TSol s, std::vector <TNode> &node, std::vector <std::vector <double> >& dist, std::vector <TVhcl> vehicle, int nbUsers,
    int nbVehicles, int periodLength, std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> > &cons)
{
    TSol tempSol;
    TSol schedSol; // A copy of the solution to schedule. Else, a solution WITH charging stations could be passed on to next LS, causing issues
    TSol bestSol = s;
    std::vector<TVecSol> tempRoute;
    int aux;
    bool feasible;

    // For each stop
    for (int veh = 0; veh < nbVehicles; veh++) {
        for (int i = 1; i < s.sol[veh].size() - 2; i++) { // Don't swap first and last element
            // Reset tempSol and tempRoute
            tempSol = s;
            tempRoute = s.sol[veh];

            // Perform the swap with the next stop (if not the corresponding drop-off)
            if (tempRoute[i + 1].stop != tempRoute[i].stop + nbUsers) {

                // Swap stops
                aux = tempRoute[i].stop;
                tempRoute[i].stop = tempRoute[i + 1].stop;
                tempRoute[i + 1].stop = aux;

                // Recalculate ET and LT for these stops
                tempRoute[i].ET = std::max(tempRoute[i - 1].ET + node[tempRoute[i - 1].stop].d + dist[tempRoute[i - 1].stop][tempRoute[i].stop], node[tempRoute[i].stop].arr);
                tempRoute[i+1].LT = std::min(tempRoute[i+2].LT - node[tempRoute[i+1].stop].d - dist[tempRoute[i+1].stop][tempRoute[i+2].stop], node[tempRoute[i+1].stop].dep);
                tempRoute[i].LT = std::min(tempRoute[i + 1].LT - node[tempRoute[i].stop].d - dist[tempRoute[i].stop][tempRoute[i + 1].stop], node[tempRoute[i].stop].dep);
                tempRoute[i+1].ET = std::max(tempRoute[i].ET + node[tempRoute[i].stop].d + dist[tempRoute[i].stop][tempRoute[i+1].stop], node[tempRoute[i+1].stop].arr);

                // If tightened window indicates infeasible, go to next swap
                if (tempRoute[i].LT < tempRoute[i].ET || tempRoute[i+1].LT < tempRoute[i+1].ET) {
                    continue;
                }
                else {
                    feasible = true;
                    // If tightened window does not indicate infeasibility, update ET and LT of other stops and check feasibility again
                    // Update LT "backwards" as long as needed
                    for (int k = i - 1; k > 0; k--) {
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
                        for (int k = i + 2; k < tempRoute.size(); k++) {
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
                    if (!feasible) continue; // If not a feasible swap, continue to next swap
                }

                tempSol.sol[veh] = tempRoute;

                // Update load capacity
                CalculateLoadCapacity(tempSol, node);

                // Apply scheduler
                schedSol = Scheduler(tempSol, node, dist, vehicle, nbUsers, nbVehicles, periodLength, cStations, maxRideTimes, weights, elecPrices, cons);

                // Calculate objective function
                tempSol.fo = objFct(schedSol, node, dist, nbVehicles, nbUsers, periodLength, maxRideTimes, weights, elecPrices);

                // If obj function has improved, copy this solution to s
                if (tempSol.fo < bestSol.fo) {
                    bestSol = tempSol;
                }
            }
        }
    }
    return bestSol;
}

TSol LS2(TSol s, std::vector <TNode>& node, std::vector <std::vector <double> >& dist, std::vector <TVhcl> vehicle, int nbUsers,
    int nbVehicles, int periodLength, std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> >& cons)
{
    // Choose (different) 2 random routes
    int route1 = irand(0, s.sol.size() - 1);
    int route2 = irand(0, s.sol.size() - 1);
    while (route1 == route2) route2 = irand(0, s.sol.size() - 1);

    // Identify all zero-load points in each route
    std::vector<int> zeroLoadsRoute1;
    std::vector<int> zeroLoadsRoute2;
    for (int i = 1; i < s.sol[route1].size() - 1; i++) { // Find all zero-load points in route 1 (except depots)
        if (s.sol[route1][i].C == vehicle[route1].C) {
            zeroLoadsRoute1.push_back(i);
        }
    }
    for (int i = 1; i < s.sol[route2].size() - 1; i++) { // Find all zero-load points in route 2 (except depots)
        if (s.sol[route2][i].C == vehicle[route2].C) {
            zeroLoadsRoute2.push_back(i);
        }
    }

    // If one of the routes does not have zero-load points, stop
    if (zeroLoadsRoute1.size() == 0 || zeroLoadsRoute2.size() == 0) {
        return s;
    }

    // Identify (potentially) feasible zero-load points combo's
    std::vector<std::vector<int>> zeroLoadCombos;
    bool feasible = false;
    for (int i : zeroLoadsRoute1) {
        for (int j : zeroLoadsRoute2) {
            feasible = s.sol[route1][i].ET + node[s.sol[route1][i - 1].stop].d + dist[s.sol[route1][i - 1].stop][s.sol[route1][i].stop] <= s.sol[route2][j + 1].LT; // Earliest possible time of arrival at new next stop <= LT of new next stop
            feasible = feasible && s.sol[route2][j].ET + node[s.sol[route2][j - 1].stop].d + dist[s.sol[route2][j - 1].stop][s.sol[route2][j].stop] <= s.sol[route1][i + 1].LT; // Idem, but route1 and 2 reversed.
            if (feasible) {
                zeroLoadCombos.push_back({ i, j });
            }
        }
    }

    // If no potentially feasible combo's, stop
    if (zeroLoadCombos.size() == 0) {
        return s;
    }

    TSol bestSol = s;
    TSol schedSol; // A copy of the solution to schedule. Else, a solution WITH charging stations could be passed on to next LS, causing issues

    // Check each feasible combo
    for (std::vector<int> zeroLoadCombo : zeroLoadCombos) {
        // Copy current solution
        TSol tempSol = s;

        // Delete all stops after zero-load point in copy
        int nbToPop = tempSol.sol[route1].size() - zeroLoadCombo[0] - 1;
        for (int i = 0; i < nbToPop; i++) {
            tempSol.sol[route1].pop_back();
        }
        nbToPop = tempSol.sol[route2].size() - zeroLoadCombo[1] - 1;
        for (int i = 0; i < nbToPop; i++) {
            tempSol.sol[route2].pop_back();
        }

        // Copy stops from after zero-load points from original solution
        // Route 1 to route 2 and vice versa
        for (int i = zeroLoadCombo[0] + 1; i < s.sol[route1].size(); i++) {
            tempSol.sol[route2].push_back(s.sol[route1][i]);
        }
        for (int i = zeroLoadCombo[1] + 1; i < s.sol[route2].size(); i++) {
            tempSol.sol[route1].push_back(s.sol[route2][i]);
        }

        // Update ET and LT
        //// Route1
        ////// Forwards ET updating
        if (!UpdateETForwards(tempSol.sol[route1], zeroLoadCombo[0], dist, node)) return s; // If infeasibility found, return s
        ////// Backwards LT updating
        if (!UpdateLTBackwards(tempSol.sol[route1], zeroLoadCombo[0], dist, node)) return s; // If infeasibility found, return s

        //// Route2
        ////// Forwards ET updating
        if (!UpdateETForwards(tempSol.sol[route2], zeroLoadCombo[1], dist, node)) return s; // If infeasibility found, return s
        ////// Backwards LT updating
        if (!UpdateLTBackwards(tempSol.sol[route2], zeroLoadCombo[1], dist, node)) return s; // If infeasibility found, return s

        // Update load capacities
        CalculateLoadCapacity(tempSol, node);

        // Schedule
        schedSol = Scheduler(tempSol, node, dist, vehicle, nbUsers, nbVehicles, periodLength, cStations, maxRideTimes, weights, elecPrices, cons);

        // Calculate obj function
        tempSol.fo = objFct(schedSol, node, dist, nbVehicles, nbUsers, periodLength, maxRideTimes, weights, elecPrices);

        // Return new solution only if better objective function
        if (tempSol.fo < bestSol.fo) {
            bestSol = tempSol;
        }
    }

    return bestSol;
}

TSol LS3(TSol s, int n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> >& cons) {

    // Choose 1 random route
    int route = irand(0, s.sol.size() - 1);

    // Iterate over all users in this route and try inserting into different route
    TSol bestSol = s;
    TSol tempSol;
    TSol schedSol; // A copy of the solution to schedule. Else, a solution WITH charging stations could be passed on to next LS, causing issues
    TSol auxSol;
    TSol auxBest;

    for (int i = 1; i < s.sol[route].size() - 1; i++) { // For each user in the route (except depots)
        tempSol = s; // Reset tempSol
        if (s.sol[route][i].stop < nbUsers) { // If pick-up node

            // Delete pick-up and drop-off from route
            RemoveRequest(tempSol.sol[route], i, nbUsers, dist, node);

            // Try inserting in different routes
            auxBest = bestSol;
            for (int veh = 0; veh < s.sol.size(); veh++) { // For each route in which to try insertion
                if (veh == route) continue;
                else {
                    auxSol = tempSol;
                    // Insert user in cheapest feasible place
                    if (BestFeasibleInsertion(auxSol.sol[veh], s.sol[route][i].stop, nbUsers, dist, node)) { // If user could be feasibly inserted in the route
                        CalculateLoadCapacity(auxSol, node);
                        schedSol = Scheduler(auxSol, node, dist, vehicle, nbUsers, nbVehicles, periodLength, cStations, maxRideTimes, weights, elecPrices, cons); // Schedule auxSol
                        auxSol.fo = objFct(schedSol, node, dist, nbVehicles, nbUsers, periodLength, maxRideTimes, weights, elecPrices); // Calculate new objective value
                        if (auxSol.fo < auxBest.fo) { // If insertion of user 'i' into route 'veh' is new best, save. 
                            auxBest = auxSol;
                        }
                    }
                }
            }

            // If best solution after trying to insert user 'i' into all routes is new best, save it.
            if (auxBest.fo < bestSol.fo) {
                bestSol = auxBest;
            }
        }
    }

    // After having tried re-inserting every user from the route into other routes, the best found solution will be returned.
    // If no feasible solution found, just the inputted solution is returned.
    return bestSol;
}

TSol LS4(TSol s, int n, std::vector <TNode> &node,
    std::vector <std::vector <double> > &dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> >& cons) {

    TSol tempSol = s;
    TSol schedSol;

    // Exception catching
    // If there is only 1 route with at least 1 request, the exchange LS operator is impossible
    int numberOfNotEmptyRoutes = 0;
    for (int i = 0; i < s.sol.size(); i++) {
        if (s.sol[i].size() > 3) { numberOfNotEmptyRoutes++; }
        if (numberOfNotEmptyRoutes > 1) { break; }
    }
    if (numberOfNotEmptyRoutes < 2) {
        return s;
    }

    // Choose 2 random routes that are not the same and have at least 1 request
    int route1 = irand(0, s.sol.size() - 1);
    while (s.sol[route1].size() < 4) {
        route1 = irand(0, s.sol.size() - 1);
    }
    int route2 = irand(0, s.sol.size() - 1);
    while (route1 == route2 || s.sol[route2].size() < 4) {
        route2 = irand(0, s.sol.size() - 1);
    }

    // Choose INDEX of 1 random request/user per route (not ID)
    int user1 = irand(1, s.sol[route1].size() - 2); // Don't pick origin of final depot
    if (s.sol[route1][user1].stop >= nbUsers) { // If you selected a drop-off
        // Find position of pick-up
        for (int i = user1 - 1; i > 0; i--) {
            if (s.sol[route1][i].stop == s.sol[route1][user1].stop - nbUsers){ // If dropoff found
                user1 = i;
                break;
            }
        }
    }
    int user2 = irand(1, s.sol[route2].size() - 2); // Don't pick origin of final depot
    if (s.sol[route2][user2].stop >= nbUsers) { // If you selected a drop-off
        // Find position of pick-up
        for (int i = user2 - 1; i > 0; i--) {
            if (s.sol[route2][i].stop == s.sol[route2][user2].stop - nbUsers) { // If dropoff found
                user2 = i;
                break;
            }
        }
    }

    // Delete requests from routes
    int user1ID = RemoveRequest(tempSol.sol[route1], user1, nbUsers, dist, node);
    int user2ID = RemoveRequest(tempSol.sol[route2], user2, nbUsers, dist, node);

    // Insert user1 into route2
    if (BestFeasibleInsertion(tempSol.sol[route2], user1ID, nbUsers, dist, node)) {
        // Only if feasible insertion found for user1 in route2, insert user2 into route1
        if (!BestFeasibleInsertion(tempSol.sol[route1], user2ID, nbUsers, dist, node)) {
            // If no feasible insertion was found for user2 in route1, return s
            return s;
        }
    }
    else { // If no feasible insertion found, return s
        return s;
    }

    CalculateLoadCapacity(tempSol, node);
    schedSol = Scheduler(tempSol, node, dist, vehicle, nbUsers, nbVehicles, periodLength, cStations, maxRideTimes, weights, elecPrices, cons);
    tempSol.fo = objFct(schedSol, node, dist, nbVehicles, nbUsers, periodLength, maxRideTimes, weights, elecPrices);
    
    if (tempSol.fo < s.fo) {
        return tempSol;
    }
    else {
        return s;
    }
}

TSol Insertion(TSol s, int n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl> vehicle,
    int nbVehicles, int nbUsers, int periodLength, int H, std::vector <int> oDepotIDs, std::vector <int> fDepotIDs,
    std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> >& cons) {

    int nbUninsertedUsers = s.fo / 1000000;
    TSol schedSol; // A copy of the solution to schedule. Else, a solution WITH charging stations could be passed on to next LS, causing issues

    sort(s.vec.begin(), s.vec.begin() + nbUsers, sortByRkLS); // need to sort again, since at the end of decoding, rk's are sorted by user again

    // Iterate over all users that weren't inserted
    for (int i = nbUsers - nbUninsertedUsers; i < nbUsers; i++) { // For each user still not inserted. i is index of place within rk vector, NOT the user index itself
        // Try inserting in different routes
        for (int veh = 0; veh < s.sol.size(); veh++) { // For each route in which to try insertion
            if (veh == floor(s.vec[i].rk * nbVehicles)) continue; // Skip route that was already tried
            else {
                // Try insertion
                if (BestFeasibleInsertion(s.sol[veh], s.vec[i].user, nbUsers, dist, node)) {
                    s.fo -= 1000000; // One less infeasibility
                    break; // If feasible insertion found, stop searching for this user
                }
            }
        }
    }

    // If all users are now inserted, apply scheduling & calculate objective value
    if (s.fo == 0) {
        CalculateLoadCapacity(s, node);
        schedSol = Scheduler(s, node, dist, vehicle, nbUsers, nbVehicles, periodLength, cStations, maxRideTimes, weights, elecPrices, cons);
        s.fo = objFct(schedSol, node, dist, nbVehicles, nbUsers, periodLength, maxRideTimes, weights, elecPrices);
    }

    sort(s.vec.begin(), s.vec.begin() + nbUsers, sortByUserLS); // To make parametric crossover work

    return s;
}


bool UpdateLTBackwards(std::vector<TVecSol>& route, int start, std::vector <std::vector <double> >& dist, std::vector <TNode>& node) {
    // Updates LT "backwards" from start as long as needed & returns false if infeasibility is encountered.

    for (int k = start; k >= 0; k--) {
        for (int i = start; i > 0; i--) {
            route[i].LT = std::min(node[route[i].stop].dep, route[i + 1].LT - node[route[i].stop].d - dist[route[i].stop][route[i + 1].stop]);
            // Check if still feasible
            if (route[i].ET > route[i].LT) {
                // If not feasible, return false
                return false;
            }
        }
    }
    return true; // If no infeasibility encountered, return TRUE
}

bool UpdateETForwards(std::vector<TVecSol>& route, int start, std::vector <std::vector <double> >& dist, std::vector <TNode>& node) {
    for (int i = start; i < route.size(); i++) {
        route[i].ET = std::max(node[route[i].stop].arr, route[i - 1].ET + node[route[i - 1].stop].d + dist[route[i - 1].stop][route[i].stop]);
        // Check if still feasible
        if (route[i].ET > route[i].LT) {
            // If not feasible, return false
            return false;
        }
    }
    // If no infeasibility encountered, return true
    return true;
}

int RemoveRequest(std::vector<TVecSol>& route, int position, int nbUsers, std::vector <std::vector <double> > &dist, std::vector <TNode> &node) {
    int requestID = route[position].stop;
    // Delete pick-up and drop-off from route
    route.erase(route.begin() + position); // Delete pick - up
    UpdateLTBackwards(route, position - 1, dist, node); // Update LT
    UpdateETForwards(route, position, dist, node); // Update ET
    // Find and delete drop-off
    for (int j = position; j < route.size() - 1; j++) { // Start searching AT position where (deleted) pick-up WAS, since now the next element is at that position
        if (route[j].stop == requestID + nbUsers) { // If this is the drop-off you're searching for
            route.erase(route.begin() + j); // Delete drop-off
            UpdateLTBackwards(route, j - 1, dist, node); // Update LT
            UpdateETForwards(route, j, dist, node); // Update ET
            break; // Stop looking further
        }
    }
    return requestID;
}

bool BestFeasibleInsertion(std::vector<TVecSol>& route, int user, int nbUsers, std::vector <std::vector <double> > &dist, std::vector <TNode> &node) {
    int bestInsert = 0;
    int j;
    float costBest;
    float costInsertion = 0;
    TVecSol precStop;
    TVecSol tempStop{};
    std::vector <TVecSol> tempRoute;
    std::vector <TVecSol> bestRoute = route;
    std::vector <TVecSol> backUp = route;
    bool feasible = true;

    for (int nodeID : {user, nbUsers + user}) { // Once for pick-up node, once for drop-off node
        costBest = INFINITO;
        // Determine cheapest feasible insertion
        // If drop-off node, start searching where you just inserted the pick-up. Else (if pick-up), start search at first position.
        j = bestInsert + 1; // If pick-up node, j will = 1. If drop-off, j will = location of pick-up +1.

        while (j < route.size()) { // For all possible insertion positions. +1 for dropoff nodes, because pickup was just inserted in bestRoute (but not in route yet)
            feasible = true;
            if (node[nodeID].arr < route[j].LT - dist[nodeID][route[j].stop] - node[nodeID].d) { // If soonest time of node you want to insert is earlier than latest time of node that would be visited next if inserted (- travel time & service time). In the other case, insertion will be infeasible anyway.
                precStop = route[j - 1];
                // Check if cheapest
                costInsertion = dist[precStop.stop][nodeID] + dist[nodeID][route[j].stop] - dist[precStop.stop][route[j].stop];
                // Only check feasibility if cheapest.
                if (costInsertion < costBest) {
                    // Create a temporary route                  
                    tempRoute = route;
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
                        feasible = UpdateLTBackwards(tempRoute, j - 1, dist, node);
                        // IF STILL FEASIBLE, update ET "forwards" as long as needed
                        if (feasible) {
                            feasible = UpdateETForwards(tempRoute, j + 1, dist, node);
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
        if (route.size() == bestRoute.size()) { // If pick-up, this tests whether bestRoute is as big as original route. If drop-off, this tests whether bestRoute is as big as original route - 1 for already inserted pickup.
            // If best route found is as long as current route, no feasible insertion was found.
            route = backUp;
            return false;
        }
        // If feasibly inserted, copy to route
        route = bestRoute;
    }
    // If pickup AND dropoff feasibly inserted return true
    return true;
}

void CalculateLoadCapacity(TSol& s, std::vector <TNode>& node) {
    for (int v = 0; v < s.sol.size(); v++) {
        for (int i = 1; i < s.sol[v].size()-1; i++) {
            s.sol[v][i].C = s.sol[v][i - 1].C - node[s.sol[v][i].stop].l;
        }
    }
}

double rand(double min, double max)
{
	return ((double)(rand()%10000)/10000.0)*(max-min)+min;
    //return uniform_real_distribution<double>(min, max)(rng);
}

int irand(int min, int max)
{
	return (int)rand(0,max-min+1.0) + min;
}
