#include "Scheduler.h"
#include "Decoder.h"

// Sort Zeroloads by ElectricityPrice
bool sortByEP(const ZL& lhs, const ZL& rhs) { return lhs.price < rhs.price; }

// Sort Zeroloads by RK
bool sortByRKsched(const ZL& lhs, const ZL& rhs) { return lhs.rk > rhs.rk; }

// Sort Zeroloads by Window
bool sortByWindow(const ZL& lhs, const ZL& rhs) { return lhs.window > rhs.window; }

//Scheduling procedure based on Bongiovanni 2020
TSol Scheduler(TSol s, std::vector <TNode>& node, std::vector <std::vector <double> >& dist, std::vector <TVhcl> vehicle, int nbUsers,
    int nbVehicles, int periodLength, std::vector <TCstat> cStations, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices, std::vector <std::vector <double> >& cons)
{
    int numChargingPriorities = 5;
    int chargingPriority = ceil(s.vec[s.vec.size() - 3].rk * numChargingPriorities + 0.000000000001);
    double tempRK;

    // Create structure to store charger availabilities
    int nbOfHours = ceil(elecPrices.size() * periodLength / 60.0);
    std::vector<unsigned long long> charger(nbOfHours, 0);
    std::vector<std::vector<unsigned long long>> chargingStation;
    std::vector<std::vector<std::vector<unsigned long long>>> availability;
    for (int i = 0; i < cStations.size(); i++) {
        for (int j = 0; j < cStations[i].cap; j++) { // add cap ampount of chargers
            chargingStation.push_back(charger);
        }
        availability.push_back(chargingStation);
        chargingStation.clear();
    }

    for (int i = 0; i < nbVehicles; i++) {
        s.sol[i][0].B = vehicle[i].B0;
        for (int j = 1; j < s.sol[i].size(); j++) {
            s.sol[i][j].B = s.sol[i][j - 1].B - cons[s.sol[i][j-1].stop][s.sol[i][j].stop];
        }
        //from here on the charging scheduling procedure really starts
        double r = vehicle[i].r;
        double Q = vehicle[i].Q;
        double tot_charg_needed = (r * Q - s.sol[i][s.sol[i].size()-1].B);        
        if (tot_charg_needed <= 0) { continue; }
        std::vector<ZL> poss(0); //possible nodes after which charging could occur
        for (int j = 0; j < s.sol[i].size()-1; j++) {
            if (s.sol[i][j].C == vehicle[i].C) {
                // Find rk of this node
                tempRK = -1.0;
                if (s.sol[i][j].stop > 2 * nbUsers) { // Charging station or origin depot => no random key. Charging stations won't be picked anyway. For origin depot, set RK = rk of veh
                    tempRK = s.vec[2*nbUsers + i].rk;
                }
                else { // Drop-offs
                    tempRK = s.vec[s.sol[i][j].stop].rk; // Take drop-off gene
                }
                // Add zero-load
                poss.push_back({ j, elecPrices[(s.sol[i][j].ET + node[s.sol[i][j].stop].d) / periodLength], s.sol[i][j+1].LT - node[s.sol[i][j].stop].d - s.sol[i][j].ET, tempRK});
            }
        }
        // Set charging priority
        switch (chargingPriority) {
            case 1: // BY ELECTICITY PRICE
                sort(poss.begin(), poss.end(), sortByEP); //sorting the possible nodes by electricity prices
                break;
            case 2: // EARLIEST FIRST
                // No sort, just by order in route  
                break;
            case 3: // LATEST FIRST
                std::reverse(poss.begin(), poss.end());
                break;
            case 4: // RK BASED
                sort(poss.begin(), poss.end(), sortByRKsched);
                break;
            case 5: // BIGGEST WINDOWS FIRST
                sort(poss.begin(), poss.end(), sortByWindow);
                break;
        }
        double charg_assigned = 0;
        while (charg_assigned < tot_charg_needed) {
            if (poss.size() == 0) {
                s.battery_infeasibles +=1;      //if there still is a need for charging, but there are no nodes available after which charging would be possible, the route/schedule is battery infeasible and the objective function of s will be penalised for this reason.
                break;
            }
            int inspect = poss[0].node;         //node position after which we'll look to add charging
            int stationIndex = -1;
            TVecSol station=findClosestCharger(s.sol[i][inspect], s.sol[i][inspect+1], cStations, dist, node, availability, stationIndex) ; //Find closest fully available charger
            if (station.stop == -1) {
                poss.erase(poss.begin());              // If no available charger found, skip this zero-load
                continue;
            }
            station.ET = 0;
            station.LT = 480;
            station.C = 3;
            station.c_station = true;
            double alpha;
            for (int s = 0; s < cStations.size(); s++) {
                if (station.stop == cStations[s].id - 1) {
                    alpha = cStations[s].alpha;
                    break;
                }
            }


            //now we have a potential station to insert. Let's update the window of poss[0] taking into accounts travel times but also time windows preceding and following the station.
            poss[0].window = poss[0].window - dist[s.sol[i][inspect].stop][station.stop] - dist[station.stop][s.sol[i][inspect + 1].stop];
            double stationET= s.sol[i][inspect].ET + node[s.sol[i][inspect].stop].d + dist[s.sol[i][inspect].stop][station.stop]; //if the station gets inserted after inspect, the ET of charging would be ET of inspect + service time + travel time. The ET of the next node would be ET of station + charging duration + travel time.
            double nextET = stationET + poss[0].window + dist[station.stop][s.sol[i][inspect + 1].stop]; //ET of next node if poss[0].window would be fully charged
            //it could be that a charging duration of current value of poss[0] window would push the ET of nodes following the station past their LT. In such a case, poss[0] window should be made smaller.
            if (nextET> s.sol[i][inspect + 1].LT) {
                poss[0].window = s.sol[i][inspect + 1].LT - dist[station.stop][s.sol[i][inspect + 1].stop] - stationET;
            }
            double stationLT = s.sol[i][inspect + 1].LT - dist[station.stop][s.sol[i][inspect + 1].stop] - poss[0].window; //LT of charging if poss[0].window has to be fully charged
            double precedingLT = stationLT - dist[s.sol[i][inspect].stop][station.stop] - node[s.sol[i][inspect].stop].d;
            //it could be that a charging duration of current value of poss[0] window would push the LT of nodes preceding the station under their ET. In such a case, poss[0] window should be made smaller.
            for (int j=inspect;j>0;j--){
                if (s.sol[i][j].ET>precedingLT) {
                    poss[0].window = poss[0].window - (s.sol[i][inspect].ET - precedingLT);
                }
                precedingLT = precedingLT - dist[s.sol[i][j-1].stop][s.sol[i][j].stop] - node[s.sol[i][j-1].stop].d - s.sol[i][j-1].w;
            }
            if (poss[0].window<= 0) {
                poss.erase(poss.begin());              //if the window is negative, no charging is added and zeroload node is removed from zeroloads list.
                continue;
            }
            if ((poss[0].window) * alpha - cons[s.sol[i][inspect].stop][station.stop] - cons[station.stop][s.sol[i][inspect + 1].stop] + cons[s.sol[i][inspect].stop][s.sol[i][inspect + 1].stop] <= 0) {
                poss.erase(poss.begin());               //if maximum charging at node considering waiting time is lower than battery consumption of detour to station, no charging is added and zeroload node is removed from zeroloads list.
                continue;
            }
            s.sol[i].insert(s.sol[i].begin() + inspect + 1, station);
            tot_charg_needed += (cons[s.sol[i][inspect].stop][s.sol[i][inspect + 1].stop] + cons[s.sol[i][inspect+1].stop][s.sol[i][inspect + 2].stop] - cons[s.sol[i][inspect].stop][s.sol[i][inspect+2].stop]);
            double cons_sum = 0;        //sum of battery consumptions from charging station to destination depot
            for (int l = inspect+1; l < s.sol[i].size()-1; l++) {
                cons_sum += cons[s.sol[i][l].stop][s.sol[i][l+1].stop];
            }
            //the charging added after the currently inspected node is bounded by four time limits: 1) charging time which would allow vehicle to reach destination depot with r*Q battery level from inspected node.   2) time to fully recharge.  3) waiting time - travel time of detour  4) charging time needed over full route not yet assigned
            double a = (r * Q + cons_sum - std::max(s.sol[i][inspect].B - cons[s.sol[i][inspect].stop][s.sol[i][inspect + 1].stop], 0.0) )/ alpha;
            double b = (Q - std::max(s.sol[i][inspect].B - cons[s.sol[i][inspect].stop][s.sol[i][inspect + 1].stop], 0.0) )/ alpha;
            double c = poss[0].window;
            double d = (tot_charg_needed - charg_assigned)/alpha;
            double charge;
            charge = std::min(a, b);
            charge = std::min(charge, c);
            charge = std::min(charge, d);      //this is the charging we'll add after the inspected node (unless charging this would mean a vehicle would have a higher battery level than Q somewhere further along the route)
            for (int l = inspect + 2; l < s.sol[i].size(); l++) {    //check if charging now doesn't imply a battery level higher than Q somewhere further along the line
                if (s.sol[i][l].B + s.sol[i][l].w*alpha + charge * alpha - (cons[s.sol[i][inspect].stop][s.sol[i][inspect + 1].stop] + cons[s.sol[i][inspect + 1].stop][s.sol[i][inspect + 2].stop] - cons[s.sol[i][inspect].stop][s.sol[i][inspect + 2].stop]) > Q) {
                    charge = (Q - (s.sol[i][l].B + s.sol[i][l].w*alpha - (cons[s.sol[i][inspect].stop][s.sol[i][inspect + 1].stop] + cons[s.sol[i][inspect + 1].stop][s.sol[i][inspect + 2].stop] - cons[s.sol[i][inspect].stop][s.sol[i][inspect + 2].stop])) ) / alpha;
                }
            }
            if (charge * alpha <= cons[s.sol[i][inspect].stop][s.sol[i][inspect + 1].stop] + cons[s.sol[i][inspect + 1].stop][s.sol[i][inspect + 2].stop] - cons[s.sol[i][inspect].stop][s.sol[i][inspect + 2].stop]) {     //if charging time ends up not being enough to compensate battery consumption of detour (+ room to charge extra), charging station is removed from route and zeroloads node is removed from zeroloads list.
                tot_charg_needed = tot_charg_needed - (cons[s.sol[i][inspect].stop][s.sol[i][inspect + 1].stop] + cons[s.sol[i][inspect + 1].stop][s.sol[i][inspect + 2].stop] - cons[s.sol[i][inspect].stop][s.sol[i][inspect + 2].stop]);
                s.sol[i].erase(s.sol[i].begin() + inspect + 1);
                poss.erase(poss.begin());
                continue;
            }
            s.sol[i][inspect + 1].ET = s.sol[i][inspect].ET + node[s.sol[i][inspect].stop].d + dist[s.sol[i][inspect].stop][s.sol[i][inspect + 1].stop] ;  //setting earliest time charging can begin
            s.sol[i][inspect+1].LT = s.sol[i][inspect+2].LT - dist[s.sol[i][inspect+1].stop][s.sol[i][inspect+2].stop] - charge;   //setting latest time charging should start so charge time can be respected
            s.sol[i][inspect].LT = std::min(s.sol[i][inspect].LT, s.sol[i][inspect + 1].LT - node[s.sol[i][inspect].stop].d - dist[s.sol[i][inspect].stop][s.sol[i][inspect + 1].stop]);   //LT of preceding drop-off is altered to be able to respect charging
            s.sol[i][inspect + 2].ET = std::max(s.sol[i][inspect + 2].ET, s.sol[i][inspect + 1].ET + charge + dist[s.sol[i][inspect + 1].stop][s.sol[i][inspect + 2].stop]);   //ET of next pick-up is altered to be able to respect charging
            updateAvailability(s.sol[i][inspect + 1].ET, s.sol[i][inspect + 1].LT + charge, stationIndex, cStations, availability); // Reserving whole time window for this vehicle at this station (since we exact charging session is only scheduled at the end)
            for (int j = inspect - 1; j >= 0; j--) { //update LT "backwards" as long as needed
                if (s.sol[i][j + 1].LT - dist[s.sol[i][j].stop][s.sol[i][j + 1] .stop] - node[s.sol[i][j].stop].d -s.sol[i][j].w < s.sol[i][j].LT) { // If LT needs to be updated
                    s.sol[i][j].LT = s.sol[i][j + 1].LT - dist[s.sol[i][j].stop][s.sol[i][j+1].stop] - node[s.sol[i][j].stop].d - s.sol[i][j].w; // Update LT
                }
                else break; // Else, you can stop backwards updating
            }
            for (int j = inspect + 3; j < s.sol[i].size(); j++) {
                if (s.sol[i][j - 1].ET + dist[s.sol[i][j - 1].stop][s.sol[i][j].stop] + node[s.sol[i][j - 1].stop].d + s.sol[i][j-1].w > s.sol[i][j].ET) { // If ET needs to be updated
                    s.sol[i][j].ET = s.sol[i][j - 1].ET + dist[s.sol[i][j - 1].stop][s.sol[i][j].stop] + node[s.sol[i][j - 1].stop].d + s.sol[i][j - 1].w; // Update ET
                }
                else break; // Else, you can stop forwards updating
            }
            s.sol[i][inspect + 1].w = charge; //the charging duration for a charging station is stored in the w variable
            s.sol[i][inspect + 1].B = s.sol[i][inspect].B - cons[s.sol[i][inspect].stop][s.sol[i][inspect + 1].stop]; //setting B at start of charging
            s.sol[i][inspect + 2].B= s.sol[i][inspect+1].B + charge*alpha - cons[s.sol[i][inspect+1].stop][s.sol[i][inspect + 2].stop]; //setting B of node just after charging station
            for (int l = inspect + 3; l < s.sol[i].size(); l++) {  //updating battery levels forward
                s.sol[i][l].B = s.sol[i][l].B + charge*alpha - (cons[s.sol[i][inspect].stop][s.sol[i][inspect + 1].stop] + cons[s.sol[i][inspect + 1].stop][s.sol[i][inspect + 2].stop] - cons[s.sol[i][inspect].stop][s.sol[i][inspect + 2].stop]);
            } 
            charg_assigned += charge*alpha;
            poss.erase(poss.begin());
            for (int p = 0; p < poss.size(); p++) {        //updating node positions in poss if position>inspect because charging station was just inserted in inspect + 1
                if (poss[p].node > inspect) {
                    poss[p].node = poss[p].node + 1;
                }
            }
        }
    }
    if (s.battery_infeasibles > 0) { 
        s.fo = objFct(s, node, dist, nbVehicles, nbUsers,periodLength, maxRideTimes, weights, elecPrices);
        return s; }

    //at this point charging stations have been inserted and charges have been determined. Now the time variables, respecting the time windows, have to be determined.
    // Based on scheduler gene: either naive "late pickup, late charging, early drop-off" or "RK-based pickup, ??? charging, early drop-off

    int numSchedulers = 2;
    int sched = ceil(s.vec[s.vec.size() - 2].rk * numSchedulers + 0.000000000001);

    switch (sched)
    {
    case 1:
        ScheduleLatePUEarlyDO(s, nbVehicles, nbUsers, node, dist, elecPrices, periodLength);
        break;
    case 2:
        ScheduleRkPUEarlyDO(s, nbVehicles, nbUsers, node, dist, elecPrices, periodLength);
        break;
    }

    s.scheduled = true;
    return s;

}

TVecSol findClosestCharger(TVecSol node1, TVecSol node2, std::vector <TCstat> cStations, std::vector <std::vector <double> >& dist, std::vector <TNode>& node, std::vector<std::vector<std::vector<unsigned long long>>>& availability, int & stationIndex) {
    TVecSol station;
    station.stop = -1; // For checking if station has been found
    int extra = 1000000;
    int stop1 = node1.stop;
    int stop2 = node2.stop;
    float start;
    float end;
    for (int s = 0; s < cStations.size(); s++) {
        start = node1.ET + node[node1.stop].d + dist[node1.stop][cStations[s].id-1]; // Soonest time charging could start at this station
        end = node2.LT - dist[cStations[s].id - 1][node2.stop];
        if (dist[stop1][cStations[s].id-1] + dist[cStations[s].id-1][stop2] < extra) { // If closest
            if (checkAvailability(start, end, s, cStations, availability)) { // AND if available
                extra = dist[stop1][cStations[s].id - 1] + dist[cStations[s].id - 1][stop2];
                station.stop = cStations[s].id - 1;
                stationIndex = s;
            }
        }
    }
    return station;
}

bool checkAvailability(float start, float end, int csIndex, std::vector <TCstat> cStations, std::vector<std::vector<std::vector<unsigned long long>>>& availability) {
    unsigned long long session = 0;
    unsigned long long available;

    // For each hour
    for (int hour = 0; hour < availability[csIndex][0].size(); hour++) {
        session = 0;
        if (start > (hour + 1) * 60) { // If session happens after this hour
            continue;
        }
        else if (start > hour * 60 && start <= (hour + 1) * 60) { // If start during this hour
            if (end > hour * 60 && end <= (hour + 1) * 60) { // If end during hour
                session = (unsigned long long)pow(2, (int)floor(end - hour * 60)) - 1 ^ (unsigned long long) pow(2, (int)floor(start - hour * 60)) - 1; // Set bitwise session
            }
            else { // Departure must be after hour if not during
                session = (unsigned long long)pow(2, 60) - 1 ^ (unsigned long long) pow(2, (int)ceil(start - hour * 60)) - 1; // Set bitwise session
            }
        }
        else { // Arrival happens before this hour
            if (end <= hour * 60) { // If departure also before this hour
                break; // Then you have already checked whole session
            }
            else if (end > hour * 60 && end <= (hour + 1) * 60) { // If departure during hour
                session = (unsigned long long)pow(2, (int)floor(end - hour * 60)) - 1; // Set bitwise session
            }
            else { // Departure must be after hour if not during
                session = (unsigned long long)pow(2, 60) - 1; // Set bitwise session
            }
        }

        // Check session against availability
        available = availability[csIndex][cStations[csIndex].cap-1][hour]; // This is the "charging plug" with the highest index. 
        //Since charging plugs with lowest index are always "filled" first, the one with the highest index only has 1's if charging station is at full capacity

        if (available & session) { // If there is at least 1 minute during the session at which capacity would be exceeded => Not allowed
            return false;
        }
    }

    return true;
}

void updateAvailability(float start, float end, int csIndex, std::vector <TCstat> cStations, std::vector<std::vector<std::vector<unsigned long long>>>& availability) {
    unsigned long long session = 0;
    unsigned long long auxULL = 0;

    // For each hour
    for (int hour = 0; hour < availability[csIndex][0].size(); hour++) {
        session = 0;
        if (start > (hour + 1) * 60) { // If session happens after this hour
            continue;
        }
        else if (start > hour * 60 && start <= (hour + 1) * 60) { // If arrival during this hour
            if (end > hour * 60 && end <= (hour + 1) * 60) { // If departure during hour
                session = (unsigned long long)pow(2, (int)floor(end - hour * 60)) - 1 ^ (unsigned long long) pow(2, (int)ceil(start - hour * 60)) - 1; // Set bitwise session
            }
            else { // Departure must be after hour if not during
                session = (unsigned long long)pow(2, 60) - 1 ^ (unsigned long long) pow(2, (int)ceil(start - hour * 60)) - 1; // Set bitwise session
            }
        }
        else { // Arrival happens before this hour
            if (end <= hour * 60) { // If departure also before this hour
                break; // Then you have already checked whole session
            }
            else if (end > hour * 60 && end <= (hour + 1) * 60) { // If departure during hour
                session = (unsigned long long)pow(2, (int)floor(end - hour * 60)) - 1; // Set bitwise session
            }
            else { // Departure must be after hour if not during
                session = (unsigned long long)pow(2, 60) - 1; // Set bitwise session
            }
        }

        // Update availability
        for (int c = 0; c < cStations[csIndex].cap; c++) {
            auxULL = availability[csIndex][c][hour] & session;
            availability[csIndex][c][hour] = availability[csIndex][c][hour] | session;
            session = auxULL; // Whaatever could not be alloted to charger is carried over to allot to next charger
        }
        // Now session should == 0
        if (session) {
            int debug = 0;
        }
    }
}

void ScheduleLatePUEarlyDO(TSol& s, int nbVehicles, int nbUsers, std::vector <TNode>& node, std::vector <std::vector <double> >& dist, std::vector <double> elecPrices, int periodLength) {
    double priceET;
    double priceLT;
    int ETPeriod;
    int LTPeriod;

    for (int i = 0; i < nbVehicles; i++) {
        s.sol[i][0].T = 0;
        for (int j = 1; j < s.sol[i].size(); j++) {
            if (s.sol[i][j].stop < nbUsers) { // Pick-up
                s.sol[i][j].T = s.sol[i][j].LT;
            }
            else if (s.sol[i][j].c_station) { // Charging station
                ETPeriod = floor(s.sol[i][j].ET / periodLength);
                LTPeriod = floor((s.sol[i][j].LT + s.sol[i][j].w) / periodLength);
                if (ETPeriod == elecPrices.size()) { ETPeriod--; } // Correction for edgecase where ET == time horizon
                if (LTPeriod == elecPrices.size()) { LTPeriod--; } // Correction for edgecase where LT+w == time horizon
                priceET = elecPrices[ETPeriod];
                priceLT = elecPrices[LTPeriod];
                if (priceLT <= priceET) { // If late is cheaper
                    s.sol[i][j].T = s.sol[i][j].LT;
                }
                else if (priceET < priceLT) { // If early is cheaper
                    s.sol[i][j].T = std::max(s.sol[i][j].ET, s.sol[i][j - 1].T + node[s.sol[i][j - 1].stop].d + dist[s.sol[i][j - 1].stop][s.sol[i][j].stop] + s.sol[i][j - 1].w);
                }
            }
            else { // Drop-off (or final depot)
                s.sol[i][j].T = std::max(s.sol[i][j].ET, s.sol[i][j - 1].T + node[s.sol[i][j - 1].stop].d +dist[s.sol[i][j - 1].stop][s.sol[i][j].stop] + s.sol[i][j-1].w);
            }
        }
    }
}

void ScheduleRkPUEarlyDO(TSol& s, int nbVehicles, int nbUsers, std::vector <TNode>& node, std::vector <std::vector <double> >& dist, std::vector <double> elecPrices, int periodLength) {

    float RK;
    double priceET;
    double priceLT;
    int ETPeriod;
    int LTPeriod;
    int precision = 1000000;

    for (int i = 0; i < nbVehicles; i++) {
        s.sol[i][0].T = 0;
        for (int j = 1; j < s.sol[i].size(); j++) {
            // Pick-up
            if (s.sol[i][j].stop < nbUsers) {
                if (s.sol[i][j - 1].c_station) { // If previous node is a charging station
                    s.sol[i][j].T = s.sol[i][j].LT; // Schedule as late as possible
                }
                else { // Else, RK-based
                    RK = s.vec[s.sol[i][j].stop].rk; // Get true random key
                    RK = (float)((int)(RK * precision) % (precision / nbVehicles)) / (float)precision; // "Filtered" random key
                    RK =  RK * (float)nbVehicles; // Scaled filtered random-key
                    s.sol[i][j].T = std::max(
                        s.sol[i][j].ET + RK * (s.sol[i][j].LT - s.sol[i][j].ET), // T based on RK
                        s.sol[i][j - 1].T + node[s.sol[i][j - 1].stop].d + dist[s.sol[i][j - 1].stop][s.sol[i][j].stop] + s.sol[i][j - 1].w // Earliest possible arrival, taking start-time + service time of previous stop into account
                    );
                    if (s.sol[i][j].T > node[s.sol[i][j].stop].dep || RK > 1) {
                        int debug = 9;
                    }
                }
            }
            else if (s.sol[i][j].c_station) { // Charging station
                ETPeriod = floor(s.sol[i][j].ET / periodLength);
                LTPeriod = floor((s.sol[i][j].LT + s.sol[i][j].w) / periodLength);
                if (ETPeriod == elecPrices.size()) { ETPeriod--; } // Correction for edgecase where ET == time horizon
                if (LTPeriod == elecPrices.size()) { LTPeriod--; } // Correction for edgecase where LT+w == time horizon
                priceET = elecPrices[ETPeriod];
                priceLT = elecPrices[LTPeriod];
                if (priceLT <= priceET) { // If late is cheaper
                    s.sol[i][j].T = s.sol[i][j].LT;
                }
                else if (priceET < priceLT) { // If early is cheaper
                    s.sol[i][j].T = std::max(s.sol[i][j].ET, s.sol[i][j - 1].T + node[s.sol[i][j - 1].stop].d + dist[s.sol[i][j - 1].stop][s.sol[i][j].stop] + s.sol[i][j - 1].w);
                }
            }
            else { // Drop-off (or final depot)
                s.sol[i][j].T = std::max(s.sol[i][j].ET, s.sol[i][j - 1].T + node[s.sol[i][j - 1].stop].d +dist[s.sol[i][j - 1].stop][s.sol[i][j].stop] + s.sol[i][j-1].w); // Earliest between ET and earliest possible arrival given T, d and dist
            }
        }
    }
}