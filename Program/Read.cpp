#include "Read.h"

// FOR TEST PURPOSES
#include <iomanip>
#include <iostream>

void ReadData(char nameTable[], int& n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl>& vehicle,
    int& nbVehicles, int& nbUsers, int& periodLength, int& H, std::vector <int>& oDepotIDs, std::vector <int>& fDepotIDs,
    std::vector <TCstat>& cStations, std::vector <int>& maxRideTimes, double(&weights)[3], std::vector <double>& elecPrices, std::vector <std::vector <double> >& cons)
{
    switch (nameTable[12]) {
    case 'a':
        ReadDataCordeau(nameTable, n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons);
        break;
    case 'u':
        ReadDataUber(nameTable, n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons);
        break;
    }

    // Time-window tightening
    for (int i = 0; i < nbUsers; i++) {
        node[i].arr = std::max(node[i].arr, node[i + nbUsers].arr - maxRideTimes[i] - node[i].d);
        node[i].dep = std::min(node[i].dep, node[i + nbUsers].dep - dist[i][i + nbUsers] - node[i].d);
    }
    for (int i = nbUsers; i < 2 * nbUsers; i++) {
        node[i].arr = std::max(node[i].arr, node[i - nbUsers].arr + dist[i - nbUsers][i] + node[i - nbUsers].d);
        node[i].dep = std::min(node[i].dep, node[i - nbUsers].dep + node[i - nbUsers].d + maxRideTimes[i - nbUsers]);
    }

    // Arc deletion
    // TO BE SEEN

    // TEST
    
    /*printf("dist\n");
    for (std::vector<double> i:dist){
        for (double j : i) {
            std::cout << std::setfill(' ') << std::setw(8) << j;
        }
        printf("\n");
    }
    printf("\n\n");

    printf("cons\n");
    for (std::vector<double> i : cons) {
        for (double j : i) {
            std::cout << std::setfill(' ') << std::setw(5) << std::fixed << std::setprecision(2)<< j;
        }
        printf("\n");
    }
    printf("\n\n");
    printf("node\n");
    for (TNode n : node) {
        std::cout << n.id << " " << n.x << " "  << n.y << " "  << n.d << " "  << n.l << " "  << n.arr << " "  << n.dep;
        printf("\n");
    }
    printf("\n\n");

    printf("nbVehicles: ");
    std::cout << nbVehicles << "\n";

    printf("nbUsers: ");
    std::cout << nbUsers << "\n";

    printf("periodLength: ");
    std::cout << periodLength << "\n";

    printf("Planning horizon H: ");
    std::cout << H << "\n";

    printf("Origin depot IDs: ");
    for (int id : oDepotIDs) { std::cout << id << " "; };
    printf("\n");

    printf("Final depot IDs: ");
    for (int id : fDepotIDs) { std::cout << id << " "; };
    printf("\n");

    printf("Charging stations: \n");
    for (TCstat cs: cStations) { std::cout << cs.id << " " << cs.alpha << " " << cs.cap << "\n"; };
    printf("\n");

    printf("Max ride times: ");
    for (int mrt : maxRideTimes) { std::cout << mrt << " "; };
    printf("\n");

    printf("Weights: ");
    for (double w : weights) { std::cout << w << " "; };
    printf("\n");

    printf("Electricity prices: ");
    for (double ep : elecPrices) { std::cout << ep << " "; };
    printf("\n");*/
}

void ReadDataCordeau(char nameTable[], int& n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl>& vehicle,
    int& nbVehicles, int& nbUsers, int& periodLength, int& H, std::vector <int>& oDepotIDs, std::vector <int>& fDepotIDs,
    std::vector <TCstat>& cStations, std::vector <int>& maxRideTimes, double(&weights)[3], std::vector <double>& elecPrices, std::vector <std::vector <double> >& cons)
{
    char name[200] = "Instances";
    strcat(name, nameTable);

    FILE* arq;
    arq = fopen(name, "r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n", name);
        getchar();
        exit(1);
    }

    // => read data

    // read instance head
    int nbODepots;
    int nbFDepots;
    int nbCStations;
    int nbRep; // Irrelevant number
    fscanf(arq, "%d %d %d %d %d %d %d %d", &nbVehicles, &nbUsers, &nbODepots, &nbFDepots, &nbCStations, &nbRep, &periodLength, &H);

    // read node informations
    int nAux = 2 * nbUsers + 2 * nbVehicles + nbCStations + 2; // !!!!! NOT GENERALISABLE; BASED ON CORDEAU INSTANCES ONLY
    node.clear();
    TNode nodeTemp;

    for (int i = 0; i < nAux; i++)
    {
        fscanf(arq, "%d %lf %lf %f %d %lf %lf", &nodeTemp.id, &nodeTemp.x, &nodeTemp.y, &nodeTemp.d, &nodeTemp.l, &nodeTemp.arr, &nodeTemp.dep);
        node.push_back(nodeTemp);
    }

    // read common origin depot id
    int trash;
    for (int i = 0; i < nbODepots; i++) {
        fscanf(arq, "%d", &trash); // Useless ID's, only artificial depots are used
    }

    // read common destination depot id
    for (int i = 0; i < nbFDepots; i++) {
        fscanf(arq, "%d", &trash); // Useless ID's, only artificial depots are used
    }

    // read (artificial) origin depots id
    int id_temp = -1;
    for (int i = 0; i < nbVehicles; i++) {
        fscanf(arq, "%d", &id_temp);
        oDepotIDs.push_back(id_temp);
    }

    // read (artificial) destination depots id
    id_temp = -1;
    for (int i = 0; i < nbVehicles; i++) {
        fscanf(arq, "%d", &id_temp);
        fDepotIDs.push_back(id_temp);
    }

    // read charging stations id
    TCstat CS_temp;
    for (int i = 0; i < nbCStations; i++) {
        fscanf(arq, "%d", &CS_temp.id);
        cStations.push_back(CS_temp); // Charging station with only ID is pushed. Alpha and cap are read/written later.
    }

    // read users maximum ride time
    int tempRT;
    for (int i = 0; i < nbUsers; i++) {
        fscanf(arq, "%d", &tempRT);
        maxRideTimes.push_back(tempRT);
    }

    // Initialize vehicle
    TVhcl vhclTemp;
    for (int i = 1; i <= nbVehicles; i++) {
        vhclTemp.id = i;
        vehicle.push_back(vhclTemp);
    }

    // read vehicles capacity
    for (int i = 0; i < nbVehicles; i++) {
        fscanf(arq, "%d", &vehicle[i].C);
    }

    // read vehicles initial battery inventory
    for (int i = 0; i < nbVehicles; i++) {
        fscanf(arq, "%lf", &vehicle[i].B0);
    }

    // read vehicles battery capacities
    for (int i = 0; i < nbVehicles; i++) {
        fscanf(arq, "%lf", &vehicle[i].Q);
    }

    // read minimum end battery ratio levels
    for (int i = 0; i < nbVehicles; i++) {
        fscanf(arq, "%lf", &vehicle[i].r);
    }

    // read recharging rates at charging stations
    for (int i = 0; i < nbCStations; i++) {
        fscanf(arq, "%lf", &cStations[i].alpha);
    }

    // read vehicles discharging rate
    double dischRate;
    fscanf(arq, "%lf", &dischRate);

    // read weight factors
    fscanf(arq, "%lf %lf %lf", &weights[0], &weights[1], &weights[2]);

    // read capacities of the charging stations
    for (int i = 0; i < nbCStations; i++) {
        fscanf(arq, "%int", &cStations[i].cap);
    }

    // read period electricity prices
    int nbPeriods = ceil((float)H / periodLength);
    double tempP;
    for (int i = 0; i < nbPeriods; i++) {
        fscanf(arq, "%lf", &tempP);
        tempP = std::max(0.0, tempP);
        elecPrices.push_back(tempP);
    }

    fclose(arq);

    // calculate the euclidean distance
    dist.clear();
    dist.resize(nAux, std::vector<double>(nAux));

    for (int i = 0; i < nAux; i++)
    {
        for (int j = i; j < nAux; j++)
        {
            dist[i][j] = dist[j][i] = sqrt((node[j].x - node[i].x) * (node[j].x - node[i].x) + (node[j].y - node[i].y) * (node[j].y - node[i].y));
        }
    }


    // Calculate the battery consumption
    cons.clear();
    cons.resize(nAux, std::vector<double>(nAux));

    for (int i = 0; i < nAux; i++)
    {
        for (int j = i; j < nAux; j++)
        {
            cons[i][j] = cons[j][i] = dischRate * dist[i][j];
        }
    }

    n = 2* nbUsers + nbVehicles + 2; // Chromosome length. 1 per pickup, 1 per dropoff, 1 for scheduler, 1 for CP, Should exclude chromosome for decoder.
} 

void ReadDataUber(char nameTable[], int& n, std::vector <TNode>& node,
    std::vector <std::vector <double> >& dist, std::vector <TVhcl>& vehicle,
    int& nbVehicles, int& nbUsers, int& periodLength, int& H, std::vector <int>& oDepotIDs, std::vector <int>& fDepotIDs,
    std::vector <TCstat>& cStations, std::vector <int>& maxRideTimes, double(&weights)[3], std::vector <double>& elecPrices, std::vector <std::vector <double> >& cons)
{
    char name[200] = "Instances";
    strcat(name, nameTable);

    FILE* arq;
    arq = fopen(name, "r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n", name);
        getchar();
        exit(1);
    }

    // => read data

    // read instance head
    int nbODepots;
    int nbFDepots;
    int nbCStations;
    int nbRep; // Irrelevant number
    fscanf(arq, "%d %d %d %d %d %d %d %d", &nbVehicles, &nbUsers, &nbODepots, &nbFDepots, &nbCStations, &nbRep, &periodLength, &H);

    // read node informations
    int nAux = 2 * nbUsers + nbVehicles + 5 + nbCStations + 2; // Nbvehicles duplicates of origin depot & 5 duplicates of final depot
    node.clear();
    TNode nodeTemp;
    float auxL;

    for (int i = 0; i < nAux; i++)
    {
        fscanf(arq, "%d %lf %lf %f %f %lf %lf", &nodeTemp.id, &nodeTemp.x, &nodeTemp.y, &nodeTemp.d, &auxL, &nodeTemp.arr, &nodeTemp.dep);
        nodeTemp.l = (int)auxL;
        node.push_back(nodeTemp);
    }

    // read common origin depot id
    int trash;
    for (int i = 0; i < nbODepots; i++) {
        fscanf(arq, "%d", &trash); // Useless ID's, only artificial depots are used
    }

    // read common destination depot id
    for (int i = 0; i < nbFDepots; i++) {
        fscanf(arq, "%d", &trash); // Useless ID's, only artificial depots are used
    }

    // read (artificial) origin depots id
    int id_temp = -1;
    for (int i = 0; i < nbVehicles; i++) {
        fscanf(arq, "%d", &id_temp);
        oDepotIDs.push_back(id_temp);
    }

    // read (artificial) destination depots id
    id_temp = -1;
    for (int i = 0; i < 5; i++) {
        fscanf(arq, "%d", &id_temp);
        fDepotIDs.push_back(id_temp);
    }

    // read charging stations id
    TCstat CS_temp;
    for (int i = 0; i < nbCStations; i++) {
        fscanf(arq, "%d", &CS_temp.id);
        cStations.push_back(CS_temp); // Charging station with only ID is pushed. Alpha and cap are read/written later.
    }

    // read users maximum ride time
    int tempRT;
    for (int i = 0; i < nbUsers; i++) {
        fscanf(arq, "%d", &tempRT);
        maxRideTimes.push_back(tempRT);
    }

    // Initialize vehicle
    TVhcl vhclTemp;
    for (int i = 1; i <= nbVehicles; i++) {
        vhclTemp.id = i;
        vehicle.push_back(vhclTemp);
    }

    // read vehicles capacity
    for (int i = 0; i < nbVehicles; i++) {
        fscanf(arq, "%d", &vehicle[i].C);
    }

    // read vehicles initial battery inventory
    for (int i = 0; i < nbVehicles; i++) {
        fscanf(arq, "%lf", &vehicle[i].B0);
    }

    // read vehicles battery capacities
    for (int i = 0; i < nbVehicles; i++) {
        fscanf(arq, "%lf", &vehicle[i].Q);
    }

    // read minimum end battery ratio levels
    for (int i = 0; i < nbVehicles; i++) {
        fscanf(arq, "%lf", &vehicle[i].r);
    }

    // read recharging rates at charging stations
    for (int i = 0; i < nbCStations; i++) {
        fscanf(arq, "%lf", &cStations[i].alpha);
    }

    // read vehicles discharging rate
    double dischRate;
    fscanf(arq, "%lf", &dischRate);

    // read weight factors
    fscanf(arq, "%lf %lf %lf", &weights[0], &weights[1], &weights[2]);

    // Read travel times
    dist.clear();
    dist.resize(nAux, std::vector<double>(nAux));
    double aux;
    for (int i = 0; i < nAux; i++)
    {
        for (int j = 0; j < nAux; j++)
        {
            fscanf(arq, "%lf", &aux);
            dist[i][j] = 2*aux;
        }
    }

    // read capacities of the charging stations
    for (int i = 0; i < nbCStations; i++) {
        fscanf(arq, "%int", &cStations[i].cap);
    }

    // read period electricity prices
    int nbPeriods = ceil((float)H / periodLength);
    double tempP;
    for (int i = 0; i < nbPeriods; i++) {
        fscanf(arq, "%lf", &tempP);
        tempP = std::max(0.0, tempP);
        elecPrices.push_back(tempP);
    }

    fclose(arq);


    // Calculate the battery consumption
    cons.clear();
    cons.resize(nAux, std::vector<double>(nAux));

    for (int i = 0; i < nAux; i++)
    {
        for (int j = 0; j < nAux; j++)
        {
            cons[i][j] = dischRate * dist[i][j];
        }
    }

    n = 2 * nbUsers + nbVehicles +2; // Chromosome length. 1 per pickup, 1 per dropoff, 1 for scheduler, 1 for CP, Should exclude chromosome for decoder.
}

void FreeMemoryProblem(std::vector <TNode> &node, std::vector <std::vector <double> > &dist, std::vector <TVhcl> &vehicle, std::vector <int> &oDepotIDs, std::vector <int>&fDepotIDs,
    std::vector <TCstat>&cStations, std::vector <int>&maxRideTimes, std::vector <double>&elecPrices, std::vector <std::vector <double> >&cons)
{
    //specific problem
    dist.clear();
    node.clear();
    vehicle.clear();
    oDepotIDs.clear();
    fDepotIDs.clear();
    cStations.clear();
    maxRideTimes.clear();
    elecPrices.clear();
    cons.clear();
}
