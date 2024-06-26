#include "BRKGA_QL.h"

/************************************************************************************
				MAIN FUNCTION AREA
*************************************************************************************/
int main()
{
	// file with test instances
	FILE *arqProblems;
	arqProblems = fopen ("arqProblems.csv", "r"); 

	if (arqProblems == NULL)
	{
		printf("\nERROR: File arqProblems.txt not found\n");
		getchar();
		exit(1);
	}

	char nameTable[100];
	fgets(nameTable, sizeof(nameTable), arqProblems); //read first line of arqProblems file

	// best solution that is saved in out file
	TSol sBest;

	// run the A-BRKGA for all test instances
	while (!feof(arqProblems))
	{
		// read the name of instances, debug mode, local search module, maximum time, maximum number of runs, maximum number of threads
		fscanf(arqProblems,"%s %d %d %d %d %d", nameTable, &debug, &ls, &MAXTIME, &MAXRUNS, &MAX_THREADS);
		strcpy(instance,nameTable);
        
		//read the informations of the instance
		ReadData(nameTable, n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons);

		double foBest = INFINITY,
		       foAverage = 0;

		float timeBest = 0,
		      timeTotal = 0;

		std::vector <double> fos;
		fos.clear();

        std::vector <bool> lss;
        lss.clear();

        std::vector <float> decs;
        decs.clear();

        std::vector <float> scheds;
        scheds.clear();

        std::vector <float> cps;
        cps.clear();

		// best solutions found in MAXRUNS
		sBest.fo = INFINITY;

		// run ABRKGA MaxRuns for each instance
		printf("\n\nInstance: %s \nRun: ", instance);
		for (int j=0; j<MAXRUNS; j++)
		{
		    // fixed seed
		    srand(j+1); 
		    //srand(time(NULL));

		    printf("%d ", j+1);
		    CPUbegin = CPUend = CPUbest = clock();
		    //gettimeofday(&Tstart, NULL);
		    //gettimeofday(&Tend, NULL);
		    //gettimeofday(&Tbest, NULL);

		    // best solution found in this run
		    bestSolution.fo = INFINITY;

		    // execute the evolutionary method
		    BRKGA();

		    CPUend = clock();
		    //gettimeofday(&Tend, NULL);

            printf(" %lf - %lf\n", (float)CPUbegin/CLOCKS_PER_SEC, (float)CPUend/CLOCKS_PER_SEC);

		    // output results
		    if (bestSolution.fo < sBest.fo)
			sBest = bestSolution;

		    // calculate best and average results
		    if (bestSolution.fo < foBest)
			foBest = bestSolution.fo;

		    foAverage += bestSolution.fo;

		    // fitness of each solution found in the runs
		    fos.push_back(bestSolution.fo);

            // Local search usage of each solution found in the runs
            lss.push_back(bestSolution.flag);

            // Decoder usage of each solution found in the runs
            decs.push_back(bestSolution.vec[n].rk);

            // Scheduling usage of each solution found in the runs
            scheds.push_back(bestSolution.vec[n - 1].rk);

            // Charging priority usage of each solution found in the runs
            cps.push_back(bestSolution.vec[n - 2].rk);


		    timeBest += (float)(CPUbest - CPUbegin)/CLOCKS_PER_SEC;
		    timeTotal += (float)(CPUend - CPUbegin)/CLOCKS_PER_SEC;
		    //timeBest += ((Tbest.tv_sec  - Tstart.tv_sec) * 1000000u + Tbest.tv_usec - Tstart.tv_usec) / 1.e6;
		    //timeTotal += ((Tend.tv_sec  - Tstart.tv_sec) * 1000000u + Tend.tv_usec - Tstart.tv_usec) / 1.e6; 
		}

		// create a .xls file with average results
		foAverage = foAverage / MAXRUNS;
		timeBest = timeBest / MAXRUNS;
		timeTotal = timeTotal / MAXRUNS;

        // Calculate TRT, ERT and CC seperately for best solution
        splitF0(sBest, node,dist,nbVehicles,nbUsers, periodLength, maxRideTimes, weights, elecPrices);

		if (!debug)
		{
			WriteSolution(sBest, n, timeBest, timeTotal, instance);
			WriteResults(foBest, foAverage, fos, lss, decs, scheds, cps, timeBest, timeTotal, instance, ls, sBest.TRT, sBest.ERT, sBest.CC);
		}
		else
		{
		    WriteSolutionScreen(sBest, n, timeBest, timeTotal, instance, nbVehicles);
		}

		// free memory of problem variables
		FreeMemory();
		FreeMemoryProblem(node, dist, vehicle, oDepotIDs, fDepotIDs,cStations, maxRideTimes, elecPrices, cons);
		//FreeMemoryProblem();
	}

    	fclose(arqProblems);
	return 0;
}


/************************************************************************************
			                  GENERAL FUNCTIONS
*************************************************************************************/

void BRKGA()
{
    // initialize Q-Table
    InitiateQTable();
    
    // initialize population
    Pop.clear();
    PopInter.clear();

    // population size
    p = sizeP[sizeof(sizeP)/sizeof(sizeP[0]) - 1]; //higher population size

    Pop.resize(p);
    PopInter.resize(p);

    // Create the initial chromossoms with random keys
    #pragma omp parallel for num_threads(MAX_THREADS)
    for (int i=0; i<p; i++)
    {
        TSol ind = CreateInitialSolutions(); 
        ind = Decoder(ind, n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons);
        Pop[i] = PopInter[i] = ind;

        // save the best solution found in this run
        updateBestSolution(Pop[i]);
    }
    
    // sort population in increase order of fitness
    sort(Pop.begin(), Pop.end(), sortByFitness);

    
    int numGenerations = 0;             // number of generations
    int bestGeneration = 0;             // generation in which found the best solution
    double bestFitness = Pop[0].fo;     // best fitness found in past generations
    double averageOffspring = 0;        // average fitness of offsprings
    double lastAvgOffspring = 0;        // last average fitness of offsprings
    double bestOffspring = 0;           // best offspring of each generation
    float currentTime = 0;              // computational time of the search process
    int numLS = 0;                      // number of local search applied in a generation

    // run the evolutionary process until stop criterion
    while(1)
    {
    	// number of generations
        numGenerations++;

        // choose a action for each parameter and update its value
        ChooseAction(numGenerations);

        // define population size
        if (Pop.size() > p)
        {
            Pop.resize(p);
            PopInter.resize(p);
        }
        else if (Pop.size() < p)
        {
            int currentP = Pop.size();

            Pop.resize(p);
            PopInter.resize(p);
            
            for (int k = currentP; k < p; k++)
            {
                TSol ind = CreateInitialSolutions();
                ind = Decoder(ind,n,node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons);
                Pop[k] = PopInter[k] = ind;

                // save the best solution found in this run
                updateBestSolution(Pop[k]);
            }
        }

        // The 'pe' best chromosomes are maintained, so we just copy these into PopInter:
        #pragma omp parallel for num_threads(MAX_THREADS)
        for (int i=0; i<(int)(p*pe); i++)
        {
            // copy the chromosome for next generation
            PopInter[i] = Pop[i];
        }

        // We'll mate 'p - pe - pm' pairs; initially, i = pe, so we need to iterate until i < p - pm:
        averageOffspring = 0;
        bestOffspring = INFINITY;

        // #pragma omp parallel for num_threads(MAX_THREADS) ** posso ter problemas devido a todas threads acessar a mesma variavel
        for (int i = (int)(p*pe); i < p - (int)(p*pm); i++)
        {
            // Parametric uniform crossover            
            PopInter[i] = ParametricUniformCrossover((int)(p*pe));
            PopInter[i] = Decoder(PopInter[i], n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons);

            // save the best solution found in this run
            updateBestSolution(PopInter[i]);

            averageOffspring += PopInter[i].fo;

            if (PopInter[i].fo < bestOffspring)
                bestOffspring = PopInter[i].fo;
        }

        // update average offspring fitness of last generation
        if (numGenerations > 1){
            lastAvgOffspring = averageOffspring;
        }

        // calculate average offspring fitness of current generation
        averageOffspring = (averageOffspring)/(p - (int)(p*pe) - (int)(p*pm));

        if (numGenerations == 1){
            lastAvgOffspring = averageOffspring;
        }
        
        // We'll introduce 'pm' mutants:
        //#pragma omp parallel for num_threads(MAX_THREADS)
        for (int i = p - (int)(p*pm) - (int)(p*pe); i < p; i++)
        {
            PopInter[i] = CreateInitialSolutions();
            PopInter[i] = Decoder(PopInter[i], n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons);

            // save the best solution found in this run
            updateBestSolution(PopInter[i]);
        }  
        
        // Update the current population
        Pop = PopInter;   

        // Sort population in increase order of fitness
        sort(Pop.begin(), Pop.end(), sortByFitness);

        // Verify if we improve the best fitness found until now and set the reward
        double gap = (bestOffspring - bestFitness)/bestFitness;
        if (gap < 0)
        {
            bestFitness = bestOffspring;
            bestGeneration = numGenerations;
            R = 1 + abs(gap);

            stagnation = 0;
        }
        else
        {
            stagnation++; 
            R = 0; 
        }
        
        // Update the values of Q-Table
        UpdateQTable();

        // ********************* LOCAL SEARCH IN COMMUNITIES *******************
        if (ls)
        {
            // number of local search runs
	    numLS = 0;
            
            //apply local search when BRKGA found a new better solution or when the BRKGA don't found a better solution in T generations
	    if (R >= 1 || stagnation > 5) 
	    {
                stagnation = 0;

	            // Find commuties in Top chromossoms
	            IC((int)(p*pe));

	            std::vector <int> promisingSol; 
                    promisingSol.clear();

	            for (int i=0; i < (int)(p*pe); i++)
	            {
	                if (Pop[i].promising == 1)
	                {
	                	promisingSol.push_back(i);
	                }
	            }

	            #pragma omp parallel for num_threads(MAX_THREADS)
		    for (unsigned int i=0; i < promisingSol.size(); i++)
		    {
			    // local search not influence the evolutionary process
			    TSol s = LocalSearch(Pop[promisingSol[i]], n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons);
                s.flag = 1; // To know afterwards whether this solution came from LS or not
                //if (s.fo < Pop[promisingSol[i]].fo) { nbTimesLSImprovedSol++; } // FOR TESTING PURPOSES
			    updateBestSolution(s);
			    numLS++;

			    // set flag as 1 to prevent new local search in the same solution
			    Pop[promisingSol[i]].flag = 1;
		     }
		     promisingSol.clear();
	        }
	    }

        // ******************************* RESTART *****************************
        if ((numGenerations - bestGeneration) > 40)
        {
            bestGeneration = numGenerations;
            #pragma omp parallel for num_threads(MAX_THREADS)
            for (int i=0; i<p; i++)
            {
                TSol ind = CreateInitialSolutions(); 
                ind = Decoder(ind, n, node, dist, vehicle, nbVehicles, nbUsers, periodLength, H, oDepotIDs, fDepotIDs, cStations, maxRideTimes, weights, elecPrices, cons);
                Pop[i] = ind;

                // save the best solution found in this run
                updateBestSolution(Pop[i]);
            }
            sort(Pop.begin(), Pop.end(), sortByFitness);
            bestFitness = Pop[0].fo;

            if (debug)
                printf("\n\nRestart...\n\n");
        }


        // print screen
        if (debug){
            printf("\nGeneration: %3d [%4d - %3d(%.2lf) (%3d) - %3d(%.2lf) - (%.2lf)] \t %.2lf  \t %.2lf [%.4lf] \t %.2lf",
                      numGenerations, p, (int)(p*pe), pe, numLS, (int)(p*pm), pm, rhoe, bestSolution.fo, bestFitness, R, averageOffspring);
            // p = Population size
            // pe = Percentage Elite chromosomes (copied to next generation)
            // pm = percentage of population to be replaced by mutants
            // rhoe = Probability that offspring inherit an allele from elite parent
            // R = reward (for Q-table?)
        }

        // terminate the evolutionary process in MAXTIME
        //gettimeofday(&Tend, NULL);
        //currentTime = ((Tend.tv_sec  - Tstart.tv_sec) * 1000000u + Tend.tv_usec - Tstart.tv_usec) / 1.e6; 

        CPUend = clock();
        currentTime = (float)(CPUend - CPUbegin)/CLOCKS_PER_SEC;
        
        // stop criterium
        if (currentTime >= MAXTIME){
            break;
        }
    }

    // print Q-Table
    if(debug) printf("\nQ-Table:");
    for (int q=0; q<Q.size() && debug; q++)
    {
        printf("\n");
        for (int j=0; j<Q[q].size(); j++)
        {
            if (q>0)
                printf("[%d %.2lf %.3lf %d] \t ", Q[q][j].S, Q[q][j].pVar, Q[q][j].q, Q[q][j].k);
            else
                printf("[%d %.0lf %.3lf %d] \t ", Q[q][j].S, Q[q][j].pVar, Q[q][j].q, Q[q][j].k);
        }
        printf("\n");
    }

    // free memory
    Pop.clear();
    PopInter.clear();
    Q.clear();
}

void updateBestSolution(TSol s)
{
    // save the best solution found in this run
    if (s.fo < bestSolution.fo)
    {
        bestSolution = s;
        CPUbest = clock();
        //gettimeofday(&Tbest, NULL);
    }
}

void InitiateQTable()
{
    Q.clear();
    Q.resize(par);

    qTotal = 0.0;

    // initiate Q-Table with random values
    for (int j=0; j<sizeof(sizeP)/sizeof(sizeP[0]); j++)
    {
        TQ aux;
        aux.S = 0;
        aux.pVar = sizeP[j];
        aux.q = 0; 
        aux.k = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }

    for (int j=0; j<sizeof(Pe)/sizeof(Pe[0]); j++)
    {
        TQ aux;
        aux.S = 1;
        aux.pVar = Pe[j];
        aux.q = 0; //randomico(1.0,1.0);
        aux.k = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }

    for (int j=0; j<sizeof(Pm)/sizeof(Pm[0]); j++)
    {
        TQ aux;
        aux.S = 2;
        aux.pVar = Pm[j];
        aux.q = 0; 
        aux.k = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }

    for (int j=0; j<sizeof(Rhoe)/sizeof(Rhoe[0]); j++)
    {
        TQ aux;
        aux.S = 3;
        aux.pVar = Rhoe[j];
        aux.q = 0; 
        aux.k = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }

    for (int j=0; j<sizeof(Epsilon)/sizeof(Epsilon[0]); j++)
    {
        TQ aux;
        aux.S = 4;
        aux.pVar = Epsilon[j];
        aux.q = 0; 
        aux.k = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }

    for (int j=0; j<sizeof(LF)/sizeof(LF[0]); j++)
    {
        TQ aux;
        aux.S = 5;
        aux.pVar = LF[j];
        aux.q = 0;
        aux.k = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }

    for (int j=0; j<sizeof(DF)/sizeof(DF[0]); j++)
    {
        TQ aux;
        aux.S = 6;
        aux.pVar = DF[j];
        aux.q = 0; 
        aux.k = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }
}

void ChooseAction(int numGeneration)
{
    // choose actions for each state from Q-Table (e-Greedy)
    for (int i=0; i<par; i++)
    {
        int a, aMax, aAux;

        // set variable a with the current action
        if (i == 0)      a = a0;
        else if (i == 1) a = a1;
        else if (i == 2) a = a2;
        else if (i == 3) a = a3;       
        else if (i == 4) a = a4;    
        else if (i == 5) a = a5;    
        else if (i == 6) a = a6;       
                

        // found actions with the highest value of Q(i,-).q 
        double bQ = -INFINITY;  
        for (int j=0; j<Q[i].size(); j++)
        {
            if (Q[i][j].q > bQ)
            {
                bQ = Q[i][j].q;
                aAux = j;
            }
            else
            if (Q[i][j].q == bQ && randomico(0,1) >= 0.5) // trie randomly
            {
                aAux = j;
            }
            aMax = aAux;
        }

        // e-greedy policy
        if (randomico(0,1) >= epsilon && numGeneration > 30) // 
        {
            // choose the action with highest Q value
            a = aAux;
        }
        else
        {
            // choose a randonly selected other action (value)
            aAux = a;
            do
            {
                a = irandomico(0,Q[i].size()-1);
            } while (a == aAux);
        }

        // set new action
        switch (i)
        {
            case 0:
                a0 = a;
                a0Max = aMax;
                break;
        
            case 1:
                a1 = a;
                a1Max = aMax;
                break;
            
            case 2:
                a2 = a;
                a2Max = aMax;
                break; 

            case 3:
                a3 = a;
                a3Max = aMax;
                break;

            case 4:
                a4 = a;
                a4Max = aMax;
                break;

            case 5:
                a5 = a;
                a5Max = aMax;
                break;

            case 6:
                a6 = a;
                a6Max = aMax;
                break;
        }
    }
        
    // update parameters with actions 
    p       = Q[0][a0].pVar;
    pe      = Q[1][a1].pVar; 
    pm      = Q[2][a2].pVar; 
    rhoe    = Q[3][a3].pVar;  
    epsilon = Q[4][a4].pVar;  
    lf      = Q[5][a5].pVar;  
    df      = Q[6][a6].pVar;  
}

void UpdateQTable()
{
    for (int s=0; s<par; s++)
    {
        int a, aMax, aP1;

        double qTermination = 0;

        switch (s)
        {
            case 0:
                a = a0;
                aP1 = a1;
                aMax = a1Max;
                break;
            
            case 1:
                a = a1;
                aP1 = a2;
                aMax = a2Max;
                break;

            case 2:
                a = a2;
                aP1 = a3;
                aMax = a3Max;
                break;

            case 3:
                a = a3;
                aP1 = a4;
                aMax = a4Max;
                break;
            
            case 4:
                a = a4;
                aP1 = a5;
                aMax = a5Max;
                break;
            
            case 5:
                a = a5;
                aP1 = a6;
                aMax = a6Max;
                break;

            case 6:
                a = a6;
                aP1 = a0;
                aMax = a0Max;
                break;
        }
        
        qTotal -= Q[s][a].q;

        // update Q-Table only if R = 1
        if (R > 0)
        {
            //Q-Learning
            if (s == par-1)
            {
                Q[s][a].q += lf*(R + df*qTermination - Q[s][a].q ); 
            }
            else
            {
                Q[s][a].q += lf*(R + df*Q[s+1][aMax].q - Q[s][a].q );
            }
        }

        Q[s][a].k++;
        qTotal += Q[s][a].q;
    }
    //printf("%.3lf \n", qTotal);
}

TSol CreateInitialSolutions()
{
	TSol s;
	TVecRk aux;

        s.vec.clear();

	// create a random-key for each allelo (consider decoder type in the n-th random-key)
	for (int j = 0; j < n+1; j++)
	{
		aux.rk  = randomico(0,1);
        if (j < 2*nbUsers) { // Pickup or dropoff gene
            aux.user = j;
        }
        else if (j < 2*nbUsers+nbVehicles) { // vehicle genes
            aux.user = j-2*nbUsers;
        }
        else {
            aux.user = -1;
        }
        s.vec.push_back(aux);
	}
    
    // Initialize empty solution
    std::vector<TVecSol> tempRoute(0);
    std::vector<std::vector<TVecSol>> tempSol(nbVehicles, tempRoute);
    s.sol = tempSol;

    // flag to control the local search memory
    s.flag = 0;

    // Set infeasibility indicators to 0
    s.battery_infeasibles = 0;

   return s;
}

TSol Perturbation(TSol s, double beta) // Not used
{
    for (int k=0; k<n*beta; k++)
    {        
        // choose two random positions
        int pos1, pos2;
        do{
            pos1 = irandomico(0, n-1);
            pos2 = irandomico(0, n-1);
        } while(pos1 == pos2);

        // swap genes pos1 e pos2
        double temp = s.vec[pos1].rk;
        s.vec[pos1].rk = s.vec[pos2].rk;
        s.vec[pos2].rk = temp;
    }
    return s;
}

TSol ParametricUniformCrossover(int Tpe)
{	
    TSol s;

    // create a new offspring
    s.vec.resize(n+1);

    // Select an elite parent:
    int eliteParent = irandomico(0,Tpe - 1);

    // Select a non-elite parent:
    int noneliteParent = Tpe + irandomico(0, p - Tpe - 1);

    // Mate:  // including decoder in the n-th rk 
    for(int j = 0; j < n+1; j++)
    {
        //copy alelos of top chromossom of the new generation
        if (randomico(0,1) < rhoe)
           s.vec[j].rk = Pop[eliteParent].vec[j].rk;
        else
           s.vec[j].rk = Pop[noneliteParent].vec[j].rk;

        // "Clean" users
        if (j < 2*nbUsers)
        	s.vec[j].user = j;
        else if (j < 2*nbUsers + nbVehicles)
            s.vec[j].user = j - 2*nbUsers;
        else
        	s.vec[j].user = -1; // User ID = -1 for the rk's relating to scheduler and decoder

        // Initialize empty solution
        std::vector<TVecSol> tempRoute(0);
        std::vector<std::vector<TVecSol>> tempSol(nbVehicles, tempRoute);
        s.sol = tempSol;
    }

    // set the flag of local search as zero
    s.flag = 0;

    return s;
}

double PearsonCorrelation(std::vector <TVecRk> X, std::vector <TVecRk> Y)
{
    double correlation = 0;
    double sumXY = 0;
    double sumX2 = 0;
    double sumY2 = 0;
    double sumX = 0;
    double sumY = 0;

    for(int j=0; j<n-2; j++)
    {
        sumX += X[j].rk;
        sumX2 += X[j].rk * X[j].rk;
        sumXY += X[j].rk * Y[j].rk;
        sumY += Y[j].rk;
        sumY2 += Y[j].rk * Y[j].rk;
    }

    //Pearson
    correlation= ((n*sumXY) - (sumX*sumY) ) / (sqrt( (n*sumX2 - sumX*sumX) * (n*sumY2 - sumY*sumY) ));
    return correlation;
}

void writeLPGraph(std::vector<std::vector<std::pair<int, double> > > &listaArestas) 
{	
	// Colocar id dos grupos entre 1 e n
	int i, j;
	int n = listaArestas.size();
	int newId = 1;
	std::vector<int> newIds(n, 0);
	for (i = 0; i < n; i++) {
		if (newIds[Pop[i].label] == 0) {
			newIds[Pop[i].label] = newId;
			newId++;
		}
	}
	for (i = 0; i < n; i++) {
		Pop[i].label = newIds[Pop[i].label];
	}

	// Dá um nome ao arquivo .json (LP-<generation>.json)
        numLP++;
	std::string jsonfile = "LP-";
	jsonfile += std::to_string(numLP);
	jsonfile += ".json";

	// Cria o arquivo do grafo (.json)
	int totalArestas = 0;
	std::string json = "{\n  \"nodes\": [\n";
	for (i = 0; i < listaArestas.size(); i++) {			
		totalArestas += listaArestas[i].size();
		if(i == listaArestas.size() - 1)
			json += "    {\"id\": \"" + std::to_string(i) + "\", \"group\": " + std::to_string(Pop[i].label) + "}\n";
		else
			json += "    {\"id\": \"" + std::to_string(i) + "\", \"group\": " + std::to_string(Pop[i].label) + "},\n";
	}
	totalArestas /= 2;
	json += "  ],\n";
	json += "  \"links\": [\n";							// Insere arestas
	int totalArestasInseridas = 0;
	for (i = 0; i < listaArestas.size(); i++) {
		for (j = 0; j < listaArestas[i].size(); j++) {
			if (i < listaArestas[i][j].first) {
				
				// Create an output string stream
				std::ostringstream peso;
				// Set Fixed -Point Notation
				peso << std::fixed;
				// Set precision to 2 digits
				peso << std::setprecision(2);
				//Add double to stream
				peso << listaArestas[i][j].second;
				// Get string from output string stream
				std::string pesoString = peso.str();
				
				json += "    {\"source\": \"" + std::to_string(i) + "\", \"target\": \"" + std::to_string(listaArestas[i][j].first) + "\", \"value\": " + pesoString + "}";
				
				totalArestasInseridas++;
				if (totalArestasInseridas != totalArestas)
					json += ",\n";
				else
					json += "\n";
			}
		}
	}
	json += "  ]\n";
	json += "}";

	// Salva arquivo .json
	std::ofstream myfile;
	myfile.open("LP/"+jsonfile);
	myfile << json;
	myfile.close();

	// Dá um nome ao arquivo .json (LP-<generation>.json)
	std::string htmlfile = "LP-";
	htmlfile += std::to_string(numLP);
	htmlfile += ".html";

	// Cria o html que inclui a biblioteca de grafos e o grafo (.json)
	std::string html = "<head>\n";
	html += "\t<style> body{ margin : 0; } </style>\n";
	html += "\t<script src = \"//unpkg.com/three\"></script>\n";
	html += "\t<script src=\"//unpkg.com/three-spritetext\"></script>\n";
	html += "\t<script src = \"//unpkg.com/3d-force-graph\"></script>\n";
	html += "</head>\n";
	html += "<body>\n";
	html += "\t<div id=\"3d-graph\"></div>\n";
	html += "\t<script>\n";
	html += "\t\tconst Graph = ForceGraph3D()\n";
	html += "\t\t(document.getElementById('3d-graph'))\n";
	html += "\t\t.jsonUrl('" + jsonfile + "')\n";
	html += "\t\t.nodeLabel('id')\n";
	html += "\t\t.nodeAutoColorBy('group')\n";
	html += "\t\t.linkThreeObjectExtend(true)\n";
	html += "\t\t.linkThreeObject(link => {\n";
	html += "\t\t\tconst sprite = new SpriteText(`${link.value}`);\n";
	html += "\t\t\tsprite.color = 'lightgrey';\n";
	html += "\t\t\tsprite.textHeight = 1.5;\n";
	html += "\t\t\treturn sprite;\n";
	html += "\t\t})\n";
	html += "\t\t.linkPositionUpdate((sprite, { start, end }) => {\n";
	html += "\t\t\tconst middlePos = Object.assign(...['x', 'y', 'z'].map(c => ({\n";
	html += "\t\t\t[c]: start[c] + (end[c] - start[c]) / 2 })));\n";
	html += "\t\t\tObject.assign(sprite.position, middlePos);\n";
	html += "\t\t});\n";
	html += "\t\tGraph.d3Force('charge').strength(-120);\n";
	html += "\t</script>\n";
	html += "</body>";
	
	// Salva arquivo .html
	myfile.open("LP/" + htmlfile);
	myfile << html;
	myfile.close();
}

void IC(int Tpe) 
{
    std::vector<std::vector<std::pair<int, double> > > listaArestas(Tpe, std::vector<std::pair<int, double> >());

	// create weighted (pearson correlation) graph
	int entrouAresta = 0;
	double pearson = 0.0;
	for (int i = 0; i < Tpe - 1; i++) {
		for (int j = i + 1; j < Tpe; j++)
		{
			pearson = PearsonCorrelation(Pop[i].vec, Pop[j].vec);
			if (pearson > 0.7) {
				entrouAresta++;
				listaArestas[i].push_back(std::make_pair(j, pearson));
				listaArestas[j].push_back(std::make_pair(i, pearson));
			}
		}
	}

	// apply clustering method
	LP(listaArestas);

	PromisingLP(Tpe);
        listaArestas.clear();
}

void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas)
{
    int nk = listaArestas.size();

	// Cria o vetor para ordem de visita
	std::vector<int> ordemVisita(nk);
	iota(ordemVisita.begin(), ordemVisita.end(), 0);

	// initialize each node with its own label
	for (int i = 0; i < nk; i++)
		Pop[i].label = i;

	int iteracao = 1;
	int labelVizinho, melhorLabel;
	double melhorPeso;
	std::map<int, double> totalLabels;
	std::map<int, double>::iterator it;

	int movimentos = 1;
	while (movimentos) {
		movimentos = 0;
		random_shuffle(ordemVisita.begin(), ordemVisita.end());
		for (auto idVertice : ordemVisita) {

			// Calcula o peso para os labels
			totalLabels.clear();
			for (auto idVizinho : listaArestas[idVertice]) {
				labelVizinho = Pop[idVizinho.first].label;
				it = totalLabels.find(labelVizinho);
				if (it != totalLabels.end()) {
					it->second += idVizinho.second;
				}
				else {
					totalLabels[labelVizinho] = idVizinho.second;
				}
			}

			// Melhor label é ele mesmo inicialmente
			melhorLabel = Pop[idVertice].label;
			melhorPeso = std::numeric_limits<double>::min();
			for (auto totais : totalLabels) {
				if (totais.second > melhorPeso) {
					melhorLabel = totais.first;
					melhorPeso = totais.second;
				}
			}

			if (melhorLabel != Pop[idVertice].label) {
				Pop[idVertice].label = melhorLabel;
				movimentos = 1;
			}
		}
		iteracao++;
	}
    ordemVisita.clear();
}

void PromisingLP(int Tpe)
{
    	std::vector <int> grupos;
	int tamanhoGrupos = 0;

	// initialize promisings solutions
	for (int i = 0; i < Tpe; i++)
		Pop[i].promising = 0;

	// save labels defined by LP in groups
	int achei;

    	for (int i = 0; i < Tpe; i++)
	{
		achei = 0;
		for (unsigned int j = 0; j < grupos.size(); j++)
		{
			if (Pop[i].label == grupos[j])
				achei = 1;
		}
		if (achei == 0)
		{
			tamanhoGrupos++;
			grupos.push_back(Pop[i].label);
		}
	}

	// find the best solution in the group (with flag = 0)
	for (unsigned int j = 0; j < grupos.size(); j++)
	{
		double menorFO = INFINITY;
		int localMenor = -1;
		int local = -1;
		for (int i = 0; i < Tpe; i++)
		{
			if (Pop[i].label == grupos[j])
			{
				// find the best solution of the group
				if (local == -1)
					local = i;

				// we not apply local search in this solution yet
                		if (Pop[i].fo < menorFO && Pop[i].flag == 0) 
				{
					menorFO = Pop[i].fo;
					localMenor = i;
				}
			}
		}

		if (localMenor == -1)
			localMenor = local;

		if (Pop[localMenor].label != -1)
			Pop[localMenor].promising = 1;
	}
}

void FreeMemory()
{
    //methods
    Pop.clear();
    PopInter.clear();
    Q.clear();
}

double randomico(double min, double max)
{
    return ((double)(rand()%10000)/10000.0)*(max-min)+min;
    //return uniform_real_distribution<double>(min, max)(rng);
}

int irandomico(int min, int max)
{
    return (int)randomico(0,max-min+1.0) + min;
}

void splitF0(TSol& s, std::vector <TNode>& node, std::vector <std::vector <double> >& dist, int nbVehicles, int nbUsers, int periodLength, std::vector <int> maxRideTimes, double weights[3], std::vector <double> elecPrices) {
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
                URTs[currStop - nbUsers] += s.sol[v][i].T; // Add start-time of service at pick-up
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
        ERT += URTs[i] - dist[i][i + nbUsers]; // Actual Ride Time - Ideal Ride Time
    }
    
    s.TRT = TRT;
    s.ERT = ERT;
    s.CC = CC;
}
