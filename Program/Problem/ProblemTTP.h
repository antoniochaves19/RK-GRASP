// *******************************************************************
//      file with specific functions to solve the TTP
//                  (Travelling Thief Problem)
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Global Variables
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     k - number of local search heuristics
extern int n;                               // size of the vector solution

//---------------------- DEFINITION OF TYPES OF PROBLEM SPECIFIC --------------------

struct Item
{
    int Index;               // Item's index
    double Profit;           // Item's profit
    double Weight;           // Item's weight
    int Selected = 0;        // False if the item is not picked in the packing plan, and true otherwise
};

struct City
{
    int Index;                           // Index of the city
    std::vector<Item> Items;             // Set of items assigned to the city
    int Distance;                        // Distance to the next city in the tour. REMARK: distance for Tour[n-1] corresponds to the distance between cities (n-1) and 0 of solution Tour
    int m;                               // Number of items in the city
    double PositionX;                    // Position X of the city
    double PositionY;                    // Position Y of the city
};

//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  ----------

static std::vector <City> Cities;            // Set of cities in the instance
static int nCities;                          // Number of cities in the instance
static int mItems;                           // Total number of items in the instance
static double Rate;                          // Renting rate
static double MinSpeed;                      // Minimal speed
static double MaxSpeed;                      // Maximal speed
static double MaxWeight;                     // Maximal knapsack's weight
static double MaxProfit;

static std::vector <std::vector <double> > dist;	    // matrix with Euclidean distance
static std::vector <std::vector <int> > KPvector;	    // matrix with possible solutions of the Knapsak problem
static int KP;                                          // number of KP possibilities


//----------------------- IMPLEMENTATION OF FUNCTIONS  -------------------------------


/************************************************************************************
 Method: ReadData
 Description: read input data of the problem
*************************************************************************************/
void ReadData(char nameTable[])
{ 
    char name[200] = "../Instances/";
    strcat(name,nameTable);

    FILE *file;
    file = fopen(name,"r");

    if (file == NULL)
    {
        printf("\nERROR: File (%s) not found!\n",name);
        getchar();
        exit(1);
    }

    // => read data

    // read instance head
    fscanf(file, "%d", &nCities);
    fscanf(file, "%d", &mItems);
    fscanf(file, "%lf", &MaxWeight);
    fscanf(file, "%lf", &MinSpeed);
    fscanf(file, "%lf", &MaxSpeed);
    fscanf(file, "%lf", &Rate);

    // printf("\n%d \n%d \n%lf \n%lf \n%lf \n%lf \n", nCities, mItems, MaxWeight, MinSpeed, MaxSpeed, Rate);

    // read node informations
    // // NODE_COORD_SECTION	(INDEX, X, Y): 
    // // Assuming City and Item structs are defined
    Cities.clear();
    Cities.resize(nCities);
    for (int i = 0; i < nCities; i++) {
        fscanf(file, "%d %lf %lf", &Cities[i].Index, &Cities[i].PositionX, &Cities[i].PositionY);
        // printf("\n%d: %d \t %lf \t %lf", i, Cities[i].Index, Cities[i].PositionX, Cities[i].PositionY);
    }

    // ITEMS SECTION	(INDEX, PROFIT, WEIGHT, ASSIGNED NODE NUMBER): 
    MaxProfit = 0;
    for (int i = 0; i < mItems; i++) {
        Item item;
        int NodeIndex = 0;

        fscanf(file, "%d %lf %lf %d", &item.Index, &item.Profit, &item.Weight, &NodeIndex);
        Cities[NodeIndex - 1].Items.push_back(item);

        MaxProfit += item.Profit;

        // printf("\n%d %lf %lf %d", item.Index, item.Profit, item.Weight, NodeIndex);
    }

    fclose(file);


    // print input data
    // if (debug){
    //     printf("\n\n");
    //     printf("\n%d \n%d \n%.0lf \n%.2lf \n%.2lf \n%.2lf \n", nCities, mItems, MaxWeight, MinSpeed, MaxSpeed, Rate);
    //     for (int i = 0; i < nCities; i++) {
    //         printf("\n%d \t %.0lf \t %.0lf", Cities[i].Index, Cities[i].PositionX, Cities[i].PositionY);
    //         for (unsigned int j=0; j<Cities[i].Items.size(); j++)
    //         {
    //             printf("\t[ %d \t %.0lf \t %.0lf ] ", Cities[i].Items[j].Index, Cities[i].Items[j].Profit, Cities[i].Items[j].Weight);
    //         }
    //     }
    // }

    // ***** gerar os padroes possiveis considerando o numero de itens por cidade ****
    KP = std::ceil((double)mItems/nCities);
    // if (debug) printf("\nK = %d\n",KP);

    // Número total de combinações possíveis para KP bits é 2^k
    int totalCombinations = std::pow(2, KP);

    KPvector.clear();
    KPvector.resize(totalCombinations, std::vector<int>(KP));

    // Iterar sobre todas as combinações
    for (int i = 0; i < totalCombinations; ++i) {
        // Converter o número atual para um vetor binário
        for (int j = 0; j < KP; ++j) {
            if (i & (1 << j)) {
                KPvector[i][KP - j - 1] = 1;  // Definir bit na posição correta
            }
        }

        // Imprimir o vetor binário
        // printf("\n");
        // for (int j=0; j<KPvector[i].size(); j++) {
        //     printf("%d ", KPvector[i][j]);
        // }
    }

    // calculate the euclidean distance
    dist.clear();
    dist.resize(nCities, std::vector<double>(nCities));

    for (int i=0; i<nCities; i++)
    {
    	for (int j=i; j<nCities; j++)
    	{
    		dist[i][j] = dist[j][i] = (floor (sqrt( (Cities[j].PositionX - Cities[i].PositionX) * (Cities[j].PositionX - Cities[i].PositionX) +
    										        (Cities[j].PositionY - Cities[i].PositionY) * (Cities[j].PositionY - Cities[i].PositionY) ) + 0.5 ) )/1.0;

            // dist[i][j] = dist[j][i] = ceil (sqrt( (Cities[j].PositionX - Cities[i].PositionX) * (Cities[j].PositionX - Cities[i].PositionX) +
    										    //   (Cities[j].PositionY - Cities[i].PositionY) * (Cities[j].PositionY - Cities[i].PositionY) ) );
    	}
    }

    // printf("\nDistancias:\n");
    // for (int i=0; i<nCities; i++)
    // {
    // 	for (int j=0; j<nCities; j++)
    // 	{
    // 		printf("%.2lf\t", dist[i][j]);
    // 	}
    //     printf("\n");
    // }

    
    n = nCities+nCities;
    // getchar();
}

/*void ReadData(char nameTable[])
{ 
    // => read data
    FILE *file;
    char FileName[] = "../Instances/";
    strcat(FileName,nameTable);
    // char input[256];
    nCities = 0;
    mItems = 0;
    Rate = 0; 
    MinSpeed = 0; 
    MaxSpeed = 0;
    MaxWeight = 0;

    file = fopen(FileName, "r");
    if (file == NULL) {
        perror("Error opening file");
        getchar();
        exit(1);
    }

    // PROBLEM NAME: 	a280-TTP
    // KNAPSACK DATA TYPE: bounded strongly corr
    // DIMENSION:	280
    // NUMBER OF ITEMS: 	279
    // CAPACITY OF KNAPSACK: 	25936
    // MIN SPEED: 	0.1
    // MAX SPEED: 	1
    // RENTING RATIO: 	5.61
    // EDGE_WEIGHT_TYPE:	CEIL_2D
    // NODE_COORD_SECTION	(INDEX, X, Y): 

    // while (fgets(input, 256, file) != NULL && !strstr(input, "NODE_COORD_SECTION")) {
    //     if (strstr(input, "DIMENSION:")) sscanf(input, "%*[^:]: %d", &nCities);
    //     if (strstr(input, "NUMBER OF ITEMS:")) sscanf(input, "%*[^:]: %d", &mItems);
    //     if (strstr(input, "CAPACITY OF KNAPSACK:")) sscanf(input, "%*[^:]: %lf", &MaxWeight);
    //     if (strstr(input, "MIN SPEED:")) sscanf(input, "%*[^:]: %lf", &MinSpeed);
    //     if (strstr(input, "MAX SPEED:")) sscanf(input, "%*[^:]: %lf", &MaxSpeed);
    //     if (strstr(input, "RENTING RATIO:")) sscanf(input, "%*[^:]: %lf", &Rate);
    // }

    fscanf(file, "%d", &nCities);
    fscanf(file, "%d", &mItems);
    fscanf(file, "%lf", &MaxWeight);
    fscanf(file, "%lf", &MinSpeed);
    fscanf(file, "%lf", &MaxSpeed);
    fscanf(file, "%lf", &Rate);

    // printf("\n%d \n%d \n%lf \n%lf \n%lf \n%lf \n", nCities, mItems, MaxWeight, MinSpeed, MaxSpeed, Rate);

    
    // fgets(input, sizeof(input), file); // skip NODE_COORD_SECTION part
    // // NODE_COORD_SECTION	(INDEX, X, Y): 
    // // Assuming City and Item structs are defined
    // Cities.clear();
    // Cities.resize(nCities);
    // for (int i = 0; i < nCities; i++) {
    //     fscanf(file, "%d %lf %lf", &Cities[i].Index, &Cities[i].PositionX, &Cities[i].PositionY);
        // printf("\n%d: %d \t %lf \t %lf", i, Cities[i].Index, Cities[i].PositionX, Cities[i].PositionY);
    // }

    // fgets(input, 256, file); // skip ITEM_SECTION part
    // fgets(input, 256, file); // skip ITEM_SECTION part
    
    // // ITEMS SECTION	(INDEX, PROFIT, WEIGHT, ASSIGNED NODE NUMBER): 

    // for (int i = 0; i < mItems; i++) {
    //     Item item;
    //     int NodeIndex = 0;

    //     fscanf(file, "%d %lf %lf %d", &item.Index, &item.Profit, &item.Weight, &NodeIndex);
    //     Cities[NodeIndex - 1].Items.push_back(item);

    //     printf("\n%d %lf %lf %d", item.Index, item.Profit, item.Weight, NodeIndex);
    // }

    fclose(file);

    // if (debug){
    //     printf("\n\n");
    //     printf("\n%d \n%d \n%.0lf \n%.2lf \n%.2lf \n%.2lf \n", nCities, mItems, MaxWeight, MinSpeed, MaxSpeed, Rate);
    //     for (int i = 0; i < nCities; i++) {
    //         printf("\n%d \t %.0lf \t %.0lf", Cities[i].Index, Cities[i].PositionX, Cities[i].PositionY);
    //         for (unsigned int j=0; j<Cities[i].Items.size(); j++)
    //         {
    //             printf("\t[ %d \t %.0lf \t %.0lf ] ", Cities[i].Items[j].Index, Cities[i].Items[j].Profit, Cities[i].Items[j].Weight);
    //         }
    //     }
    //     getchar();
    // }


    // // ***** gerar os padroes possiveis considerando o numero de itens por cidade ****
    // KP = std::ceil((double)mItems/nCities);
    // // if (debug) printf("\nK = %d\n",KP);

    // // Número total de combinações possíveis para KP bits é 2^k
    // int totalCombinations = std::pow(2, KP);

    // KPvector.clear();
    // KPvector.resize(totalCombinations, std::vector<int>(KP));

    // // Iterar sobre todas as combinações
    // for (int i = 0; i < totalCombinations; ++i) {
    //     // Converter o número atual para um vetor binário
    //     for (int j = 0; j < KP; ++j) {
    //         if (i & (1 << j)) {
    //             KPvector[i][KP - j - 1] = 1;  // Definir bit na posição correta
    //         }
    //     }

    //     // Imprimir o vetor binário
    //     printf("\n");
    //     for (int j=0; j<KPvector[i].size(); j++) {
    //         printf("%d ", KPvector[i][j]);
    //     }
    // }
    // // getchar();

    // // calculate the euclidean distance
    // dist.clear();
    // dist.resize(nCities, std::vector<double>(nCities));

    // for (int i=0; i<nCities; i++)
    // {
    // 	for (int j=i; j<nCities; j++)
    // 	{
    // 		dist[i][j] = dist[j][i] = (floor (sqrt( (Cities[j].PositionX - Cities[i].PositionX) * (Cities[j].PositionX - Cities[i].PositionX) +
    // 										        (Cities[j].PositionY - Cities[i].PositionY) * (Cities[j].PositionY - Cities[i].PositionY) ) + 0.5 ) )/1.0;
    // 	}
    // }

    // printf("\nDistancias:\n");
    // for (int i=0; i<nCities; i++)
    // {
    // 	for (int j=0; j<nCities; j++)
    // 	{
    // 		printf("%.2lf\t", dist[i][j]);
    // 	}
    //     printf("\n");
    // }
    // getchar();
    

    // define the random-key vector size
    n = nCities+nCities;

    printf("\n n = %d", n);
    // getchar();
}*/

/************************************************************************************
 Method: CalculateFitness
 Description: calculate the fitness of a chromossom s
*************************************************************************************/
double CalculateFitness(TSol s, int print)
{
    // calculate objective function
    s.ofv = 0;

    // output
    std::vector<int> cityVector;
    std::vector<int> itemVector;

    double CollectedWeight = 0;                                 // Collected weight while traveling along Tour
    double SpeedCoef = (MaxSpeed - MinSpeed) / MaxWeight;       // constant

    double fp = 0;
    double ft = 0;

    // foreach (City in Tour)
    for (int i=0; i<nCities; i++)
    {
        int city = s.vec[i].sol;
        int packing = s.vec[nCities+i].sol;

        // output
        cityVector.push_back(city+1);

        // foreach (Item in city)
        for (unsigned int j=0; j<KPvector[packing].size(); j++)
        {
            // if (item.Selected)
            if (KPvector[packing][j] == 1 && city > 0)
            {
                if (Cities[city].Items[j].Weight + CollectedWeight <= MaxWeight)
                {
                    CollectedWeight += Cities[city].Items[j].Weight;
                    fp += Cities[city].Items[j].Profit;

                    // output
                    itemVector.push_back(Cities[city].Items[j].Index);
                }
            }
        }
        
        ft += (dist[s.vec[i%nCities].sol][s.vec[(i+1)%nCities].sol]) / (MaxSpeed - CollectedWeight * SpeedCoef);
    }

    s.ofv =  (fp - Rate * ft) * -1;

    // output
    if (print)
    {
        printf("\nOFV: %lf \n[",s.ofv*-1);
        for (unsigned int i = 0; i < cityVector.size(); i++){
            if (i < cityVector.size()-1)
                printf("%d, ", cityVector[i]);
            else
                printf("%d", cityVector[i]);
        }
        printf("]");

        printf("\n[");
        for (unsigned int i = 0; i < itemVector.size(); i++){
            if (i < itemVector.size()-1)
                printf("%d, ", itemVector[i]);
            else
                printf("%d", itemVector[i]);
        }
        printf("]\n");
    }
    

    return s.ofv;
}

/************************************************************************************
 Method: Encode
 Description: encode a solution s into a random key vector
*************************************************************************************/
void Encode(TSol &s)
{
    // // encode solution sol into random key
    // double delta = 1.0 / n;         // Size of each chunk
    // double X_bar = delta / 2.0;     // Center of the first chunk

    // // Loop through the input sequence π
    // for (int i = 0; i < n; i++) {
    //     // Assign center values to appropriate positions in s
    //     s.vec[s.vec[i].sol].rk = X_bar;

    //     // Move to the center of the next chunk
    //     X_bar += delta;
    // }

    // // Randomize the protocol by adding uniform noise to each element in s
    // for (int i = 0; i < n; i++) {
    //     double delta_i = (double)rand() / RAND_MAX * delta - delta / 2.0;
    //     s.vec[i].rk += delta_i;
    // }

    // s.vec[n].rk = 0.0;
}

/************************************************************************************
 Method: Decoders
 Description: users need to implement at least one decoder, DecK (K = [1,2,3,4,5])
*************************************************************************************/

void Dec1A(TSol &s) // sort
{
    // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < nCities)
		    s.vec[j].sol = j;
        else
        if (j < nCities+nCities)
            s.vec[j].sol = 0;
        else
            s.vec[j].sol = -1;
	}

    // sort random-key vector with nCities
    sort(s.vec.begin()+1, s.vec.begin() + nCities-1, sortByRk);   

    // assign the KP combination to each city
    for (int i=0; i<nCities; i++)
    {
        s.vec[nCities+i].sol = floor(s.vec[nCities+i].rk * (int)KPvector.size());
    }

    // printf("\n");
    // for (int i = 0; i < n+1; i++)
    // {
    //     if (i == nCities)
    //         printf(" | ");
    //     printf("%d ", s.vec[i].sol);
    // }
    // getchar();
    

    s.ofv = CalculateFitness(s,0);
}

void Dec1V(TSol &s) // k-Nearest Insertion
{
    double CollectedWeight = 0;                                 // Collected weight while traveling along Tour
    double SpeedCoef = (MaxSpeed - MinSpeed) / MaxWeight;       // constant
    std::vector<double> WeightCity;

     // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < nCities)
		    s.vec[j].sol = j;
        else
        if (j < nCities+nCities)
            s.vec[j].sol = 0;
        else
            s.vec[j].sol = -1;
	}

    // sort random-key vector with nCities
    sort(s.vec.begin()+1, s.vec.begin() + nCities-1, sortByRk);   


    // assign the KP combination to each city
    WeightCity.resize(nCities,0);
    for (int i=0; i<nCities; i++){
        s.vec[nCities+i].sol = floor(s.vec[nCities+i].rk * (int)KPvector.size());

        // foreach (Item in city)
        int packing = s.vec[nCities+i].sol;
        for (unsigned int j=0; j<KPvector[packing].size(); j++)
        {
            // if (item.Selected)
            if (KPvector[packing][j] == 1 && i > 0)
            {
                WeightCity[i] += Cities[i].Items[j].Weight;
            }
        }
    }

    // list of candidates
    std::vector <int> sC;
    sC.resize(nCities-1);

    for (int i = 1; i < nCities; i++){
        sC[i]= s.vec[i].sol;
    }


    // (dist[s.vec[i%nCities].sol][s.vec[(i+1)%nCities].sol]) / (MaxSpeed - CollectedWeight * SpeedCoef);

    // construct a solution with k Nearest insertion
    int pos = 0;
    while (!sC.empty())
    {
        // find the point i nearest from the route into the k first points of sC
        int i = 0;
        double costNearest = INFINITO;

        for (unsigned int k=0; k<3 && k<sC.size(); k++)
        {
            if (dist[s.vec[pos].sol][sC[k]] / (MaxSpeed - (CollectedWeight + WeightCity[sC[k]]) * SpeedCoef) < costNearest) //&&(CollectedWeight + WeightCity[sC[k]]) < MaxWeight)
            {
                costNearest = dist[s.vec[pos].sol][sC[k]] / (MaxSpeed - (CollectedWeight + WeightCity[sC[k]]) * SpeedCoef);
                i = k;
            }
        }

        // insert the i-th point in the position pos+1
        s.vec[pos+1].sol = sC[i]; 
        
        // update the collected weight
        CollectedWeight += WeightCity[sC[i]];

        // erase the i-th point of the sC list
        sC.erase(sC.begin()+i);
    }

    s.ofv = CalculateFitness(s,0);
}


void Dec1(TSol &s) // Cheapest Insertion
{
     // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < nCities)
		    s.vec[j].sol = j;
        else
        if (j < nCities+nCities)
            s.vec[j].sol = 0;
        else
            s.vec[j].sol = -1;
	}

    // sort random-key vector with nCities
    sort(s.vec.begin()+1, s.vec.begin() + nCities-1, sortByRk);   

    // assign the KP combination to each city
    for (int i=0; i<nCities; i++){
        s.vec[nCities+i].sol = floor(s.vec[nCities+i].rk * (int)KPvector.size());
    }

    // list of candidates
    std::vector <int> sC;
    sC.resize(nCities-1);

    for (int i = 3; i < nCities; i++){
        sC[i]= s.vec[i].sol;
    }

    std::vector <int> sol;
    sol.resize(nCities);
    sol[0] = s.vec[0].sol;
    sol[1] = s.vec[1].sol;
    sol[2] = s.vec[2].sol;

    // construct a solution with cheapest insertion
    int sizeAtual = 3;
    for (int i = 3; i<nCities; i++)
    {
        // find the cheapest position to insert the i-th point of sC
        int bestPosition = 0;
        float costBest = INFINITO;
        float costInsertion = 0;
        for (unsigned int j = 0; j<sizeAtual; j++)
        {
            if (j == sizeAtual-1)
            {
                // cost to insert between i-1 and 0
                costInsertion = dist[sol[j]][sC[i]] + dist[sC[i]][sol[0]] - dist[sol[j]][sol[0]];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
            else
            {
                // cost to insert between i and i+1
                costInsertion = dist[sol[j]][sC[i]] + dist[sC[i]][sol[j+1]] - dist[sol[j]][sol[j+1]];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
        }

        // insert the i-th point in the cheapest position
        sol.insert(sol.begin()+bestPosition+1,sC[i]);
        sizeAtual++;
    }

    for (int i = 0; i < nCities; i++)
    {
        s.vec[i].sol = sol[i];
    }
    
    s.ofv = CalculateFitness(s,0);
}


void Dec2(TSol &s) {}

void Dec3(TSol &s) {}

void Dec4(TSol &s) {}

void Dec5(TSol &s) {}


/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem()
{
    //specific problem
    dist.clear();
    KPvector.clear();
    Cities.clear();
}

#endif