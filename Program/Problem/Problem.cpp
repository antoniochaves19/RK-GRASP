// *******************************************************************
//      arquivo com as funcoes especificas para resolver o TSP 
// *******************************************************************


#include "Problem.h"

// Sort TSol by random-keys
bool sortByRk(const TVecSol &lhs, const TVecSol &rhs) { return lhs.rk < rhs.rk; }

void ReadData(char nameTable[])
{ 
    char name[200] = "../Instances/";
    strcat(name,nameTable);

    FILE *arq;
    arq = fopen(name,"r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n",name);
        getchar();
        exit(1);
    }

    // => read data

    // read instance head
    char temp[100];
    fgets(temp, sizeof(temp), arq);
    fgets(temp, sizeof(temp), arq);
    fgets(temp, sizeof(temp), arq);
    fgets(temp, sizeof(temp), arq);
    fgets(temp, sizeof(temp), arq);
    fgets(temp, sizeof(temp), arq);

    // read node informations
    int nAux = 0;
    node.clear();
    TNode nodeTemp;

    while (!feof(arq))
    {
    	fscanf(arq, "%d %lf %lf", &nodeTemp.id, &nodeTemp.x, &nodeTemp.y);
    	node.push_back(nodeTemp);

    	nAux++;
    }
    fclose(arq);

    // calculate the euclidean distance
    dist.clear();
    dist.resize(nAux, std::vector<double>(nAux));

    for (int i=0; i<nAux; i++)
    {
    	for (int j=i; j<nAux; j++)
    	{
    		dist[i][j] = dist[j][i] = (floor (sqrt( (node[j].x - node[i].x) * (node[j].x - node[i].x) +
    										        (node[j].y - node[i].y) * (node[j].y - node[i].y) ) + 0.5 ) )/1.0;
    	}
    }
    
    n = nAux;
}

void Decoder(TSol &s)
{
    // copy the random-key sequence of current solution 
    TSol temp = s;

    // define decoder function based in the random-key of position n+1
    int dec = 1;//floor(s.vec[n].rk*numDecoders)+1;
    switch (dec)
    {
        case 1: 
            Dec1(s);
            break;

        case 2: 
            Dec2(s);
            break;
        
        case 3: 
            Dec3(s);
            break;

        case 4: 
            Dec4(s);
            break;

        case 5: 
            Dec5(s);
            break;

        default:
            break;
    }

    // return initial random-key sequence and maintain the solution sequence
    for (int i=0; i<n; i++){
        s.vec[i].rk = temp.vec[i].rk;
    }
}

void LocalSearch(TSol &s)
{
    // ***** we use a Random Variable Neighborhood Descent (RVND) as local search ****

    // current neighborhood
	int k = 1;

    // predefined number of neighborhood moves
    std::vector <int> NSL;
    std::vector <int> NSLAux;
    
    for (int i=1; i<=numLS; i++)
    {
        NSL.push_back(i);
        NSLAux.push_back(i);
    }

	while (!NSL.empty())
	{
        // current objective function
        double foCurrent = s.ofv;

        // randomly choose a neighborhood
        int pos = rand() % NSL.size();
        k = NSL[pos];

        switch (k)
        {
            case 1: 
                LS1(s); 
                break;

            case 2:
                LS2(s); 
                break;

            case 3:
                LS3(s); 
                break;

            case 4:
                LS4(s); 
                break;

            default:
                break;
        }

        // we better the current solution
        if (s.ofv < foCurrent)
        {
            // refresh NSL
            NSL.clear();
            NSL = NSLAux;
        }
        else
        {
            // Remove N(n) from NSL
            NSL.erase(NSL.begin()+pos);
        }
	} //end while

    // encode solution sol into random key
    double delta = 1.0 / n;         // Size of each chunk
    double X_bar = delta / 2.0;     // Center of the first chunk

    // Loop through the input sequence π
    for (int i = 0; i < n; i++) {
        // Assign center values to appropriate positions in s
        s.vec[s.vec[i].sol].rk = X_bar;

        // Move to the center of the next chunk
        X_bar += delta;
    }

    // Randomize the protocol by adding uniform noise to each element in s
    for (int i = 0; i < n; i++) {
        double delta_i = (double)rand() / RAND_MAX * delta - delta / 2.0;
        s.vec[i].rk += delta_i;
    }

    s.vec[n].rk = 0.0;
}

double CalculateFitness(TSol s)
{
    // calculate objective function
    s.ofv = 0;
    for (int i=0; i<n; i++){
        s.ofv += dist[s.vec[i%n].sol][s.vec[(i+1)%n].sol];
    }

    return s.ofv;
}

void Dec1(TSol &s) // sort
{
    // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < n)
		    s.vec[j].sol = j;
        else
            s.vec[j].sol = -1;
	}

    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);   

    s.ofv = CalculateFitness(s);
}

void Dec2(TSol &s) // 2-Opt
{
    // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < n)
		    s.vec[j].sol = j;
        else
            s.vec[j].sol = -1;
	}

    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);

    int t = 0, i = 0, j = 0, Mi1= 0, Mj = 0;

    float foOpt = 0;

    t = n; // use a circular list
    for (i=0; i < t; i++)
    {
        j = i + 2;
        while (((j+1)%t) != i)
        {
        int vi  = s.vec[i].sol;
        int vi1 = s.vec[(i+1)%t].sol;
        int vj  = s.vec[j%t].sol;
        int vj1 = s.vec[(j+1)%t].sol;

        foOpt = - dist[vi][vi1]
                - dist[vj][vj1]
                + dist[vi][vj]
                + dist[vi1][vj1];

        if (foOpt < 0)
        {
            // first improvement strategy
            Mi1 = (i+1)%t;
            Mj  = j%t;

            int inicio = Mi1,
                fim = Mj;

            int tam, p1, p2, aux;

            if(inicio > fim)
                tam = t - inicio + fim + 1;
            else
                tam = fim - inicio + 1;

            p1=inicio;
            p2=fim;

            for(int k=0; k < tam/2; k++)
            {
                aux = s.vec[p1%t].sol;
                s.vec[p1%t].sol = s.vec[p2%t].sol;
                s.vec[p2%t].sol = aux;

                p1 = (p1==t-1)?0:p1+1;
                p2 = (p2 == 0)?t-1:p2-1;
            }
        }
        j++;
        }//while
    }//for

    s.ofv = CalculateFitness(s);
}

void Dec3(TSol &s) // Cheapest Insertion
{
    // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < n)
		    s.vec[j].sol = j;
        else
            s.vec[j].sol = -1;
	}

    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);

    TVecSol aux = s.vec[n];

    // order list of candidates
    TSol sC = s;

    // partial route with three points
    s.vec.resize(3);

    // construct a solution with cheapest insertion
    for (int i = 3; i<n; i++)
    {
        // find the cheapest position to insert the i-th point of sC
        int bestPosition = 0;
        float costBest = INFINITO;
        float costInsertion = 0;
        for (unsigned int j = 0; j<s.vec.size(); j++)
        {
            if (j == s.vec.size()-1)
            {
                // cost to insert between i-1 and 0
                costInsertion = dist[s.vec[j].sol][sC.vec[i].sol] + dist[sC.vec[i].sol][s.vec[0].sol] - dist[s.vec[j].sol][s.vec[0].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
            else
            {
                // cost to insert between i and i+1
                costInsertion = dist[s.vec[j].sol][sC.vec[i].sol] + dist[sC.vec[i].sol][s.vec[j+1].sol] - dist[s.vec[j].sol][s.vec[j+1].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
        }

        // insert the i-th point in the cheapest position
        s.vec.insert(s.vec.begin()+bestPosition+1,sC.vec[i]);
    }

    // last RK
    s.vec.push_back(aux);

    s.ofv = CalculateFitness(s);
}

void Dec4(TSol &s) // k-Farthest Insertion
{
    // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < n)
		    s.vec[j].sol = j;
        else
            s.vec[j].sol = -1;
	}

    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);

    // order list of candidates
    std::vector <TVecSol> sC = s.vec;
    sC.erase(sC.begin()+n); // apagar ultima chave de sC
    TVecSol aux = s.vec[n]; // copiar ultima chave de rk

    //TSol temp = s;

    // partial route with one point
    s.vec.clear();
    s.vec.push_back(sC[0]);
    sC.erase(sC.begin());

    // construct a solution with k farthest insertion
    while (!sC.empty())
    {
        // find the point i farthest from the partial route into the k first points of sC
        int i = 0;
        double costFarthest = -INFINITO;

        for (unsigned int k=0; k<3 && k<sC.size(); k++)
        {
            for (unsigned int j = 0; j<s.vec.size(); j++)
            {
                if (dist[s.vec[j].sol][sC[k].sol] > costFarthest)
                {
                    costFarthest = dist[s.vec[j].sol][sC[k].sol];
                    i = k;
                }
            }
        }

        // find the cheapest position to insert the point i into the partial route
        int bestPosition = 0;
        float costBest = INFINITO;
        float costInsertion = 0;
        for (unsigned int j = 0; j<s.vec.size(); j++)
        {
            if (j == s.vec.size()-1)
            {
                // cost to insert between n-1 and 0
                costInsertion = dist[s.vec[j].sol][sC[i].sol] + dist[sC[i].sol][s.vec[0].sol] - dist[s.vec[j].sol][s.vec[0].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
            else
            {
                // cost to insert between j and j+1
                costInsertion = dist[s.vec[j].sol][sC[i].sol] + dist[sC[i].sol][s.vec[j+1].sol] - dist[s.vec[j].sol][s.vec[j+1].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
        }

        // insert the i-th point in the cheapest position
        s.vec.insert(s.vec.begin()+bestPosition+1,sC[i]); //

        // erase the i-th point of the sC list
        sC.erase(sC.begin()+i);
    }

    // last random-key
    s.vec.push_back(aux);

    s.ofv = CalculateFitness(s);
}

void Dec5(TSol &s) // k-Nearest Insertion
{
    // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < n)
		    s.vec[j].sol = j;
        else
            s.vec[j].sol = -1;
	}

    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);

    // order list of candidates
    std::vector <TVecSol> sC = s.vec;
    sC.erase(sC.begin()+n); // apagar ultima chave de sC
    TVecSol aux = s.vec[n]; // copiar ultima chave de rk

    //TSol temp = s;

    // partial route with one point
    s.vec.clear();
    s.vec.push_back(sC[0]);
    sC.erase(sC.begin());

    // construct a solution with k farthest insertion
    while (!sC.empty())
    {
        // find the point i nearest from the partial route into the k first points of sC
        int i = 0;
        double costNearest = INFINITO;

        for (unsigned int k=0; k<3 && k<sC.size(); k++)
        {
            for (unsigned int j = 0; j<s.vec.size(); j++)
            {
                if (dist[s.vec[j].sol][sC[k].sol] < costNearest)
                {
                    costNearest = dist[s.vec[j].sol][sC[k].sol];
                    i = k;
                }
            }
        }

        // find the cheapest position to insert the point i into the partial route
        int bestPosition = 0;
        float costBest = INFINITO;
        float costInsertion = 0;
        for (unsigned int j = 0; j<s.vec.size(); j++)
        {
            if (j == s.vec.size()-1)
            {
                // cost to insert between n-1 and 0
                costInsertion = dist[s.vec[j].sol][sC[i].sol] + dist[sC[i].sol][s.vec[0].sol] - dist[s.vec[j].sol][s.vec[0].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
            else
            {
                // cost to insert between j and j+1
                costInsertion = dist[s.vec[j].sol][sC[i].sol] + dist[sC[i].sol][s.vec[j+1].sol] - dist[s.vec[j].sol][s.vec[j+1].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
        }

        // insert the i-th point in the cheapest position
        s.vec.insert(s.vec.begin()+bestPosition+1,sC[i]); //

        // erase the i-th point of the sC list
        sC.erase(sC.begin()+i);
    }

    // last random-key
    s.vec.push_back(aux);

    s.ofv = CalculateFitness(s);
}

void LS1(TSol &s) // 2-Opt
{
    int t = n, // use a circular list
        i = 0,
        j = 0,
        Mi1= 0,
        Mj = 0;

    double foOpt = 0;

    if (t > 4)
    {
        for (i=0; i < t; i++)
        {
            j = i+2;
            while (((j+1)%t) != i)
            {
                int vi  = s.vec[i].sol;
                int vi1 = s.vec[(i+1)%t].sol;
                int vj  = s.vec[j%t].sol;
                int vj1 = s.vec[(j+1)%t].sol;

                foOpt = - dist[vi][vi1]
                        - dist[vj][vj1]
                        + dist[vi][vj]
                        + dist[vi1][vj1];

                if (foOpt < 0)
                {
                    // first improvement strategy
                    Mi1 = (i+1)%t;
                    Mj  = j%t;

                    int inicio = Mi1,
                    fim = Mj;

                    int tam, p1, p2, aux;

                    if(inicio > fim)
                        tam = t - inicio + fim + 1;
                    else
                        tam = fim - inicio + 1;

                    p1=inicio;
                    p2=fim;

                    for(int k=0; k < tam/2; k++)
                    {
                        aux = s.vec[p1%t].sol;
                        s.vec[p1%t].sol = s.vec[p2%t].sol;
                        s.vec[p2%t].sol = aux;

                        p1 = (p1==t-1)?0:p1+1;
                        p2 = (p2 == 0)?t-1:p2-1;
                    }
                    s.ofv = s.ofv + foOpt;
                }
                j++;
            }//while
        }//for
    }//if t > 4
}

void LS2(TSol &s) // NodeInsertion
{
    int i = 0,
        j = 0;
    
    double foOpt;

    for (i=0; i < n; i++)
    {
        j = (i+1)%n;
        while ( ((j+1)%n) != i )
        {
            int vi  = s.vec[i].sol;
            int viP = s.vec[(i+1)%n].sol;
            int viM = 0;
            if (i == 0)
                viM = s.vec[n-1].sol;
            else
                viM = s.vec[i-1].sol;
            
            int vp  = s.vec[j%n].sol;
            int vq = s.vec[(j+1)%n].sol;

            foOpt = - dist[vp][vq]
                    - dist[viM][vi]
                    - dist[vi][viP]
                    + dist[vp][vi]
                    + dist[vi][vq]
                    + dist[viM][viP];

            if (foOpt < 0)
            {
                // first improvement strategy
                TVecSol aux;
                
                aux.sol = s.vec[i].sol;
                aux.rk = s.vec[i].rk;
                
                s.vec.insert(s.vec.begin()+((j+1)%n), aux);
                
                if (i < ((j+1)%n))
                    s.vec.erase(s.vec.begin()+i);
                else
                    s.vec.erase(s.vec.begin()+(i+1));

                s.ofv = s.ofv + foOpt; 
            }
            j++;
        }
    }
}

void LS3(TSol &s) // NodeExchange
{
    int i = 0,
        j = 0;
    
    double foOpt;

    for (i=0; i < n-1; i++)
    {
        j = i+2;
        while (j < n)
        {
            if (i != 0 && j != n-1) //no exchange edge of the tour
            {
                int vi  = s.vec[i].sol;
                int viP = s.vec[i+1].sol;
                int viM = 0;
                if (i == 0)
                    viM = s.vec[n-1].sol;
                else
                    viM = s.vec[i-1].sol;
                
                int vj  = s.vec[j].sol;
                int vjM = s.vec[j-1].sol;
                int vjP = 0;
                if (j < n-1)
                    vjP = s.vec[j+1].sol;
                else
                    vjP = s.vec[0].sol;

                foOpt = - dist[viM][vi]
                        - dist[vi][viP]
                        - dist[vjM][vj]
                        - dist[vj][vjP]
                        + dist[viM][vj]
                        + dist[vj][viP]
                        + dist[vjM][vi]
                        + dist[vi][vjP];

                if (foOpt < 0)
                {
                    // first improvement strategy
                    TVecSol aux;

                    // exchange i and j
                    aux = s.vec[i];
                    s.vec[i] = s.vec[j];
                    s.vec[j] = aux;
                    
                    s.ofv = s.ofv + foOpt; 
                }
            }
            j++;
        }
    }
}

void LS4(TSol &s) // OrOpt2
{
    int i = 0,
        j = 0;
    
    double foOpt;

    for (i=0; i < n-1; i++)
    {
        j = i+2;
        while ( ((j+1)%n) != i )
        {
            int vi  = s.vec[i].sol;
            int viP1 = s.vec[(i+1)%n].sol;
            int viP2 = s.vec[(i+2)%n].sol;
            int viM = 0;
            if (i == 0)
                viM = s.vec[n-1].sol;
            else
                viM = s.vec[i-1].sol;
            
            int vp  = s.vec[j%n].sol;
            int vq = s.vec[(j+1)%n].sol;

            foOpt = - dist[vp][vq]
                    - dist[viM][vi]
                    - dist[viP1][viP2]
                    + dist[vp][vi]
                    + dist[viP1][vq]
                    + dist[viM][viP2];

            if (foOpt < 0) 
            {
                // first improvement strategy
                TVecSol aux1, aux2;
                
                aux1 = s.vec[i];
                aux2 = s.vec[(i+1)%n];

                if (i < ((j+1)%n))
                {
                    if (j%n == n-1)
                    {
                        s.vec.push_back(aux1);
                        s.vec.push_back(aux2);
                    }
                    else
                    {
                        s.vec.insert(s.vec.begin()+((j+1)%n), aux2); // add vi+1
                        s.vec.insert(s.vec.begin()+((j+1)%n), aux1); // add vi
                    }

                    s.vec.erase(s.vec.begin()+((i+1)%n));
                    s.vec.erase(s.vec.begin()+i);
                }
                else if (((j+1)%n) < i)
                {
                    s.vec.erase(s.vec.begin()+((i+1)%n));   // drop vi+1
                    s.vec.erase(s.vec.begin()+(i));         // drop vi

                    s.vec.insert(s.vec.begin()+((j+1)%n), aux2); // add vi+1
                    s.vec.insert(s.vec.begin()+((j+1)%n), aux1); // add vi
                }

                //update fitness
                s.ofv = s.ofv + foOpt; 
            }
            j++;
        }
    }
}

void FreeMemoryProblem()
{
    //specific problem
    dist.clear();
    node.clear();
}