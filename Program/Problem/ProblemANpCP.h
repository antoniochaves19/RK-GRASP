// *******************************************************************
//      file with specific functions to solve the alpha-NpCP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Variables declared in main.cpp
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     > k - number of local search heuristics
extern int n;                               // size of cromossoms

//------ DEFINITION OF TYPES OF PROBLEM SPECIFIC --------


//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------

static std::vector <std::vector <int> > dist;	// matrix with Euclidean distance

static int numNodes;								// number of nodes
static int numMedians;								// number of medians
static int alphaN;									// number of alpha neighbors




//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------


// Sort TSol by random-keys
// bool sortByRk(const TVecSol &lhs, const TVecSol &rhs) { return lhs.rk < rhs.rk; }

// Sort facilities in terms of distances
// bool sortBydist(const nodeinfo1 &lhs, const nodeinfo1 &rhs) { return lhs.xdistance < rhs.xdistance;}

static int pcenter=0, alphap=2,nNodes=0, nAux = 0, Narcs0 = 0;


void ReadData(char nameTable[], int &n, int& nNode)
{ 

 //  std::cout << " enter 1" << std::endl;
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


// read node informations
    
   // int nAux = 0,Narcs1=0,Narcs0=0;
    //node.clear();
   // TNode nodeTemp;
    static int  Narcs1 = 0;
   // read instance head
  //  char temp[100];
   // 6linhas ignoradas no tsp if (fgets(temp, sizeof(temp), arq) == NULL) exit(1);;
   //read first instance line
    if (fscanf(arq, "%d %d %d", &nAux, &Narcs0, &pcenter) != 3) {
        // Handle error, e.g., by setting default values or reporting the issue
    }
   //  std::cout<< nAux<< " , "<< Narcs0 << " , "<< pcenter<<std::endl;
 //pcenter=10;
  //  alphap=5;
    // calculate the euclidean distance
   // dist.clear();
  //  dist.resize(nAux, std::vector<double>(nAux,10000000));
    for (int ii = 0; ii < nAux; ii++) {
        dist.push_back(std::vector<int>(nAux, 1000000));
        dist[ii][ii] = 0;
    }
   
  int nodexx=0,nodexy=0,nodexc=0;
  // while (!feof(arq))
  while (fscanf(arq, "%d %d %d", &nodexx, &nodexy, &nodexc) == 3)
    {  
        //int nodexx=0,nodexy=0,nodexc=0;
        //fscanf(arq, "%d %d %d", &nodexx, &nodexy, &nodexc);
    	//if (fscanf(arq, "%d %d %d", &nodexx, &nodexy, &nodexc) == 0) 
       //{
          //  break;
        //}
        // std::cout<< nodexx<< " , "<< nodexy << " , "<< nodexc<<std::endl;
            dist[nodexx-1][nodexy-1]=nodexc;
            dist[nodexy-1][nodexx-1]=nodexc;
          //cout<< dist[nodex-1][nodey-1] << "  ";
	
       // Narcs1++;
      // if (Narcs1==Narcs0)
          //  break;
    }
    fclose(arq);
 // std::cout << " out 1" << std::endl;
// completa o grafo usando o Algoritmo

    for (int kk = 0; kk < nAux; kk++) {
       
        // Pick all vertices as source one by one
        for ( int ii = 0; ii < nAux; ii++) {
            // Pick all vertices as destination for the
            // above picked source
            for (int jj = 0; jj < nAux; jj++) {
                // If vertex k is on the shortest path from
                // i to j, then update the value of
                // dist[i][j]
                if (dist[ii][jj] > (dist[ii][kk] + dist[kk][jj]))
                    dist[ii][jj] = dist[ii][kk] + dist[kk][jj];

                
            }
           
        }
       
    }
    
    int maior_dist=0;
    for (int ii = 0; ii < nAux; ii++) {
     //   std::cout << "  indice " << ii << "  ";
        for (int jj = 0; jj < nAux; jj++) {
            if(dist[ii][jj]>maior_dist)
              maior_dist=dist[ii][jj];
           // std::cout << dist[ii][jj] << "  ";
        }
      //  std::cout << std::endl;
    }
  //std::cout << std::endl<< " dsd "<< maior_dist<<std::endl;
     v_dists.resize(maior_dist+1,1);

        for (int ii=0; ii<nAux; ii++)
    {
   
        std::vector<nodeinfo1> aux01;
            nodeinfo2 aux03;
    	for (int jj=0; jj<nAux; jj++)
    	{
           // cout<< dist[i][j] << "  ";
           // if(jj!=ii){
               nodeinfo1 aux02;
               aux02.node_s=jj;
               aux02.xdistance = dist[ii][jj];
               v_dists[dist[ii][jj]]++;
               aux01.push_back(aux02);
           // }
    	}
        sort(aux01.begin(), aux01.end(), sortBydist);
        
        for (int jj=0; jj<nAux; jj++)
    	{
            
                    aux03.dist_s.push_back(aux01[jj]);
                  
            
    	}
        
       //  cout<<endl<<endl;
        vec_min.push_back(aux03);
    }
  
 
    n =pcenter;
    nNodes=nAux;
    nNode = nNodes;
}

void Decoder(TSol &s, int n, int nDec)
{
    // copy the random-key sequence of current solution 
   // TSol temp = s;

    // define decoder function based in the random-key of position n+1
    int dec = floor(s.vec[n].rk*nDec)+1;
    //std::cout<< " dec 11  "<<std::endl;
    if(dec==1){
       
            Dec1(s, n);  
           //  std::cout<< " dec 22  "<<std::endl;       
    }
     
    // return initial random-key sequence and maintain the solution sequence
  //  for (int i=0; i<n; i++){
       // s.vec[i].rk = temp.vec[i].rk;
   // }
    
}

void LocalSearch(TSol &s, int n, int nLS)
{
    // ***** we use a Random Variable Neighborhood Descent (RVND) as local search ****

    // current neighborhood
	int k = 1;

    // predefined number of neighborhood moves
    std::vector <int> NSL;
    std::vector <int> NSLAux;
    
    for (int i=1; i<=nLS; i++)
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
                LS1(s, n); 
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
}

double CalculateFitness(TSol& s, int n)
{
    double sofv=0;
    s.desempate = 0;
   //s.ofv = 0;
     std::vector <int> d_dist;
     std::vector <bool> check_facility;
     check_facility.resize(nNodes, 0);
     for (int ip = 0; ip < n; ip++)
         check_facility[s.vec[ip].sol] = 1;
     d_dist.resize(v_dists.size(), 0);
     int Max_dist = 0, max_index = -1;
    // for each node find the alpha nearest medians 
    for (int j=0; j<nNodes; j++)
    {
        std::vector<int> medians(n,0);
        int minDist;
        int minIndex;
        

        // find the k nearest median
         //bool check_facility=0;
                // for(int ip=0;ip<n;ip++)
                  //  if(s.vec[ip].sol==j)
                     //  check_facility=1;
      if(check_facility[j] == 0) {
        for (int k=0; k<alphap; k++)
        {
            minDist = 99999999;
            minIndex = 0;
            Max_dist=0;
    

            for (int i=0; i<n; i++)
            {
                if ( (dist[s.vec[i].sol][j] < minDist) && (medians[i] == 0) )
                {
                    minDist = dist[s.vec[i].sol][j];
                    minIndex = i;
                }
            }

            // insert the alpha nearest median in the solution
            // solution[j][k] = s.vec[minIndex].sol;

            medians[minIndex] = 1;
            d_dist[minDist]++;
            if (minDist > Max_dist) {
                Max_dist = minDist;
  
            }
               //sofv=std::max(sofv,minDist);
          //  s.ofv += minDist;

            // printf("\n");
            // for (int i=0; i<numMedians; i++){
            //     printf("%d \t", medians[i]);
            // }
        }
            if(Max_dist>sofv)
                sofv=Max_dist;
       }
        // getchar();
    }
    // desempate
    double delta_chagas=0;
     for(int ii=0;ii<=sofv;ii++){
        if(d_dist[ii]>0){
            delta_chagas= double(delta_chagas+ d_dist[ii])/ (v_dists[ii]);
           // std::cout <<" delta "<< delta_chagas << std::endl;
        }
     }
    // if(sofv<25)
       // std::cout <<" objetivo fit  "<< sofv<<" vertices extremos " << d_dist[sofv] << std::endl;

        // s.desempate = delta_chagas;
     s.desempate = d_dist[sofv];

     s.ofv = sofv;

    return s.ofv;
}

void Dec1(TSol &s, int n) // sort
{

     // create an initial solution of the problem
    std::vector<int> nodesList;

    // list of candidate nodes
    for (int j = 0; j < nNodes; j++){
        nodesList.push_back(j);
	}

    // define the p medians
    for (int i=0; i<n; i++)
    {
        // define the index of the nodeList
        int index = floor(s.vec[i].rk * (int)nodesList.size());

        // define the node that is a median
        s.vec[i].sol = nodesList[index];

        // remove this node of the nodeList
        nodesList.erase(nodesList.begin()+index);
    }
     

   // s.ofv = CalculateFitness(s,n);
     s.ofv=CalculateFitness(s, n);

   //  std::cout << " dec 222  " << s.ofv<< std::endl;
}
void LS1(TSol& s, int n) // 2-Opt
{


    // std::cout << "Run 1 " << std::endl;
    int runs = int(s.desempate);
    double OBJ1 = s.ofv, OBJ2 = 0;
   // std::cout << std::endl << " objetivo Run " << s.ofv << " vertices extremos " << s.desempate << std::endl;
    for (int ite = 0; ite < runs; ite++) {
         OBJ2 = OBJ1;
      //  std::cout << std::endl << "iteracao " << ite << std::endl;
        // inicializar o vetor de nos
        std::vector<double> nodes;
       
        std::vector<int> inversoes;
        double foOpt = 0;
        TSol sViz = s;

        std::vector <nodeinfo2> vec_solutions;					// vector of TSP nodes
        vec_solutions.resize(n);
        inversoes.resize(n, 0);

        std::vector <int> d_dist;
        d_dist.resize(v_dists.size(), 0);

        std::vector <bool> check_facility;
        check_facility.resize(nNodes, 0);
        for (int ip = 0; ip < n; ip++)
            check_facility[s.vec[ip].sol] = 1;

        int sofv = 0;
        int Max_dist = 0, max_index = -1;
        int mindex = 0, tail_node = 0;
        // for each node find the alpha nearest medians 
        for (int j = 0; j < nNodes; j++)
        {
            std::vector<int> medians(n, 0);
            int minDist;
            int minIndex;


            // find the k nearest median
         //   bool check_facility = 0;
           // for (int ip = 0; ip < n; ip++)
              //  if (s.vec[ip].sol == j)
                  //  check_facility = 1;
           if (check_facility[j] == 0) {


               // nodes.push_back(j);

                for (int k = 0; k < alphap; k++)
                {
                    minDist = 99999999;
                    minIndex = 0;
                    // Max_dist = 0;


                    for (int i = 0; i < n; i++)
                    {
                        if ((dist[s.vec[i].sol][j] < minDist) && (medians[i] == 0))
                        {
                            minDist = dist[s.vec[i].sol][j];
                            minIndex = i;
                        }
                    }

                    // insert the alpha nearest median in the solution
                    // solution[j][k] = s.vec[minIndex].sol;
                    nodeinfo1 solution_aux;
                    solution_aux.node_s = j;
                    vec_solutions[minIndex].dist_s.push_back(solution_aux);
                    medians[minIndex] = 1;
                    d_dist[minDist]++;

                    if (minDist > Max_dist) {
                        Max_dist = minDist;
                        max_index = minIndex;
                        tail_node = j;
                    }
                    //sofv=std::max(sofv,minDist);
               //  s.ofv += minDist;

                 // printf("\n");
                 // for (int i=0; i<numMedians; i++){
                 //     printf("%d \t", medians[i]);
                 // }
                }
               // if (Max_dist > sofv)
                   sofv = Max_dist;
            }
            // getchar();
        }
        OBJ1 = sofv;

       // std::cout << std::endl << " fitness  " << OBJ2<< " "<< OBJ1<< " novos extremos "<<d_dist[OBJ1]<< std::endl;
        for (int ip = 0; ip < n; ip++) {
            nodeinfo1 solution_aux;
            solution_aux.node_s = tail_node;
            vec_solutions[ip].dist_s.push_back(solution_aux);
        }
        mindex = max_index;
        for (int ip = 0; ip < vec_min[tail_node].dist_s.size(); ip++)
            if (vec_min[tail_node].dist_s[ip].xdistance < sofv && check_facility[vec_min[tail_node].dist_s[ip].node_s]==0)
                nodes.push_back(vec_min[tail_node].dist_s[ip].node_s);
     //   std::cout << std::endl;
      //   for (int i = 0; i < n; i++)
      //   {
           //  std::cout << s.vec[i].sol << " ";
            // for (int ii = 0; ii < vec_solutions[i].dist_s.size(); ii++)
                // std::cout << " node " << vec_solutions[i].dist_s[ii].node_s << " dist " << dist[s.vec[i].sol][vec_solutions[i].dist_s[ii].node_s] << " || ";
      //  }
         // adiciona sol[max_index] ao vetor de distancias
         //nodeinfo1 solution_aux;
         //solution_aux.node_s = s.vec[max_index].sol;
          // vec_solutions[max_index].dist_s.push_back(solution_aux);
       //  std::cout << std::endl << " vertice extremo  " << s.vec[max_index].sol << " max tail " << tail_node << std::endl;
        //find candidate solution
        int melhorou = 0, search_cont = 0, tentativa =0;
        while (nodes.size() > 0)
        {
          //  std::cout << std::endl << " err  " << std::endl;
           // std::cout << std::endl << " node size " << nodes.size() << std::endl;
            int tentativa = 0;
            int K = irandomico2(0, nodes.size() - 1);

            // define the index of the nodeList

           // std::cout << std::endl;
            // define the cadidate node to open
              int node_i = -1;
            node_i = nodes[K];
           // std::cout << nodes.size() << " " << nodes[K] << " || ";

            //find facility to close
            std::vector<int> facilidades1(n,0);

            for (int it = 0; it < n; it++)
                facilidades1[it] = it;

            while (facilidades1.size() > 0) {
                search_cont = 0;
              //  std::cout << std::endl << " facility size " << facilidades1.size() << std::endl;
                int KK = 0;
                if (tentativa == 0) {
                    KK = mindex;
                    tentativa = 1;
                }
                else {
                    KK = irandomico2(0, facilidades1.size() - 1);
                   
                }
                max_index = KK;

                std::vector<int> medians(n, 0);
                int minDist;
                int minIndex;

               
                for (int k = 0; k < alphap; k++)
                {
                    minDist = 99999999;
                    minIndex = 0;
                   

                    for (int i = 0; i < max_index; i++)
                    {
                     
                        if ((dist[s.vec[i].sol][s.vec[max_index].sol] < minDist) && (medians[i] == 0))
                        {
                            minDist = dist[s.vec[i].sol][s.vec[max_index].sol];
                            minIndex = i;
                        }
                    }

                    for (int i = max_index; i < n; i++)
                    {
                         
                        if ((dist[s.vec[i].sol][s.vec[max_index].sol] < minDist) && (medians[i] == 0))
                        {
                            minDist = dist[s.vec[i].sol][s.vec[max_index].sol];
                            minIndex = i;
                        }
                    }

                    // insert the alpha nearest median in the solution
                    medians[minIndex] = 1;

                    if (minDist < s.ofv) {
                        search_cont++;
                    }

                }

                medians.clear();

               
                if (search_cont >= alphap) {

                    for (int ii = 0; ii < vec_solutions[max_index].dist_s.size(); ii++) {

                        int node_j = vec_solutions[max_index].dist_s[ii].node_s;
                        //  if (node_j != node_i) {
                        melhorou = 0;
                        if (dist[node_i][node_j] >= sofv) {
                            int nearestalpha = 0;
                            for (int jj = 0; jj < max_index; jj++)
                                if (dist[s.vec[jj].sol][node_j] < sofv) {
                                    nearestalpha++;
                                    //melhorou = 1;
                                   // std::cout << ii << " | " << dist[node_i][node_j] << " || " << dist[s.vec[jj].sol][node_j] << " |||| ";
                                    //jj = max_index;
                                }
                            // if (melhorou == 0) {
                            for (int jj = max_index + 1; jj < n; jj++)
                                if (dist[s.vec[jj].sol][node_j] < sofv) {
                                    nearestalpha++;
                                    // melhorou = 1;
                                   //  std::cout << ii << " | " << dist[node_i][node_j] << " || " << dist[s.vec[jj].sol][node_j] << " |||| ";
                                    // jj = n;
                                }
                            //}
 
                                if (nearestalpha >= alphap)
                                    melhorou = 1;
                            
                        }
                        else {
                            melhorou = 1;
                           //   std::cout << ii << " | " << dist[node_i][node_j] << " |||| ";
                        }

                        if (melhorou == 0)
                            ii = vec_solutions[max_index].dist_s.size(); // terminate loop process
                        // }
                    }

                    // std::cout << "node " << node_i << " " << melhorou << std::endl;
                    if (melhorou == 0) {
                        // remove this node of the nodeList
                    //    std::cout << "Removing facility at index: " << KK << " from facilidades1" << std::endl;
                        facilidades1.erase(facilidades1.begin() + KK);
                     //   std::cout << "Facility removed, new size: " << facilidades1.size() << std::endl;
                    }
                    else
                        break;

                    // if (nodes.size() == 0)
                      //   break;
                }
                else {
                   // std::cout << "Removing facility at index: " << KK << " from facilidades1" << std::endl;
                    facilidades1.erase(facilidades1.begin() + KK);
                   // std::cout << "Facility removed, new size: " << facilidades1.size() << std::endl;
                }
            }
            //  std::cout << std::endl << " saiu com " << melhorou << "  " << vec_solutions[max_index].dist_s.size() << std::endl;
           
            // create inversion vector
            if (melhorou == 1) {

                //  std::cout << std::endl << s.vec[max_index].sol << " trocou por " << node_i << " antiga fo " << s.ofv << "  max dist " << sofv << std::endl;
                  // for (int ii = 0; ii < n; ii++)
                     //  std::cout << s.vec[ii].sol << "  ";
               // std::cout << std::endl << " no lugar de " << s.vec[max_index].sol << std::endl;
                s.vec[max_index].sol = node_i;
              //  std::cout << std::endl << " adicionar" << node_i << std::endl;
                for (int ii = max_index; ii < n; ii++) {
                    for (int tt = 0; tt < ii; tt++)
                        if (s.vec[ii].sol > s.vec[tt].sol) {
                            inversoes[ii]++;
                        }
                     //std::cout << ii<< " sol "<<s.vec[ii].sol << " " << inversoes[ii] << std::endl;
                }
                for (int ii = max_index; ii < n; ii++) {
                    s.vec[ii].rk = (double(s.vec[ii].sol + 1 - inversoes[ii]) / double(nNodes + 1 - ii));
                    // if(s.vec[ii].rk>1)
                      //   std::cout << s.vec[ii].rk << "  "<< s.vec[ii].sol<<" "<< inversoes[ii] << " " << nNodes << " " << ii << std::endl;
                }
               // std::cout << std::endl << " err 1 " << std::endl;
               // s.ofv = CalculateFitness(s, n);
              //  std::cout << std::endl << " nova fo " << s.ofv << " " << ite << std::endl;
                break;
                // vec_solutions.clear();
                // vec_solutions.resize(n);
               // sofv = 0; Max_dist = 0, max_index = -1;


                // s.ofv = sofv;
                //std::cout << std::endl << " apos o encoder virou " << std::endl;
                // for (int ii = 0; ii < n; ii++)
                   //  std::cout << s.vec[ii].sol << "  ";


               // std::cout << std::endl << " nova fo " << s.ofv << " " << s.vec[max_index].sol << std::endl;

                // for (int i = 0; i < n; i++)
                 //{
                    // std::cout<<std::endl << " facility " << s.vec[i].sol << std::endl;
                    // for (int ii = 0; ii < vec_solutions[i].dist_s.size(); ii++)
                     //    std::cout << " node " << vec_solutions[i].dist_s[ii].node_s << " dist " << dist[s.vec[i].sol][vec_solutions[i].dist_s[ii].node_s] << " || ";
                 //}

            }
            else
            {
              //  std::cout << std::endl << " err 2 " << std::endl;
              //  std::cout << "Removing node at index: " << K << " from nodes" << std::endl;
                nodes.erase(nodes.begin() + K);
              //  std::cout << "Node removed, new size: " << nodes.size() << std::endl;
              //  std::cout << std::endl << " err 3 " << std::endl;
            }
            //facilidades1.clear();
           // std::cout << std::endl << " err 4 " << std::endl;
           
        }

       // std::cout << std::endl << " err 5 " << std::endl;

        vec_solutions.clear();
        nodes.clear();
        inversoes.clear();
        check_facility.clear();
        
    }
   
    s.ofv = CalculateFitness(s, n);
   // std::cout << std::endl << " nova fitness " << s.ofv<< std::endl;
    
   // std::cout << std::endl << " Out " << std::endl;

}

void FreeMemoryProblem()
{
    std::cout << " enter 11" << std::endl;
    //specific problem
    dist.clear();
    v_dists.clear();
    node.clear();
    vec_min.clear();
    //pcenter = 0; alphap = 3; nNodes = 0; nAux = 0; Narcs0 = 0;
    std::cout << " out 11 " << std::endl;
}
