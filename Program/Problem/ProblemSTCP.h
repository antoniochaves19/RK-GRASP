// *******************************************************************
//      file with specific functions to solve the STCP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Global Variables
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     > k - number of local search heuristics
extern int n;                               // size of cromossoms

//---------------------- DEFINITION OF TYPES OF PROBLEM SPECIFIC --------------------


//---------- DEFINITION OF CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  ------------

static std::vector <std::vector <int> > triples;	// matrix with triples
static std::vector <std::vector <int> > varCover;	// matrix with variables cover
static int numVariables;							// number of variable
static int numTriples;								// number of triples


//----------------------- IMPLEMENTATION OF FUNCTIONS  -------------------------------

/************************************************************************************
 Method: ReadData
 Description: read input data of the problem
*************************************************************************************/
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
    fscanf(arq, "%d %d", &numVariables, &numTriples);

    // read graph informations
    triples.clear();
    triples.resize(numTriples, std::vector<int>(3));

    varCover.clear();
    varCover.resize(numVariables);

    for (int i=0; i<numTriples; i++)
    {
    	fscanf(arq, "%d %d %d", &triples[i][0], &triples[i][1], &triples[i][2]);

        // corrigir o indice
        triples[i][0] = triples[i][0] - 1;
        triples[i][1] = triples[i][1] - 1;
        triples[i][2] = triples[i][2] - 1;

        // inserir a informacao que a tripla i eh atendida pela variavel k
        varCover[triples[i][0]].push_back(i);
        varCover[triples[i][1]].push_back(i);
        varCover[triples[i][2]].push_back(i);
    }

    fclose(arq);
    
    n = numVariables;

    // imprimir
    // printf("%d %d \n", numVariables, numTriples);

    // for (int i=0; i<numTriples; i++)
    // {
    // 	printf("%d %d %d\n", triples[i][0], triples[i][1], triples[i][2]);
    // }

    // for (int i=0; i<numVariables; i++)
    // {   
    //     printf("\n");
    //     for (int j=0; j<varCover[i].size(); j++)
    //     {
    //         printf("%d ", varCover[i][j]);
    //     }
    // }
    // printf("\n");

    // getchar();
}

/************************************************************************************
 Method: Decoders
 Description: users need to implement at least one decoder, DecK (K = [1,2,3,4,5])
*************************************************************************************/

void Dec1(TSol &s) // best + redudance
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

    // usar as variaveis em ordem enquanto nao cobrir todas os pontos
    int numCobertura = 0;
    std::vector<int> cover;
    cover.clear();
    cover.resize(numTriples);

    std::vector<int> usedColumn;
    usedColumn.clear();
    usedColumn.resize(n,0);

    int k = 0;              // coluna corrente
    int numVarUsed = 0;     // numero de colunas utilizadas
    int sumCover = 0;       // soma do numero de linhas cobertas (contando sobreposicao)
    while (numCobertura < numTriples)
    {
        int used = 0;

        // para cada linha coberta pela coluna k
        for (int i=0; i<(int)varCover[s.vec[k].sol].size(); i++)
        {
            // verificar se a linha ainda nao foi coberta
            if (cover[varCover[s.vec[k].sol][i]] == 0)
            {
                // atualizar a cobertura da linha
                cover[varCover[s.vec[k].sol][i]] = 1;
                numCobertura++;

                used = 1;  
            }
            else
            {
                // atualizar o numero de vezes que a linha eh coberta
                cover[varCover[s.vec[k].sol][i]]++;
                used = 1;
                sumCover++;
            }
        }

        // se a coluna k foi usada, aumentar o numero de colunas 'abertas'
        if (used == 1){
            numVarUsed++;     
            usedColumn[s.vec[k].sol] = 1;       
        }
        k++;
    }

    // verificar se posso retirar alguma coluna redundante
    for (int i=0; i<k; i++)
    {
        int redundante = 1;
        for (int j=0; j<(int)varCover[s.vec[i].sol].size(); j++)
        {
            if (cover[varCover[s.vec[i].sol][j]] == 1)
            {
                redundante = 0;
                break;
            }
        }

        // remover a coluna i
        if (redundante)
        {
            numVarUsed--;
            usedColumn[s.vec[i].sol] = 0;       

            // atualizar o vetor de cobertura
            for (int j=0; j<(int)varCover[s.vec[i].sol].size(); j++)
            {
                cover[varCover[s.vec[i].sol][j]]--;
                sumCover--;
            }
        }
    }

    // set solution vector
    for (int i=0; i<n; i++){
        s.vec[i].sol = usedColumn[i];
    }

    // set objective function value
    s.ofv = numVarUsed;
}

void Dec2(TSol &s) // best
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

    // usar as variaveis em ordem enquanto nao cobrir todas os pontos
    int numCobertura = 0;
    std::vector<int> cover;
    cover.clear();
    cover.resize(numTriples);

    int k = 0;  // coluna corrente
    int numVarUsed = 0;
    while (numCobertura < numTriples)
    {
        int used = 0;

        // para cada linha coberta pela coluna k
        for (int i=0; i<(int)varCover[s.vec[k].sol].size(); i++)
        {
            // verificar se a linha ainda nao foi coberta
            if (cover[varCover[s.vec[k].sol][i]] == 0)
            {
                // atualizar a cobertura da linha
                cover[varCover[s.vec[k].sol][i]] = 1;
                numCobertura++;

                used = 1;  
            }
        }

        // se a coluna k foi usada, aumentar o numero de colunas 'abertas'
        if (used == 1){
            numVarUsed++;            
        }
        k++;
    }

    s.ofv = numVarUsed;
}

void Dec3(TSol &s) // binary representation
{
    // create a initial solution of the problem
    s.ofv = 0;

    int numVarUsed = 0;             // numero de colunas utilizadas
    for (int j = 0; j < n; j++)
	{
        if (s.vec[j].rk >= 0.5)
        {
		    s.vec[j].sol = 1;
            numVarUsed++;
        }
        else
            s.vec[j].sol = 0;
	}

    int numCobertura = 0;           // numero de linhas cobertas

    std::vector<int> cover;
    cover.clear();
    cover.resize(numTriples);

    // verificar se a solucao eh viavel
    for (int k=0; k<n; k++)
    {
        if (s.vec[k].sol == 1)
        {
            // para cada linha coberta pela coluna k
            for (int i=0; i<(int)varCover[k].size(); i++)
            {
                // verificar se a linha ainda nao foi coberta
                if (cover[varCover[k][i]] == 0)
                {
                    // atualizar a cobertura da linha
                    cover[varCover[k][i]] = 1;
                    numCobertura++;
                }
                else
                {
                    // atualizar o numero de vezes que a linha eh coberta
                    cover[varCover[k][i]]++;
                }
            }
        }
    }

    // solucao nao eh viavel
    while (numCobertura < numTriples)
    {
        // encontrar a coluna k nao aberta que cobre o maior numero de linhas nao cobertas
        int bestCover = -1;
        int bestK = 0;
        for (int k=0; k<n; k++)
        {
            if (s.vec[k].sol == 0)
            {
                int numCoverK = 0;
                 // para cada linha coberta pela coluna k
                for (int i=0; i<(int)varCover[k].size(); i++)
                {
                    // verificar se a linha ainda nao foi coberta
                    if (cover[varCover[k][i]] == 0)
                    {
                        // somar o numero de linhas novas cobertas
                        numCoverK++;
                    }
                }

                if (numCoverK > bestCover)
                {
                    bestCover = numCoverK;
                    bestK = k;
                }
            }
        }

        // inserir a coluna bestK na solucao   
        s.vec[bestK].sol = 1;    
        numVarUsed++;   

        // atualizar as linhas cobertas
        for (int i=0; i<(int)varCover[bestK].size(); i++)
        {
            // verificar se a linha ainda nao foi coberta
            if (cover[varCover[bestK][i]] == 0)
            {
                // atualizar a cobertura da linha
                cover[varCover[bestK][i]] = 1;
                numCobertura++;
            }
            else
            {
                // atualizar o numero de vezes que a linha eh coberta
                cover[varCover[bestK][i]]++;
            }
        }
        // printf("\nPreso aqui: %d - %d", numCobertura, numVarUsed);
        // getchar();
    }

    // verificar se posso retirar alguma coluna redundante
    for (int i=0; i<n; i++)
    {
        if (s.vec[i].sol == 1)
        {
            int redundante = 1;
            for (int j=0; j<(int)varCover[i].size(); j++)
            {
                if (cover[varCover[i][j]] == 1)
                {
                    redundante = 0;
                    break;
                }
            }

            // remover a coluna i
            if (redundante)
            {
                numVarUsed--;
                s.vec[i].sol = 0;       

                // atualizar o vetor de cobertura
                for (int j=0; j<(int)varCover[i].size(); j++)
                {
                    cover[varCover[i][j]]--;
                }
            }
        }
    }


    // encode
    for (int i=0; i<n; i++){
        if (s.vec[i].sol == 1 && s.vec[i].rk < 0.5)
            s.vec[i].rk = 1 - s.vec[i].rk;

        if (s.vec[i].sol == 0 && s.vec[i].rk >= 0.5)
            s.vec[i].rk = 1 - s.vec[i].rk;
    }

    // set objective function value
    s.ofv = numVarUsed;

    // printf("\n%lf",s.ofv);
}
void Dec4(TSol &s){}
void Dec5(TSol &s){}

/************************************************************************************
 Method: Local Search Heuristics
 Description: users need to implement local searchs heuristics if he/she use the RVND,
 LS_K (K = [1,2,3,4,5])
*************************************************************************************/

void LS1(TSol &s) 
{
    TSol sViz = s;
    // trocar a random-key de uma coluna usada por uma não usada
    for (int i = 0; i < n-1; i++)
    {
        for (int j = i+1; j < n; j++)
        {
            if (s.vec[i].sol != s.vec[j].sol)
            {
                // trocar as random keys das colunas i e j
                double temp = sViz.vec[i].rk;
                sViz.vec[i].rk = sViz.vec[j].rk;
                sViz.vec[j].rk = temp;

                Decoder(sViz);

                // atualiza a solucao corrente se melhorou e continua a busca a partir de sViz
                if (sViz.ofv < s.ofv){
                    s = sViz;
                }
                // retorna para a solucao corrente
                else{
                    sViz = s;
                    // temp = sViz.vec[hi].rk;
                    // sViz.vec[hi].rk = sViz.vec[pj].rk;
                    // sViz.vec[pj].rk = temp;
                }
            }
        }
    }  
}

void LS2(TSol &s){}
void LS3(TSol &s){}
void LS4(TSol &s){}
void LS5(TSol &s){}

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem()
{
    //specific problem
    triples.clear();
}

#endif