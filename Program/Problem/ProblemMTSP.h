// *******************************************************************
//      file with specific functions to solve the MTSP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Global Variables
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     > k - number of local search heuristics
extern int n;                               // size of cromossoms


//---------------------- DEFINITION OF TYPES OF PROBLEM SPECIFIC --------------------

//estruturas para representar o grafo
struct TGrafo
{
	int task;
	int tool;
};

struct TRelacao
{
	static std::vector <std::vector <TGrafo> > relacao;
};

//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  ----------

static int numTasks;  			   // numero de tarefas
static int numTools; 			   // numero de ferramentas
static int cap;     	  		   // capacidade da maquina

//Armazena os dados dos arquivos - ferramentas e tarefas
static std::vector <std::vector <int> > dados;

static std::vector <std::vector <TRelacao> >  grafo;



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
        // getchar();
        exit(1);
    }

    // => read data

    fscanf(arq, "%d", &numTasks);
	fscanf(arq, "%d", &numTools);
	fscanf(arq, "%d", &cap);

    dados.clear();
    dados.resize(numTools, std::vector<int>(numTasks));

	for (int i=0; i<numTools; i++)
	{
		for (int j=0; j<numTasks; j++)
		{
			fscanf(arq, "%d", &dados[i][j]);
			//printf("%d ", dados[i][j]);
		}
		//printf("\n");
	}

    fclose(arq);
    
    n = numTasks;

    // imprimir
    // printf("%d %d %d\n", numTasks, numTools, cap);

    // for (int i=0; i<numTools; i++)
    // {   
    //     printf("\n");
    //     for (int j=0; j<numTasks; j++)
    //     {
    //         printf("%d ", dados[i][j]);
    //     }
    // }
    // printf("\n");

    // getchar();
}

void verificarMemoria(void *p)
{
	 if (p == NULL)
	 {
		 printf("\n>>>Memoria insuficiente para %s\n","Erro memoria");
		 exit(1);
	 }
	 return;
}

void liberarMatriz(int **matriz, int nlinhas)
{
	int i;

	for (i=nlinhas-1; i >= 0; i--)
		free((int *) matriz[i]);
	free((int *) matriz);
}

/*****************************************************************************
Metodo: CALCULAR TROCAS FERRAMENTAS
Funcao: Determinacao do numero minimo de trocas de ferramentas dada uma
		sequencia de tarefas. Baseado na politica KTNS (Keep Tool Needed
		Soonest), de Tang e Denardo, 1988.

		Njobs - numero de tarefas sendo sequenciadas  - n
		Ljobs - lista de tarefas sendo sequenciadas   - s.seq
******************************************************************************/
int calcularTrocasFerramentas(int Njobs, TSol s)
{
	 int i,j,k,m1;  //m -> m1
	 int soma,maior,tem;
	 int Ntrocas;
	 int *usada;
	 int **L;

	 int tools = numTools;

	 if (Njobs == 1) return 0;

	 // L[i][k] = primeira tarefa, a partir da k-esima tarefa (inclusive),
	 // que requer a ferramenta i

	 L = (int **)calloc(tools, sizeof(int *));
	 verificarMemoria(L);

	 for (i = 0; i < tools; i++)
	 {
		 L[i] = (int *)calloc(Njobs, sizeof(int));
		 verificarMemoria(L[i]);
	 }

	 usada = (int *)calloc(tools, sizeof(int));
	 verificarMemoria(usada);

	 for (k = 0; k < Njobs; k++)
	 {
		 for (i = 0; i < tools; i++)
		 {
			 tem = 0;
			 for (j = k; ((j < Njobs) && (tem == 0)); j++)
			 {
				 if (dados[i][s.vec[j].sol] == 1)
				 {
					 L[i][k] = j;
					 tem = 1;
				 }
			 }
			 if (tem == 0)
				 L[i][k] = Njobs + 1;
		 }
	 }

	 for (i = 0; i < tools; i++)
		 usada[i] = 0;

	 soma = 0;
	 for (j = 0; ((soma < cap) && (j < Njobs)); j++)
	 {
		 for (i = 0; ((soma < cap) && (i < tools)); i++)
		 {
			 if ((usada[i] == 0) && (L[i][0] == j))
			 {
				 usada[i] = 1;
				 soma++;
			 }
		 }
	 }

	 Ntrocas = 0;
	 for (m1 = 0; m1 < Njobs; m1++)
	 {
		 for (i = 0; i < tools; i++)
		 {
			 if ((L[i][m1] == m1) && (usada[i] == 0))
			 {
				 maior = 0;
				 for (j = 0; j < tools; j++)
				 {
					 if (usada[j] == 1)
					 {
						 if (L[j][m1] > maior)
						 {
							 k = j;
							 maior = L[j][m1];
						 }
					 }
				 }
				 usada[i] = 1;
				 usada[k] = 0;
				 Ntrocas++;
			 }
		 }
	 }
	 free(usada);
	 liberarMatriz(L,tools);

	 return Ntrocas;
}

/************************************************************************************
 Method: CalculateFitness
 Description: calculate the fitness of a chromossom s
*************************************************************************************/
double CalculateFitness(TSol s)
{
    // calculate objective function
    s.ofv = 0;
    // for (int i=0; i<n; i++){
    //     s.ofv += dist[s.vec[i%n].sol][s.vec[(i+1)%n].sol];
    // }

    return s.ofv;
}

/************************************************************************************
 Method: Decoders
 Description: users need to implement at least one decoder, DecK (K = [1,2,3,4,5])
*************************************************************************************/

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

    s.ofv = calcularTrocasFerramentas(n, s);
}

void Dec2(TSol &s){}
void Dec3(TSol &s){}
void Dec4(TSol &s){}
void Dec5(TSol &s){}

/************************************************************************************
 Method: Local Search Heuristics
 Description: users need to implement local searchs heuristics if he/she use the RVND,
 LSK (K = [1,2,3,4,5])
*************************************************************************************/

void LS1(TSol &s) 
{
	TSol sViz = s;
    // trocar as random-key de duas tarefas
    for (int i = 0; i < n-1; i++)
    {
        for (int j = i+1; j < n; j++)
        {
			// trocar as random keys das tarefas i e j
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
	dados.clear();
	grafo.clear();
}

#endif