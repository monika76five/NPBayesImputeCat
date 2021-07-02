#include "CData.h"
#include "CParam.h"
#include "margin_conditions.h"
CParam::~CParam() {
  delete aux_dirCumJK_ND;
  delete  psiJKL_ND;
  delete xIJ_ND;
  if (nZeroMC > 0) {
    delete MCZ_ND;
    delete x2_NMax_J_ND;
  }
  
}

void CParam::class_construct(int J, int K, int L, int *levelsJ, int *cumLevelsJ, int n){
	//init members	
	this->J=J; this->K=K; this->n = n; 
	this->levelsJ = levelsJ;
	this->cumLevelsJ = cumLevelsJ;
	this->L = L;

	//allocate
	zI = new int[n];
	nuK = new double[K];
	log_nuK = new double[K];
	countK = new int[K];
  
  aux_dirCumJK_ND = CArrayND<int>::CreateArray(2, cumLevelsJ[J], K);
  aux_dirCumJK = (int **)aux_dirCumJK_ND->data;
  psiJKL_ND = CArrayND<double>::CreateArray(2, cumLevelsJ[J], K);
  psiJKL = (double **)psiJKL_ND->data;
}

void CParam::class_construct(int Nmis_max, int** MCZ_, int nZeroMC, int **x){
	//this assumes that base class has been initialized
	//initialize members.
	this->nZeroMC = nZeroMC;
	if (nZeroMC > 0) {
		//this->Nmis = Nmis;
		this->N_mis_max = Nmis_max;
		//allocate extra.
		pZeroMC_I = new double[nZeroMC];
		z2_Nmax = new int[Nmis_max];
		count_partition = new unsigned int[nZeroMC];
    
    MCZ_ND = CArrayND<int>::CreateArray(2, this->nZeroMC, J);
    MCZ = (int **)MCZ_ND->data;
    x2_NMax_J_ND =CArrayND<int>::CreateArray(2, Nmis_max, J);
    x2_NMax_J = (int**)x2_NMax_J_ND->data;
		//copy MCZ
		std::copy(MCZ_[0], MCZ_[0] + this->nZeroMC * J, this->MCZ[0]); 
	} else {
		this->Nmis = 0;
		this->N_mis_max = 0;
	}
	xIJ_ND = CArrayND<int>::CreateArray(2, n, J);
	xIJ = (int **)xIJ_ND->data; 
	//fill x with data
	std::copy(x[0], x[0] + this->n * J, xIJ[0]); //NOTE THAT HERE X INCLUDES MISSING VALUES!!! HAVE TO INITIALIZE.
	
}

void CParam::UpdateX(CData* dat,MTRand & mt) {
  //fill x with data
  std::copy(dat->x[0], dat->x[0] + n * J, xIJ[0]); 
  
  //Need to add the code below to reinitlize xIJ, etc when there are misisng values and MCZ
  //Discuss with Nicole later
  if (this->nZeroMC > 0)  {
    std::vector<double> p(this->L); //reserve at least the maximum number of levels.
    std::fill(p.begin(), p.end(), 1.0);
    //Initialize values for xIJ.
    bool badvalue;
    for(int i = 0; i < n; ++i){
      std::vector<int> x_working(xIJ[i], xIJ[i]+J);
      do{
        for (int j = 0; j < J; ++j){
          if (xIJ[i][j] == -1) {
            x_working[j] =SpecialFunctions::discreterand(this->levelsJ[j], &(p[0]),mt);
          }
        }
        //check
        badvalue = false;
        if (this->nZeroMC > 0) {
          for (int ii  = 0; ii < this->nZeroMC; ++ii){
            if (!check_x_notin_mu(x_working.begin(), x_working.end(), MCZ[ii])) {
              badvalue = true; 
              break; 
            }
          }
        }
      } while ( badvalue );
      std::copy<std::vector<int>::iterator, int*>(x_working.begin(), x_working.end(), xIJ[i]);
    }  
  }
}
void CParam::initizalize(MTRand& mt){
	//Initializes the parameters for the chain.
	int j,k,l;
	for (k =0; k < K; k++){
		this->nuK[k] = 1.0 / double(K);
	}

	for(j = 0; j < J; j++){
		for(k = 0; k < K; k++){
			for(l = 0; l < this->levelsJ[j]; l++){
				psiJKL[cumLevelsJ[j]+l][k] = 1.0 / double(levelsJ[j]);
			}
		}
	}

	//just the extra.
	double logk = -log(double(K));
	for( k = 0; k < K; k++) { log_nuK[k] = logk;}
	alpha = 1 ; //initialize at full spread: uniform prior.


	std::vector<double> p(this->L); //reserve at least the maximum number of levels.
	std::fill(p.begin(), p.end(), 1.0);

	prob_zero = 0;
	Nmis = 0;
	for (int i = 0; i < this->nZeroMC; i++){
		count_partition[i] = 0;
	}
	memset(this->countK, 0, sizeof(int)*this->K);
	//Initialize values for xIJ.
	bool badvalue;
	for(int i = 0; i < n; ++i){
		std::vector<int> x_working(xIJ[i], xIJ[i]+J);
		do{
			for (int j = 0; j < J; ++j){
				if (xIJ[i][j] == -1) {
					x_working[j] =SpecialFunctions::discreterand(this->levelsJ[j], &(p[0]),mt);
				}
			}
			//check
			badvalue = false;
      if (this->nZeroMC > 0) {
  			for (int ii  = 0; ii < this->nZeroMC; ++ii){
  				if (!check_x_notin_mu(x_working.begin(), x_working.end(), MCZ[ii])) {
  					badvalue = true; 
  					break; 
  				}
  			}
      }
		} while ( badvalue );
		std::copy<std::vector<int>::iterator, int*>(x_working.begin(), x_working.end(), xIJ[i]);
	}
}

