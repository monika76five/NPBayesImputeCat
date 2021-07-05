/*
 * Copyright (C) 2007-2014 Daniel Manrique-Vallier
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 * Modified by Quanli Wang, 2014
 */ 
#include <cstring>
#include <cmath>

#include "CData.h"
#include "CLcm.h"
#include "margin_conditions.h"
#include "R.h" //for ISNAN function

//--------------------------------------------------------------------------------
// Implementation of class CLcm
//--------------------------------------------------------------------------------


void CLcm::Update() {
	//regular stuff
	sam_z();
	sam_psi();
	sam_nu();
	//truncation!
	compute_probs_miss();
	sam_Nmis();
	sam_Z2_X2();
	//and regularization!
	sam_alpha();
	sam_x();
}

// Samplers for the full conditional distributions.
//  each function update the values on the *par object directly.
void CLcm::sam_x(){
	//imputation step as gibbs sampler.
	//goes through each missing component. Optimize this precomputing the locations of missing data.
	std::vector<int> indexes;
	std::vector<double> probs(data->L); 
	for (int i = 0; i < data->n; i++){
		int k = par->zI[i];
		for (int j = 0 ; j < data->J; ++j){
			if (data->x[i][j] == -1){ //if component is missing, impute component by component, conditional on the previous draws.
				get_valid_levels_for_j(par->xIJ[i], data->ZeroMC_IJ, j, data->levelsJ[j], data->J, data->nZeroMC, indexes);
				for (unsigned int c = 0; c < indexes.size(); ++c){
					probs[c] = par->psiJKL[par->cumLevelsJ[j]+indexes[c]][k];
				}
				par->xIJ[i][j] = indexes[SpecialFunctions::discreterand(indexes.size(), &(probs[0]),mt)];
			}
		}
	}
}

inline
void CLcm::sam_psi(){
	int i, j, k, l;
	std::fill_n(par->aux_dirCumJK[0], data->cumLevelsJ[par->J] * par->K, 0);
	//the data part.
	for (i = 0; i < data->n; ++i){
		int k = par->zI[i];
		for (j = 0; j < par->J; ++j){
			++par->aux_dirCumJK[par->cumLevelsJ[j] + par->xIJ[i][j]][k];
		}
	}
	//the imputed part.
	for (i = 0; i < par->Nmis; ++i){
		k = par->z2_Nmax[i];
		for (j = 0; j < par->J; ++j){
			++par->aux_dirCumJK[par->cumLevelsJ[j] + par->x2_NMax_J[i][j]][k];
		}
	}
	
	//add the prior (+1) and sample
	for(j = 0; j < data->J; ++j){
		for(k = 0; k < par->K; ++k){
			//add prior and cast.
			double acc = 0.0;
			for(l = 0; l < data->levelsJ[j]; ++l){
				acc += par->psiJKL[par->cumLevelsJ[j]+l][k] = SpecialFunctions::gammarand_int(1 + par->aux_dirCumJK[par->cumLevelsJ[j] + l][k],mt);// +1 = uniform prior
			}
			acc = 1.0 / acc;
			//normalize
			for (l = 0; l < data->levelsJ[j]; ++l){
				par->psiJKL[par->cumLevelsJ[j]+l][k] *= acc;
			}
		}
	}
}
inline
void CLcm::sam_z(){
	//THIS IS THE CRITICAL FUNCTION IN THIS SAMPLER. for K=50 the smapler spends 50% of the time executing this function. 
	// any little optimization in the inner loop has a huge impact.
	//Here we also count how many individuals belong to each class for i <= n.
	int i, j, k;
	int K = par->K;
	double *probs = new double[K];
  if (par->Nmis == 0) { //only needed for no-strutual zeros
	  memset(par->countK, 0, sizeof(int)*par->K);// <- Now this is cleared in the X2 section. MAKE SURE OF INITIALIZING IT TO ZEROS.
  }
	//Z1
	for (i = 0; i < data->n; ++i){
		std::copy(par->nuK, par->nuK + K, probs);
    
		int *xiJ = par->xIJ[i];
		for (j = 0; j < par->J; ++j){
			int l = xiJ[j];
			for (k = 0; k < par->K; k++){
				probs[k] *= par->psiJKL[par->cumLevelsJ[j]+l][k];
			}
		}
		//sample zi
		par->zI[i] = k = SpecialFunctions::discreterand(K, probs, mt);
		par->countK[k]++;
	}
	par->k_star = K - std::count(par->countK, par->countK + K, 0);
	//We're sampling Z2 in sam_Z2_X2.
  delete [] probs;
}
inline
void CLcm::sam_Z2_X2(){
	if (par->N_mis_max ==0) { return;}
	//updates:
	//-number of elements in each partition (count_partition)
	//-all individuals' responses for i > n (X2)
	//-all individuals' membership for i >n (Z2)
	//-Number of individuals in each class for i>n (countK)
	int c, k, j;
	int K = par->K;
	double *prob = new double[K];
	
	
	//1) sample total number of individuals in each partition.
	SpecialFunctions::multinomialrand(data->nZeroMC, par->Nmis, par->pZeroMC_I, par->count_partition,mt);
	//2)sample individuals within each partition
	int index = 0; //individual counter for X2.
	std::fill(par->countK, par->countK + K, 0);
	double * temp = new double[data->L];

	for (c = 0; c < data->nZeroMC; ++c){
		//compute the *unnormalized posterior probability of being in class k given
		// that individual is in partition c.
		std::copy(par->nuK, par->nuK + K, prob);
		for (j = 0; j < data->J; ++j){
			int l = data->ZeroMC_IJ[c][j];
			for (k = 0; k < K; ++k){
				if (l != -1) prob[k] *= par->psiJKL[par->cumLevelsJ[j]+l][k];
			}
		}
		//setup a sampler for the partition's zs:
		for (unsigned int i = 0; i < par->count_partition[c]; ++i, ++index){
			//sample class
			par->z2_Nmax[index] = k = SpecialFunctions::discreterand(K,prob,mt);
			++(par->countK[k]); //increment the class counter.
			for (j = 0; j < data->J; ++j){
				if (data->ZeroMC_IJ[c][j] == -1) {
					for (int l = 0; l < data->levelsJ[j]; l++) {
						temp[l] = par->psiJKL[par->cumLevelsJ[j]+l][k];
					}
					par->x2_NMax_J[index][j] =SpecialFunctions::discreterand_norm(data->levelsJ[j], temp, 1.0, mt);
				} else {
					par->x2_NMax_J[index][j] = data->ZeroMC_IJ[c][j];
				}
			}
		}
	}
	delete [] temp;
  delete [] prob;
}

inline
void CLcm::sam_nu(){
  int k;
  double b = 0.0, a = 0.0;
  int n_acc = 0;
  double l_acc_prod = 0.0;
	for (k = 0; k < par->K - 1; ++k){
		  n_acc += par->countK[k];
		  a = double(1 + par->countK[k]);
		  b = par->alpha + double(par->Nmis + data->n - n_acc);
		  double lgamma1 = SpecialFunctions::log_gamma_rand(a,mt);
		  double lgamma2 = SpecialFunctions::log_gamma_rand(b,mt);
		  double lsumgamma = SpecialFunctions::log_sum(lgamma1, lgamma2);
		  par->log_nuK[k] = (lgamma1 - lsumgamma + l_acc_prod);
      if (ISNAN(par->log_nuK[k])) {
        par->log_nuK[k] = -50;
      }
      if (par->log_nuK[k] < -50) par->log_nuK[k] = -50;//approx 1e-22
      
		  l_acc_prod += lgamma2 - lsumgamma; //log(1-V)
		  par->nuK[k] = exp(par->log_nuK[k]);
	}
	par->log_nuK[par->K -1] = l_acc_prod;
  if (ISNAN(par->log_nuK[par->K-1])) {
        par->log_nuK[par->K-1] = -50;
  }
  if (par->log_nuK[par->K-1] < -50) par->log_nuK[par->K - 1] = -50;//approx 1e-44
    
  par->nuK[par->K - 1] = exp(par->log_nuK[par->K-1]);
}

/*
inline
void CLcm::sam_nu(){
  int k;
	double b = 0.0, a = 0.0, lgamma1, lgamma2, lsumgamma;
	int n_acc = 0;
	double l_acc_prod = 0.0;
	for (k = 0; k < par->K - 1; ++k){
		  n_acc += par->countK[k];
		  a = double(1 + par->countK[k]);
		  b = par->alpha + double(par->Nmis + data->n - n_acc);
		  lgamma1 = SpecialFunctions::log_gamma_rand(a,mt);
		  lgamma2 = SpecialFunctions::log_gamma_rand(b,mt);
		  lsumgamma = SpecialFunctions::log_sum(lgamma1, lgamma2);
		  par->log_nuK[k] = (lgamma1 - lsumgamma + l_acc_prod);
      if (par->log_nuK[k] < -50) par->log_nuK[k] = -50;//approx 1e-22
		  l_acc_prod += lgamma2 - lsumgamma; //log(1-V)
		  par->nuK[k] = exp(par->log_nuK[k]);
	}
	par->log_nuK[par->K -1] = l_acc_prod;
  if (par->log_nuK[par->K-1] < -50) par->log_nuK[par->K - 1] = -50;//approx 1e-44
	par->nuK[par->K - 1] = exp(par->log_nuK[par->K-1]);
}
*/

inline
void CLcm::sam_alpha(){
	par->alpha = SpecialFunctions::gammarand(par->a_alpha + double(par->K) - 1.0, 1.0 / (par->b_alpha - par->log_nuK[par->K-1]),mt);
}

typedef std::pair<double, int> _mypair;
bool _comparator (const _mypair& l, const _mypair& r){
	return (l.first < r.first);
}

void CLcm::Initializes_no_MCZ(){
  
	par->initizalize(mt);
	par->alpha = 1.0;
  
  for (int k =0; k < par->K; k++){
		par->nuK[k] = 1.0 / double(K);
	}

	for (int j =0; j < par->J; j++){
		for (int l = 0; l < par->levelsJ[j]; l ++) 
			for (int k = 0; k < par->K; k++)
				par->psiJKL[par->cumLevelsJ[j]+l][k] = 1.0 / double(par->levelsJ[j]);
	}
  sam_z();
  sam_psi();
	sam_nu();
}

void CLcm::Initializes(int nwarming){
	double g1, g2, gs;
	par->initizalize(mt);
	par->alpha = 1.0;
  NmisOverflow = 0;
	double lacc = 0.0;
	for (int k = 0; k < par->K - 1; k++){
		g1 = SpecialFunctions::log_gamma_rand(1.0,mt);
		g2 = SpecialFunctions::log_gamma_rand(1.0,mt);
		gs = SpecialFunctions::log_sum(g1, g2);
		par->log_nuK[k] = g1 - gs + lacc;
		lacc += g2 - gs;
	}
	par->log_nuK[par->K -1] = lacc;
	par->nuK[par->K-1] = exp(lacc);
	std::vector<std::vector<int> > v(par->J);
	std::vector<int> denom(par->J);
	int den;
	for (int j =0; j < par->J; j++){
		den = 1;
		for( int l = 0; l < par->levelsJ[j]; l++) v[j].push_back(1);
		for (int i = 0; i < par->n; i++){
			if (data->x[i][j] != -1){
				++v[j][ data->x[i][j] ];
				++den;
			}
		}
		for (int l = 0; l < par->levelsJ[j]; l ++) 
			for (int k = 0; k < par->K; k++)
				par->psiJKL[par->cumLevelsJ[j]+l][k] = double(v[j][l]) / double(den);
	}
	par->alpha = 0.3 / double(K);
	for(int k =0; k < 4; k++){
		par->nuK[k] = 0.25 * 0.9;
	}
	for(int k = 4; k < par->K; k++){
		par->nuK[k] = 0.1 / double (par->K -4);
	}
	for (int i = 0; i < nwarming; ++i){
		sam_z();
		sam_psi();
		sam_nu();
		//truncation!
		compute_probs_miss();
		sam_Nmis();
		sam_Z2_X2();
	}
  
	//reorder decreasingly according to nuK. 
	std::vector<_mypair> s_index(par->K);
	for (int k = 0; k < par->K; k++){
		s_index[k].first = 1.0 - par->nuK[k];
		s_index[k].second = k;
	}
	std::sort(s_index.begin(), s_index.end(), _comparator);


	par_temp = new CParam(data, par->K, par->N_mis_max, par->a_alpha, par->b_alpha);
	for (int k = 0; k < par->K; ++k){
		int indx = s_index[k].second;
		par_temp->nuK[k] = par->nuK[indx];
		par_temp->log_nuK[k] = par->log_nuK[indx];
		par_temp->countK[k] = par->countK[indx];
	}
	for (int k = 0; k < par->K; ++k){
		par->nuK[k] = par_temp->nuK[k];
		par->log_nuK[k] = par_temp->log_nuK[k];
		par->countK[k] = par_temp->countK[k];
	}

	//switch columns (indexed by k]
	for (int k = 0; k < par->K; ++k){
		int index = s_index[k].second;
		for (int j = 0; j < par->cumLevelsJ[J]; j++){
			par_temp->psiJKL[j][k] = par->psiJKL[j][index];
		}
	} 

	//copy tmp to par
	
	
	std::copy(par_temp->psiJKL[0],par_temp->psiJKL[0] + par->cumLevelsJ[J] * K , par->psiJKL[0]);
	for (int i = 0;  i < data->n; ++i) {par_temp->zI[i] =0;}
	for (int i = 0;  i < par->Nmis; ++i) {par_temp->z2_Nmax[i] =0;}

	for (int k = 0; k < par->K; ++k){
		for (int i = 0; i < data->n; ++i){
			if (k == s_index[k].second) par_temp->zI[i] = k;
		}
		for (int i = 0; i < par->Nmis; ++i){
			if (k == s_index[k].second) par_temp->z2_Nmax[i] = k;
		}
	}
	for (int i = 0; i < data->n; ++i){
			par->zI[i] = par_temp->zI[i];
		}
	for (int i = 0; i < par->Nmis; ++i){
		par->z2_Nmax[i] = par_temp->z2_Nmax[i];
	}
	
}



inline
void CLcm::compute_probs_miss(){
	if (par->nZeroMC == 0) { return;}
	
	//computes the probability of each partition and the aggregate
	int c, k, j;
	par->prob_zero = 0.0;
	std::fill_n(par->pZeroMC_I, par->nZeroMC, 0.0);

	for (c = 0; c < data->nZeroMC; ++c){
		for (k = 0; k < par->K; ++k){
			double prod = par->nuK[k];
			for (j = 0; j < data->J; ++j){
				int l = data->ZeroMC_IJ[c][j];
				if (l != -1) prod *= par->psiJKL[par->cumLevelsJ[j]+l][k];
			}
			par->pZeroMC_I[c] += prod;
		}
		par->prob_zero += par->pZeroMC_I[c];
	}
}

inline
void CLcm::sam_Nmis(){
	NmisOverflow = 0;
	if (par->N_mis_max == 0) {
		par->Nmis = 0;
	} else {
		int tries = 1000;
		par->Nmis = SpecialFunctions::negative_binomial_rand(1.0 - par->prob_zero, data->n,mt);
		int count = 0;
		while (count < tries && par->Nmis >= par->N_mis_max) { //try a thousand times
			par->Nmis = SpecialFunctions::negative_binomial_rand(1.0 - par->prob_zero, data->n,mt);
			count++;
		}
		if (count==tries) { //truncate here
			par->Nmis = par->N_mis_max - 1;
		}
		if (count > 0) {
				NmisOverflow = 1;
		}
	}
}

