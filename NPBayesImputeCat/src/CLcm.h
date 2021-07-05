/*
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
 * Modified by Quanli Wang, 2014, 2019
 */ 
 
#ifndef _CLcm_H
#define _CLcm_H

#include <cstdio>
#include <cstdlib>
#include "CParam.h"
#include "SpecialFunctions.h"

class CLcm {
public:
	//Constructors
	CLcm(CData *_data, CParam *_par) : par(_par), data(_data){
		class_construct(data, par);
	}
	CLcm(CData* _data, int K, int Nmis_max,  double a_alpha, double b_alpha, int seed)
			: randomseed(seed), data(_data) {
				par = new CParam(data->J, K, data->L, data->levelsJ, data->cumLevelsJ,data->n, Nmis_max, data->ZeroMC_IJ, data->nZeroMC, a_alpha, b_alpha,data->x);
		class_construct(data, par);
	}
	void Update();
	
	double predict(int* xJ);
	double predict_renorm(int* xJ);
	void Initializes(int nwarming); //a more careful initialization than the defalt.
  void Initializes_no_MCZ();
	
	CParam *par;
	CParam *par_temp; //temp solution, only to be used during initialization

	MTRand mt;
	int NmisOverflow;
	int randomseed;
private:
	CData* data;
	void class_construct(CData* data, CParam* par){
		//r = dan_initRandom(); //initializes a random number generator
		//mt.seed(1234);
		mt.seed(randomseed);
		current_iteration = 0;
		burnin=1;
		thining=1;
		iterations = 99;

		K = par->K;
		J = data->J; 
		n = data->n; 

	}
protected:
	int J, n, K; //local copies. bad idea...
	//Full conditional samplers
	void sam_x(); //imputation step.
	void sam_nu(); 
	void sam_psi();
	void sam_z();
	void sam_Z2_X2();
	void sam_Nmis();
	void sam_alpha();
	void compute_probs_miss();

protected:
	int iterations;
	int current_iteration;
	int thining;
	int burnin;
};


#endif  //_CLcm_H
