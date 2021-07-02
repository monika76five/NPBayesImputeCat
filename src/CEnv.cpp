/*
 * Copyright (C) 2014 Quanli Wang
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
 */ 
#include <Rcpp.h>
#include <R_ext/Utils.h>
#include <set>
#include "CData.h"
#include "CLcm.h"
#include "CTrace.h"
#include "CEnv.h"

typedef std::vector<int> intvec;

/**
* Do the actual check for an interrupt.
* @attention This method should never be called directly.
* @param[in] dummy Dummy argument.
*/
static inline void check_interrupt_impl(void* /*dummy*/) {
    R_CheckUserInterrupt();
}

/**
* Call this method to check for user interrupts.
* This is based on the results of a discussion on the
* R-devel mailing list, suggested by Simon Urbanek.
* @attention This method must not be called by any other
* thread than the master thread. If called from within
* an OpenMP parallel for loop, make sure to check
* for omp_get_thread_num()==0 before calling this method!
* @return True, if a user interrupt has been detected.
*/
inline bool check_interrupt() {
    return (R_ToplevelExec(check_interrupt_impl, NULL) == FALSE); }

CEnv::CEnv(Rcpp::IntegerMatrix x_, Rcpp::IntegerMatrix mcz_, int K, int Nmis_max,  
           double a_alpha, double b_alpha, int seed) {
	SetData(x_, mcz_);
	SetModel(K, Nmis_max, a_alpha, b_alpha, seed);
	mnburnin = 0;
	mniters = 0;
	mncurrnetburnin = 0;
	mncurrentiter = 0;
	mnsaved = 0;
	mnthinning = 1;
	t = new CTrace(m);
	mbEnableTracer = true;
	NmisOverflow = 0;
	mbsilent = true;
}

CEnv::~CEnv(void)
{
}

void CEnv::EnableTracer(bool enable) {
	mbEnableTracer= enable;
}

bool CEnv::GetTracerStatus() {
	return mbEnableTracer;
}

void CEnv::Update() {
	m -> Update();
	if (m->NmisOverflow >0) {
		NmisOverflow++;
		Rprintf( "Warning %d: maximum Nmis exceeded.\n",NmisOverflow);
		if (NmisOverflow >= 100) {
			//Rprintf("Error: Maximum Nmis has been exceeded too many times. Abort.");
			throw interrupt_exception("Maximum Nmis has been exceeded too many times. Abort.");
		}
	}
}

void CEnv::Initialize(int nWarming) {
  if (m->par->N_mis_max > 0) {
    Rprintf( "Run model with structural zeros.\n");
    m->Initializes(nWarming);
  } else {
    Rprintf( "Run model without structural zeros.\n");
    m->Initializes_no_MCZ();
  }
}

Rcpp::List CEnv::GetTrace() {
	Rcpp::List result = Rcpp::List::create();
	std::vector< std::string > list_ = t->GetTracedList();
	for (unsigned int p = 0; p < list_.size(); p++) {
		if (list_[p] == "index") {
			result("index") = std::vector<double>(t->trace[p],t->trace[p]+t->mnindex);
		}
		if (list_[p] == "alpha") {
			result("alpha") = std::vector<double>(t->trace[p],t->trace[p]+t->mnindex);
		}
		if (list_[p] == "k_star") {
			result("k_star") = std::vector<double>(t->trace[p],t->trace[p]+t->mnindex);
		}
		if (list_[p] == "Nmis") {
			result("Nmis") = std::vector<double>(t->trace[p],t->trace[p]+t->mnindex);
		}
		if (list_[p] == "nu") {
			Rcpp::NumericMatrix NU(t->mnindex, m->par->K);
			for ( int i = 0; i < t->mnindex; i++) {
				for (int j = 0; j < m->par->K; j++) {
					NU(i,j) = t->trace[p][i*m->par->K+ j];
				}
			}
			result("nu") = NU;
		}
		if (list_[p] == "z") {
			Rcpp::IntegerMatrix Z(t->mnindex, m->par->n);
			for ( int i = 0; i < t->mnindex; i++) {
				for (int j = 0; j < m->par->n; j++) {
					Z(i,j) = (int)t->trace[p][i*m->par->n+ j];
				}
			}
			result("z") = Z;
		}
		if (list_[p] == "ImputedX") {
			int size = m->par->n*m->par->J;
			Rcpp::IntegerMatrix ImputedX(t->mnindex,size);
			for (int i = 0; i < t->mnindex; i++) {
				for (int j = 0; j < size; j++) {
					ImputedX(i,j) = t->trace[p][i*size+ j];
				}
			}
			result("ImputedX") = ImputedX;
		}
		if (list_[p] == "psi") {
			int size = m->par->K * m->par->cumLevelsJ[m->par->J];
      
			Rcpp::Dimension d(m->par->J*m->par->K*m->par->L,t->mnindex);
			Rcpp::NumericVector psi4d(d);

			int size_new = m->par->J*m->par->K*m->par->L;
			for (int s = 0; s < t->mnindex; s++) {
        for (int k = 0; k < m->par->K; k++) {
  			  for (int  j= 0; j < m->par->J; j++) {
					  for (int l = 0; l < m->par->levelsJ[j]; l++) {
						  psi4d[l + k * m->par->L + j * m->par->L * m->par->K + s*size_new] = t->trace[p][size * s + (m->par->cumLevelsJ[j]+l) * m->par->K + k];
					  }
					  for (int  l = m->par->levelsJ[j]; l <m->par->L;l++) {
						  psi4d[l + k * m->par->L + j * m->par->L * m->par->K + s*size_new] = NA_REAL;
					  }
				  }
			  }
			}
			
			Rcpp::NumericVector d4(4);
			d4[0] = m->par->L; d4[1] = m->par->K; d4[2]= m->par->J; d4[3] =t->mnindex;
			psi4d.attr("dim") = d4;
			result("psi") = psi4d;
		}
	}
	return result;
}
Rcpp::IntegerMatrix CEnv::GetMCZ() {
	Rcpp::IntegerMatrix MCZ(data->nZeroMC, data->J);
	for (int i = 0; i < data->nZeroMC; i++) {
		for (int j = 0; j < data->J; j++) {
			MCZ(i,j) = data->ZeroMC_IJ[i][j];
		}
	}
	return MCZ;
}
Rcpp::List CEnv::GetParameters(std::vector< std::string > list_) {
	Rcpp::List result = Rcpp::List::create();
	for (unsigned int p = 0; p < list_.size(); p++) {
		if (list_[p] == "alpha") {
			result("alpha") = m->par->alpha;
		}
		if (list_[p] == "k_star") {
			result("k_star") = m->par->k_star;
		}
		if (list_[p] == "Nmis") {
			result("Nmis") = m->par->Nmis;
		}
		if (list_[p] == "nu") {
			result("nu") = std::vector<double>(m->par->nuK,m->par->nuK+m->par->K);
		}
		if (list_[p] == "z") {
			result("z") = std::vector<int>(m->par->zI,m->par->zI+ m->par->n);
		}
		if (list_[p] == "ImputedX") {
			Rcpp::IntegerMatrix ImputedX(m->par->n, m->par->J);
			for (int i = 0; i < m->par->n; i++) {
				for (int j = 0; j < m->par->J; j++) {
					ImputedX(i,j) = m->par->xIJ[i][j];
				}
			}
			//result("ImputedX") = mX_df;
			result("ImputedX") = ImputedX;
		}
		if (list_[p] == "psi") {
			Rcpp::Dimension d(m->par->L,m->par->K,m->par->J);
			Rcpp::NumericVector psi3d(d);
			for (int k = 0; k < m->par->K; k++) {
				for (int  j= 0; j < m->par->J; j++) {
					for (int l = 0; l < m->par->levelsJ[j]; l++) {
						psi3d[l + k * m->par->L + j * m->par->L * m->par->K] = m->par->psiJKL[m->par->cumLevelsJ[j]+l][k];
					}
					for (int  l = m->par->levelsJ[j]; l <m->par->L;l++) {
						psi3d[l + k * m->par->L + j * m->par->L * m->par->K] = NA_REAL;
					}
				}
			}
			result("psi") = psi3d;
		}
	}
	return result;
}
int CEnv::GetCurrentIter() {
	return mncurrentiter;
}
Rcpp::List CEnv::GetData() {
	return GetParameters(t->GetParameterList());
}

void CEnv::Resume() {
	if (mniters == 0) {
		Rprintf( "Run method has to be called first. Ignored.\n");
		return;
	}
	if (mncurrnetburnin < mnburnin) {
		Rprintf( "Resuming burnin at %d\n", mncurrnetburnin);
		for (; mncurrnetburnin < mnburnin; mncurrnetburnin++) {
			Update();
			if (check_interrupt()) {
				throw interrupt_exception("The burnin stage was interrupted.");
			}
		}
	}
	if (mncurrentiter < mniters) {
		Rprintf( "Resuming mcmc at %d\n", mncurrentiter);
		for (; mncurrentiter < mniters; mncurrentiter++) {
			Update();
		  if (!mbsilent) {
			    Rprintf( "iter = %d  kstar = %d alpha = %g Nmis = %d\n", mncurrentiter, m->par->k_star,m->par->alpha, m->par->Nmis ) ;
		  }
			if ( mbEnableTracer && (mncurrentiter+1) % mnthinning == 0) {
				if (t->Trace(mnsaved,mncurrentiter)) {
					mnsaved++;
				} else {
					Rprintf( "Tracer is full.\n") ;
				}
			}
			if (check_interrupt()) {
				throw interrupt_exception("The mcmc iteration was interrupted.");
			}
		}
	} else{
		Rprintf( "The last run was finished.\n") ;
	}
}
void CEnv::Run(int burnin, int iter, int thining, bool silent) {
	mnburnin = burnin;
  mbsilent = silent;
	if (mniters  == 0) { //only do Initialize once if mniters was not set
    Rprintf( "Initializing...\n");
    if (burnin == 1) {
      Initialize(1); // for short testing     
    } else {
      Initialize(500);    
    }
		
		t->PrepareTrace();
    if (!mbsilent) {
		  Rprintf( "iter = %d  kstar = %d alpha = %g Nmis = %d\n", mniters, m->par->k_star,m->par->alpha, m->par->Nmis ) ;
    }
		mnsaved = 0; //just for consistance 
	} else {
	  if (!mbsilent) {
		  Rprintf( "Continuing MCMC from previous run(s)...\n");
	  }
	}
	mniters = mncurrentiter + iter;

	//prepare tracing
	if (thining <=0) {	thining = 1; }
	mnthinning = thining;
	
	for (mncurrnetburnin = 0; mncurrnetburnin < mnburnin; mncurrnetburnin++) {
		Update();
		if (check_interrupt()) {
			throw interrupt_exception("The burnin stage was interrupted.");
		}
	}
	
	for (; mncurrentiter < mniters; mncurrentiter++) {
		Update();
	  if (!mbsilent) {
		  Rprintf( "iter = %d  kstar = %d alpha = %g Nmis = %d\n", mncurrentiter, m->par->k_star,m->par->alpha, m->par->Nmis ) ;
	  }
		if ( mbEnableTracer && (mncurrentiter+1) % mnthinning == 0) {
			if (t->Trace(mnsaved,mncurrentiter)) {
				mnsaved++;
			} else {
				Rprintf( "Tracer is full.\n") ;
			}
		}
		if (check_interrupt()) {
			throw interrupt_exception("The mcmc iteration was interrupted.");
		}
	}
}

void CEnv::Iterate(int iters) {
	for (int i = 0; i < iters; i++){
		Update();
	  if (!mbsilent) {
		    Rprintf( "iter = %d  kstar = %d alpha = %g Nmis = %d\n", i, m->par->k_star,m->par->alpha, m->par->Nmis ) ;
	  }
	}
}

void CEnv::SetData(Rcpp::IntegerMatrix x_, Rcpp::IntegerMatrix mcz_) {
	int J = x_.nrow();
	int n = x_.ncol();
	int nZeroMC = mcz_.ncol();
	if (mcz_.nrow() != J) { // no mcz
		nZeroMC = 0; 
	}
  	
	intvec x = Rcpp::as<intvec>(x_);
	intvec mcz = Rcpp::as<intvec>(mcz_);
  
	Rcpp::List something(x_.attr("levels"));
	intvec levels(something.size());
	for (unsigned int i = 0; i < levels.size(); i++) {
		SEXP exp = something[i];
		Rcpp::CharacterVector v(exp);
		levels[i] = v.size();
		//Rprintf( "%d\n", levels[i]) ;
	}
	SetData(x, J, n, mcz, nZeroMC, levels);	
}

void CEnv::SetData(std::vector<int>& x_flat, int J, int n, std::vector<int>& ZeroMC_flat, int nZeroMC, std::vector<int>& levels) {
	data = new CData;
	data->SetData(x_flat, J, n, ZeroMC_flat, nZeroMC, levels);
}

void CEnv::SetModel(int K, int Nmis_max,  double a_alpha, double b_alpha, int seed) {
	m = new CLcm(data, K, Nmis_max,a_alpha,b_alpha, seed);
}

void CEnv::SetTrace(std::vector< std::string > list_, int size) {
	t -> SetTrace(list_,size);
	if (mniters > 0) {
		Rprintf( "Tracer has been reset.\n");
		t->PrepareTrace();
		mnsaved = 0;
	} 
}

std::vector<std::string>  CEnv::traced(){
	return t->GetTracedList();
}

std::vector<std::string>  CEnv::traceable(){
	return t->GetParameterList();
}

void CEnv::SetXAsDataframe(Rcpp::DataFrame X_df) {
	mX_df = X_df;
}

void CEnv::UpdateX(Rcpp::IntegerMatrix x_) {
  intvec x = Rcpp::as<intvec>(x_);
  data->UpdateX(x);
  m->par->UpdateX(data,m->mt);
  //Rprintf( "Data updated.\n");
}
