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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 */ 
#pragma once
#if !defined(_CENV_H)
#define _CENV_H

/**
* Base class for interrupt exceptions thrown when user
* interrupts are detected.
*/
class interrupt_exception : public std::exception {
public:
    /**
	* Constructor.
	* @param[in] message A description of event that
	* caused this exception.
	*/
    interrupt_exception(std::string message): detailed_message(message){};

    /**
	* Virtual destructor. Needed to avoid "looser throw specification" errors.
	*/
    virtual ~interrupt_exception() throw() {};
    
	/**
	* Obtain a description of the exception.
	* @return Description.
	*/
    virtual const char* what() const throw() { return detailed_message.c_str();
    }
    
	/**
	* String with details on the error.
	*/
    std::string detailed_message;
};

class CEnv
{
public:
	CEnv(Rcpp::IntegerMatrix x_, Rcpp::IntegerMatrix mcz_, int K, int Nmis_max,  double a_alpha, double b_alpha, int seed);
	void SetXAsDataframe(Rcpp::DataFrame X_df);
	void UpdateX(Rcpp::IntegerMatrix x_);
	~CEnv(void);
	void Update();
	void Initialize(int nWarming);
	Rcpp::List GetTrace();
	void SetTrace(std::vector< std::string > list_, int size);
	Rcpp::List GetParameters(std::vector< std::string > list_);
	Rcpp::List GetData();
	Rcpp::IntegerMatrix GetMCZ();
	std::vector<std::string> traceable();
	std::vector<std::string> traced();
	void Run(int burnin, int iter, int thining, bool silent = true);
	void EnableTracer(bool enable);
	bool GetTracerStatus();
	void Resume();
	int GetCurrentIter();
	void Iterate(int iters);
	CLcm*  m;
	CTrace* t;
	CData* data;
	int NmisOverflow;
private:
	void SetData(Rcpp::IntegerMatrix x_, Rcpp::IntegerMatrix mcz_);
	void SetData(std::vector<int>& x_flat, int J, int n, std::vector<int>& ZeroMC_flat, int nZeroMC, std::vector<int>& levels);
	void SetModel(int K, int Nmis_max,  double a_alpha, double b_alpha, int seed);
	
	
	Rcpp::DataFrame mX_df;
	int mnburnin;
	int mniters;
	int mncurrnetburnin;
	int mncurrentiter;
	int mnsaved;
	int mnthinning;
	
	bool mbsilent;
	bool mbEnableTracer;

};
#endif  //_CENV_H
