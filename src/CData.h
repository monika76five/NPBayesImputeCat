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

/*
	CData class declaration.
	This function holds the data for the computations.
	It includes accessors and a function for reading the data from files
*/

#if !defined(_CDATA_H)
#define _CDATA_H
#include <vector>
class CData {
public:
	CData(); //Constructor
	~CData(); //destructor

	//data storage. These should be private, but this is going to be more efficient.
	int n; //sample size
	int J;
	int L; //max #of levels (for dimentions)
	int **x; //responses vector
	
	int *levelsJ; // (levelsJ[j] = #levels of variable j)
	int *cumLevelsJ; 

	int** ZeroMC_IJ;
	int nZeroMC;	

	int M; //number of patterns
	int *cell_countM; //count of the cell.
	
	
	void Process_MC();
	void SetData(int *x_flat, int J, int n, int* ZeroMC_flat, int nZeroMC, int *levels);
	void SetData(std::vector<int>& x_flat, int J, int n, std::vector<int>& ZeroMC_flat, int nZeroMC, std::vector<int>& levels);
	void UpdateX(std::vector<int>& x_flat); // only for Nicole Dalzell 
};

#endif  //_CDATA_H
