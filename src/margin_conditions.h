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
 */ 

#ifndef _MARGIN_CONDITIONS_H
#define _MARGIN_CONDITIONS_H
#include <vector>
#include <list>

int** MC_to_MCpartition( int** mc, const int&J, const int &nCond, const int* levels, int &SizeDest);
bool check_disjoint_MC(int** cond, int nCond, int nVars);

template<typename a_iterator, typename b_iterator>
inline bool check_disjoint(const a_iterator &a_start, const a_iterator &a_end, const b_iterator &b_start){
	bool res = false;
	b_iterator it_b = b_start;
	for (a_iterator it_a = a_start; it_a != a_end; it_a++, it_b++){
		if (*it_a != -1 && *it_b != -1 && *it_a != *it_b){
			res = true;
			break;	
		}
	}
	return(res);
}
template <typename T_Iterator>
inline
std::vector<int> MC_Intersection(const T_Iterator &X_start, const T_Iterator &X_end, const T_Iterator &Y_start){
	std::vector<int> res;
	for (T_Iterator it1 = X_start, it2 = Y_start; it1 != X_end; ++it1, ++it2){
		if (*it1 == *it2) res.push_back( *it1 );
		else if(*it1 == -1) res.push_back (*it2);
		else if(*it2 == -1) res.push_back (*it1);
		else {
			res.clear();
			break;
		}
	}
	return(res);
}

template <typename iterator_x, typename iterator_mu>
inline
bool check_x_notin_mu (iterator_x x_start, iterator_x x_end, iterator_mu mu_start){
	//checks if an *incomplete* vector x can be determined to be contained into the cells
	//defined by a margin condition.
	//Note that for x, -1 means NA. For mu, -1 means any level.
	bool notin = false;
	iterator_mu itmu = mu_start;
	for (iterator_x itx = x_start; itx != x_end; ++itx, ++itmu){
		if (*itmu != -1 && *itx != *itmu ){
			notin = true;
			break;
		}
	}
	return(notin);
}

inline
void get_valid_levels_for_j(int* x, int** MCZs, int j, int levels_J, int J, int n_mus, std::vector<int> &result){
	//returns which levels of component j of x would avoid x belonging to any of the margin conditions in MCZs.
	result.clear();
	bool ok;
	std::vector<int> v(x, x+J);
	for (int l = 0; l < levels_J; ++l){
		v[j] = l;
		ok = true;
		for (int i = 0; i < n_mus; ++i){
			if (!check_x_notin_mu(v.begin(), v.end(), MCZs[i])){
				ok = false;
				break; //being in one mu is enough.
			}
		}
		if (ok) result.push_back(l);
	}
}



#endif
