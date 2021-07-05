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
#include<set>
#include<algorithm>
#include<vector>
#include<list>
#include "margin_conditions.h"
#include "CArrayND.h"

using namespace std;
template<typename T_Vector>
inline
void get_relevant_cols(const T_Vector &condition, std::vector<int> &lista, const std::vector<int> &irrelevant, const int &length){
	std::set<int> EX(irrelevant.begin(), irrelevant.end());
	lista.clear();
	for (int i = 0; i < length; i++){
		if (condition[i] != -1 && EX.count(i) ==0){
			lista.push_back(i);
		}
	}
}

template<typename T>
inline
int next_arr(T dat, const std::vector<int> &seq, const int* &levelsJ, const int &order = 0){
	int j = order;
	int nCol = seq.size();
	bool carry;
	do{
		carry = false;
		dat[seq[j]]++;
		if(dat[seq[j]] >= levelsJ[seq[j]]) { 
			dat[seq[j]] = 0; 
			j ++;
			carry =true;
		}
	} while(carry && j < nCol);
	return (j);
}

template<typename T_Container>
void get_seq_exploration(const T_Container &excl_vec, std::vector<int> &seq, const int &J, const std::vector<int> &irrelevant){
	std::set<int> already(irrelevant.begin(), irrelevant.end());
	seq.clear();
	//int n = excl_vec.size();
	for (typename T_Container::const_iterator it = excl_vec.begin(); it != excl_vec.end(); it++){
		for (int j = 0; j < J; j++){
			if (((*it)[j]!= -1) && (already.count(j) == 0)){
				already.insert(j);
				seq.push_back(j);
			}
		}
	}
} 

void expand_lite(std::list<std::vector<int> > &excl_list, const int* curr_cond, const int *levels, const int &J){
	std::list<std::vector<int> > compar_list;
	vector<int> curr_cond_v(curr_cond, curr_cond + J);
	//1) test one by one and only keep those that are not disjoint.
	for (list<vector<int> >::iterator it = excl_list.begin(); it != excl_list.end(); it++){
		if (!check_disjoint(it->begin(), it->end(), curr_cond)){
			compar_list.push_back(*it);
		}
	}
	//Now compar_list is a vector of pointers to the working conditions (vectors).
	// if it's empty, the condition can enter as is
	if (compar_list.size() == 0){
		excl_list.push_back(curr_cond_v);
		return;
	}
	vector<int> cols; //sequence of column expansion
	vector<int> fixed_cols; //fixed columns
	get_relevant_cols(curr_cond, fixed_cols, cols, J);
	get_seq_exploration(compar_list, cols, J, fixed_cols);
	//special case: if cols is empty, then 
	if (cols.size() == 0){
		excl_list.push_back(curr_cond_v);
		return;
	}
	vector<int> current_expanded(curr_cond, curr_cond + J);
	for(vector<int>::iterator it = cols.begin(); it != cols.end(); it++){
		current_expanded[*it] = 0;
	}
	//bool done = false;
	unsigned int order;
	do{
		//check that is disjoint with *every* element in the list
		int  disj = 1;
		for (list<vector<int> >::iterator it = compar_list.begin(); it != compar_list.end(); it++){
			disj *= check_disjoint(it->begin(), it->end(), current_expanded.begin());
		}
		if (disj) excl_list.push_back(current_expanded);
		order = next_arr(&(current_expanded[0]), cols, levels);
	}while(order < cols.size());
}
typedef std::pair<double, int*> mypair;
bool comparator (const mypair& l, const mypair& r){
	return (l.first < r.first);
}

int** MC_to_MCpartition( int** mc, const int&J, const int &nCond, const int* levels, int &SizeDest){
	list<vector<int> > exclusive;
	//reorder according to set size.
	std::vector<mypair> sizes(nCond);
	for (int i = 0; i < nCond; i++){
		sizes[i].second = mc[i];
		sizes[i].first = 1.0;
		for (int j = 0; j < J; j++){
			if (mc[i][j] != -1) sizes[i].first *= double(levels[j]);
		}
	}
	std::sort(sizes.begin(), sizes.end(), comparator);
	exclusive.push_back(vector<int>(sizes[0].second, sizes[0].second + J));
	for (int i = 1; i < nCond; i++){
		expand_lite(exclusive, sizes[i].second, levels, J);
	}
	int C = exclusive.size();
	int **it_dest, **dest;


	it_dest = dest = (int**) CArrayND<int>::flat2arrayND(new int[J*C], sizeof(int),2,C, J);


	for(list<vector<int> >::iterator it = exclusive.begin(); it != exclusive.end(); it++){
		copy(it->begin(), it->end(), (it_dest++)[0]);
	}
	SizeDest = C;
	return(dest);
}

bool check_disjoint_MC(int** cond, int nCond, int nVars){
	for (int** it1 = cond; it1 != cond + nCond - 1; it1++){
		for (int** it2 = it1 + 1; it2 != cond + nCond; it2++){
			if (!check_disjoint(*it1, *it1 + nVars, *it2)){
				return(false);
			}
		}
	}
	return (true);
}

template<typename itT_A, typename itT_B>
inline bool check_condA_incl_condB(const itT_A &A_begin, const itT_A &A_end, const itT_B &B_begin){
	// A \ninc B iff a!=-1 and b!=a OR A=-1 and B!=-1
	itT_B itb = B_begin;
	for(itT_A ita = A_begin; ita != A_end; ita++, itb++){
		if((*ita != -1 && *itb != *ita) || (*ita == -1 && itb != -1)) return(false);
	}
	return(true);
}


/*
void test_MC_conversion(){
	int ln[] = {4,10};
	int data[] = {
		 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1, 1, 1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1, 1, 1,
	};
	int levels[] ={10, 10, 2, 6, 5, 5, 3, 10, 2, 2};

	int** arr = (int**)CArrayND<int>::flat2arrayND(data, sizeof(int), 2, ln[0],ln[1]);
	
	bool disj = check_disjoint_MC(arr, ln[0], ln[1]);
	//cout << "disjoint = " << (disj ? "yes" : "no") << endl ;
	int ndest;
	int** final = MC_to_MCpartition(arr, ln[1], ln[0], levels, ndest);
	//cout << "number of partitions " << ndest << endl;
	disj = check_disjoint_MC(final, ndest, ln[1]);
	//cout << "disjoint = " << (disj ? "yes" : "no") << endl ;
}
*/

//int main(){
//	test_MC_conversion();
//	getchar();
//}
