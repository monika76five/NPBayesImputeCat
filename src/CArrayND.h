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

#ifndef _CARRAYND_H
#define _CARRAYND_H

#include <cstdarg>
#include <iostream>
#include <cstdlib>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <string>
#include <cstdio>
#include <algorithm>
#include <cstring>

template<typename T>
class CArrayND{
public:
	int dims; //number of dimensions
	int lengths[20]; //size of each dimension. MAX = 20.
	void* data; //pointer to the data structure (THE ARRAY OF POINTERS). HAS TO BE RECASTED IN THE CALLING CLASS.
	void* data_base; //the base address of the data (the actual memory obtained through malloc)
	bool allocated;
	CArrayND(int dims,  int* lengths){
		construct_class();
		int size_elem = sizeof(T);
		alloc_ln(size_elem, dims, lengths);
	}
	~CArrayND(){
		if (this->allocated) {
			free(data_base);
			if (this->dims >1) free(data);
		}
	}
	
	static CArrayND<T>* CreateArray(int dims, ...) {
		int lengths[20], i;
		va_list args;
		va_start(args, dims);
		for (i = 0; i < dims; i++){
			lengths[i] = va_arg(args, int);
		}
		va_end(args);
		return (new CArrayND<T>(dims, lengths));
	}
	
	static void *flat2arrayND(void* data, int size_elem, int dims, ...){
		//read the variable input...
		int lengths[20];
		int i;

		va_list args;
		va_start(args,dims);
		for (i = 0; i < dims; i++){
			lengths[i] = va_arg(args, int);
		}
		va_end(args);
		return flat2arrayND_ln(data, size_elem, dims, lengths);
	}

private:
	static void *flat2arrayND_ln(void* data, int size_elem, int dims, int *lengths){
		int size[20]; //max 20 dimensions (!). More would be too much, I think...
		char *indexes, *current, *next;
		int i, j, aux;
		
		if (dims == 1) return data;
		//compute the length of each level of indexes and its sum
		size[0] = lengths[0];
		aux = lengths[0];
		for (i = 1; i < dims -1; i ++){//is dims-1
			size[i] = size[i-1]*lengths[i];
			aux += size[i];
		}
		//allocate memory for the array of indexes. The last level (data) is already allocated!!
		if ((indexes =(char*)malloc(sizeof(void*)*aux)) == NULL){
			//std::cerr << "Failed to allocate memory" << std::endl;
			//exit(1);
			return NULL;
		}
		current = indexes; //current level to fill (leftmost index)
		
		//fill all levels except the one pointing directly to the data
		for (i = 0; i < dims -2; i++){
			next = current + size[i]*sizeof(void*); //next points to the next level of indexes
			//fill each element of the current level of indexes with the addresses
			for (j = 0; j < size[i]; j++){
				// of every lengths[i+1] elements of the next level.
				(*( (void**)current + j)) = next + lengths[i+1]*j*sizeof(void*);
			}
			current = next;
		}
		//fill the last level (pointers to the real data!)
		for (j = 0; j < size[dims -2]; j++){
			*( (void**)current + j) = (char*)data + lengths[dims -1]*j*size_elem;
		}
		return (void*)indexes;
	}
	void alloc_ln(int size_elem, int dims, int *lengths) {
		int size=1, i;
		if (dims < 1 || allocated) return;
		for (i = 0; i < dims; i++) size *= lengths[i];
		size *= size_elem;
		if ((this->data_base = malloc(size)) == NULL){
			//std::cerr << "Failed to allocate memory" << std::endl;
			//exit(1);
		}
		this->data = flat2arrayND_ln(data_base, size_elem, dims, lengths);
		this->dims = dims;
		memcpy(this->lengths, lengths, sizeof(int)*dims);
		this->allocated = true;
	}

	void construct_class(){		
		dims  = 0;
		data_base = data = NULL;
		allocated  = false;
	}
};

#endif
