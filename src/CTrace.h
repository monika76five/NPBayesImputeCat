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
#if !defined(_CTRACE_H)
#define _CTRACE_H
class CTrace
{
public:
	CTrace(CLcm* pm);
	~CTrace(void);
	
	void SetTrace(std::vector<std::string> list, int size);
	std::vector<std::string> GetParameterList();
	std::vector<std::string> GetTracedList();
	void PrepareTrace();
	bool Trace(int index,int currentiter);
	void ClearTrace();
	
	double ** trace;
	
	int mnindex;
	int mnsize;
private:
	CLcm*  m;
	std::vector<std::string> TracedParameters;
};
#endif  //_CENV_H
