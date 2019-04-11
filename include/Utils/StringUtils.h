//*****************************************************
//*** AsFem: a simple finite element method program ***
//*** Copyright (C) 2018 walkandthinker             ***
//*** Contact: walkandthinker@gmail.com             ***
//*****************************************************
//* This file is part of the ASFEM framework
//* All rights reserved, see COPYRIGHT for full restrictions
//* Licensed under GPL 3.0, please see LICENSE for details
//******************************************************
//
// Created by walkandthinker on 16.08.18.
// several utils for string operation

#ifndef ASFEM_STRINGUTILS_H
#define ASFEM_STRINGUTILS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

string RemoveSpace(string instr);
string StrToLower(string instr);
vector<string> StrVecToLower(vector<string> instrvec);

string StrToUpper(string instr);
vector<string> StrVecToUpper(vector<string> instrvec);

vector<string> SplitStr(string instr,char symbol);

vector<double> SplitNum(string instr); //
vector<double> SplitNumAfter(string instr,int pos);

bool IsBracketMatch(ifstream &in,string &bracketstr,int &startline);
bool IsSubBracketMatch(ifstream &in,string &bracketstr,int &startline);

string SplitStrFromBracket(string &inputstr);

void GotoLine(ifstream &in,int linenum);

#endif //ASFEM_STRINGUTILS_H
