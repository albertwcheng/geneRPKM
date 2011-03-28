/***************************************************************************
 Copyright 2011 Wu Albert Cheng <albertwcheng@gmail.com>
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 *******************************************************************************/ 


#ifndef _ADV_GET_OPT_H
#define _ADV_GET_OPT_H

#include <vector>
#include <string>
#include <iostream>
#include <map>

using namespace std;

class OptStruct{
public:
	string opname;
	string opvalue;
	inline OptStruct(const string& _opname,const string& _opvalue):opname(_opname),opvalue(_opvalue){}
};

class EasyAdvGetOptOut{
public:
	string programName;
	vector<OptStruct> opts;
	vector<string> args;
	bool success;
	inline EasyAdvGetOptOut():success(false){}
	inline EasyAdvGetOptOut(string _programName,vector<OptStruct>& _opts,vector<string>& _args,bool _success):programName(_programName),opts(_opts),args(_args),success(_success){}
	inline void print(ostream& os=cout){
		for(vector<OptStruct>::iterator i=opts.begin();i!=opts.end();i++)
		{
			os<<"option "<<i->opname<<" : "<<i->opvalue<<endl;
		}
		
		for(unsigned int i=0;i<args.size();i++)
		{
			os<<"args "<<i<<" : "<<args[i]<<endl;
		}
		
	}
};

string argv2vectorOfString(vector<string>& vectorOfString,int argc,char* argv[]); //return programName
bool preprocessFileLoadableArgs(vector<string>& args,vector<string>& processedArgs);
bool getopt(vector<OptStruct>& out_opts,vector<string>& out_args,vector<string> &in_args, string options,vector<string>* long_options=NULL);
EasyAdvGetOptOut easyAdvGetOpt(int argc,char* argv[],string options,vector<string>* long_options=NULL);
void parseOptsIntoMap(vector<OptStruct>& opts,map<string,string>& optmap);
void parseOptsIntoMultiMap(vector<OptStruct>& opts,multimap<string,string>& optmap);
bool hasOpt(map<string,string>& optmap,const string& key);
bool hasOpt(multimap<string,string>& optmap,const string& key);
string getOptValue(map<string,string>& optmap,const string& key,const string& defaultValue="");
string getOptValue(multimap<string,string>& optmap,const string& key,const string& defaultValue="");
bool getOptValues(vector<string>& values, multimap<string,string>& optmap,const string& key);
bool checkRequiredOpts(map<string,string>& optmap,vector<string>& requiredOpts,const char*message="Error: option %s not specified.\n");
bool checkRequiredOpts(multimap<string,string>& optmap,vector<string>& requiredOpts,const char*message="Error: option %s not specified.\n");

#endif /*_ADV_GET_OPT_H*/