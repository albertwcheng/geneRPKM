/***************************************************************************
 Copyright 2010 Wu Albert Cheng <albertwcheng@gmail.com>
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


#include "AdvGetOpt.h"
#include <string.h>
#include <map>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

#define FILE_BUFFER_LENGTH 10240
#define FILE_ARGS_SPLIT "\t"

string preprocessFileLoadableArgs_inner_formkeystring(char* key,string extraPrefix="")
{
	
string prefix="";
	
if(key[0]!='-')
{
	if (strlen(key)==1)
	{
		prefix="-";
	}
	else
	{
		prefix="--";
	}
}
	
char* ckeyptr=key;
while(*ckeyptr!='\0'){
	if(*ckeyptr=='_')
		*ckeyptr='-';
	else if(*ckeyptr==' ')
		*ckeyptr='-';
	
	ckeyptr++;
}
	
	return prefix+extraPrefix+key;
}

string argv2vectorOfString(vector<string>& vectorOfString,int argc,char* argv[])
{
	for(int i=1;i<argc;i++)
		vectorOfString.push_back(string(argv[i]));		
	
	return string(argv[0]);
}

bool preprocessFileLoadableArgs(vector<string>& args,vector<string>& processedArgs)

{
	//cerr<<"size="<<args.size()<<endl;
	char buffer[FILE_BUFFER_LENGTH];
	
	int i=0;
	processedArgs.clear();
	int argLength=args.size();
	while(i<argLength)
	{
		string& avalue=args[i];
		if (avalue=="--@import-args")
		{
			if(i<argLength-1)
			{
				string& filename=args[i+1];
				ifstream fil(filename.c_str());
				if(!fil.good()){
					//file not good
					cerr<<"error: arg file "<<filename<<" not good"<<endl;
					return 1;
				}
				
				while(fil.good())
				{	
					strcpy(buffer,"");
					fil.getline(buffer,FILE_BUFFER_LENGTH);
					
					if(strlen(buffer)<1)
						break; //no more to load
					
					if(buffer[0]=='#'){
						//ignore this commented line
						continue;
					}
					
					char* pch;
					pch=strtok(buffer,FILE_ARGS_SPLIT);
					int j=0;
					while(pch!=NULL)
					{
						if(j==0){
							string argname=preprocessFileLoadableArgs_inner_formkeystring(pch);
							processedArgs.push_back(argname);
						}else {
							processedArgs.push_back(string(pch));
						}

						j++;
						pch=strtok(NULL,FILE_ARGS_SPLIT);
					}
					
					
				}
				
				fil.close();
			}else {
				cerr<<"error: --@import-args not followed by a filename"<<endl;
				return false;
			}
			
			i++;

		}
		else if(avalue=="--@import-cfg")
		{
			//not supported yet
			cerr<<"--@import-cfg not supported yet"<<endl;
			i++;
			return false;
		}else{
			processedArgs.push_back(avalue);
		}
		
		i++;
	}
				
	return true;
	
}



bool getopt(vector<OptStruct>& out_opts,vector<string>& out_args,vector<string> &in_args, string options,vector<string>* long_options)
{
	//build arg list	
	map<string,bool> optionlist; //optionlist[optname]=<whether it consumes an argvalue>
	//process short option first "ab:c" -> "-a",F "-b",T "c",F
	
	unsigned int optionLength=options.length();
 	for(unsigned int i=0;i<optionLength;i++){
		string optname=string("-")+options.substr(i,1);
		bool requireArgValue=false;
		if(i<optionLength-1 && options[i+1]==':'){
			requireArgValue=true;
			i++;
		}
		optionlist.insert(map<string,bool>::value_type(optname,requireArgValue));
		
	}
	
	//now process long options "a-1=", "a-2" => "--a-1",T "--a-2",F
	if(long_options)
	{
		for(vector<string>::iterator i=long_options->begin();i!=long_options->end();i++)
		{
			string optname=string("--")+(*i);
			int optnameLength=optname.length();
			bool requireArgValue=false;
			if(optname[optnameLength-1]=='='){
				requireArgValue=true;
				optname.erase(optnameLength-1);
			}
			
			optionlist.insert(map<string,bool>::value_type(optname,requireArgValue));
		}
	}
	
	int i=0;
	int argLength=in_args.size();
	//cerr<<"size="<<argLength<<endl;
	//now process in_args
	while(i<argLength){
		string& avalue=in_args[i];
		if(avalue[0]=='-')
		{
			map<string,bool>::iterator findOptI=optionlist.find(avalue);
			if (findOptI==optionlist.end()) {
				cerr<<"Error: Option "<<avalue<<" not found"<<endl;
				return false;
			}
			
			string& optname=avalue;
			string optvalue="";
			
			bool requireArgValue=findOptI->second;
			
			if(requireArgValue)
			{
				if(i>=argLength-1){
					cerr<<"Error: Option "<<optname<<" not followed by a required optvalue"<<endl;
					return false;
				}
				optvalue=in_args[i+1];
				i++;
			}
			
			out_opts.push_back(OptStruct(optname,optvalue));
		}else{
			out_args.push_back(avalue);
		}
		i++;
	}
		   
	return true;
	
}

EasyAdvGetOptOut easyAdvGetOpt(int argc,char* argv[],string options,vector<string>* long_options){
	vector<string> preprocessed_in_args;
	vector<string> processed_in_args;
	EasyAdvGetOptOut output;
	output.programName=argv2vectorOfString(preprocessed_in_args,argc,argv);
	output.success=preprocessFileLoadableArgs(preprocessed_in_args,processed_in_args);
	if(!output.success)
		return output;
	output.success=getopt(output.opts,output.args,processed_in_args, options,long_options);
	
	return output;
}

void parseOptsIntoMap(vector<OptStruct>& opts,map<string,string>& optmap){
	for(vector<OptStruct>::reverse_iterator i=opts.rbegin();i!=opts.rend();i++) //the latter values have more precedence!
		optmap.insert(map<string,string>::value_type(i->opname,i->opvalue));
}
void parseOptsIntoMultiMap(vector<OptStruct>& opts,multimap<string,string>& optmap){
	for(vector<OptStruct>::iterator i=opts.begin();i!=opts.end();i++)
		optmap.insert(multimap<string,string>::value_type(i->opname,i->opvalue));
}
bool hasOpt(map<string,string>& optmap,const string& key){
	return !(optmap.find(key)==optmap.end());
}
bool hasOpt(multimap<string,string>&optmap,const string& key){
	return !(optmap.find(key)==optmap.end());	
}
string getOptValue(map<string,string>& optmap,const string& key,const string& defaultValue){
	map<string,string>::iterator findI=optmap.find(key);
	if(findI==optmap.end())
		return defaultValue;
	
	return findI->second;
}

string getOptValue(multimap<string,string>& optmap,const string& key,const string& defaultValue){
	
	pair<multimap<string,string>::iterator,multimap<string,string>::iterator> findI=optmap.equal_range(key);
	
	if(findI.first==optmap.end())
		return defaultValue;
	
	//for(multimap<string,string>::iterator i=findI.first;i!=findI.second;i++)
	//	values.push_back(i->second);
	return findI.first->second;
}


bool getOptValues(vector<string>& values, multimap<string,string>& optmap,const string& key){
	values.clear();
	pair<multimap<string,string>::iterator,multimap<string,string>::iterator> findI=optmap.equal_range(key);
	if(findI.first==optmap.end())
		return false;
	
	for(multimap<string,string>::iterator i=findI.first;i!=findI.second;i++)
		values.push_back(i->second);
	
	return true;
}

bool checkRequiredOpts(map<string,string>& optmap,vector<string>& requiredOpts,const char*message){
	bool ok=true;
	for(vector<string>::iterator i=requiredOpts.begin();i!=requiredOpts.end();i++)
	{
		string optname=*i;
		if (optname[0]!='-') {
			if(optname.length()==1){
				optname="-"+optname;
			}
			else{
				optname="--"+optname;
			}
		}
		if(optname[optname.length()-1]=='=')
			optname.erase(optname.length()-1);
		
		if(!hasOpt(optmap,optname)){
			ok=false;
			fprintf(stderr,message,optname.c_str());
		}
	}
	
	return ok;
}
bool checkRequiredOpts(multimap<string,string>& optmap,vector<string>& requiredOpts,const char*message){
	bool ok=true;
	for(vector<string>::iterator i=requiredOpts.begin();i!=requiredOpts.end();i++)
	{
		string optname=*i;
		if (optname[0]!='-') {
			if(optname.length()==1){
				optname="-"+optname;
			}
			else{
				optname="--"+optname;
			}
		}
		if(optname[optname.length()-1]=='=')
			optname.erase(optname.length()-1);
		
		if(!hasOpt(optmap,optname)){
			ok=false;
			fprintf(stderr,message,optname.c_str());
		}
	}
	
	return ok;
}