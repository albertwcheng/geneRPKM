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

#include <iostream>
#include <fstream>
#include <set>
#include <Gff.h>
#include <BamUtil.h>
#include <libgen.h>
#include <SystemUtil.h>
#include "AdvGetOptCpp/AdvGetOpt.h"
using namespace std;
using namespace Gff;

#include <sam.h>

#define TABS "\n  "

void outArgsHelp(string argname,string help){
	cerr<<argname<<TABS<<help<<endl;
}



class OptionStruct {
public:
	string bamfile;
	string outfile;
	//bool bestQual;
	unsigned int maxHits;
	bool addNH;
	bool useNHFlag;
	string printStatFile;
};



void printUsage(string programName){
	
	cerr<<"Usage: "<<programName<<" [options] --in in.bam --out out.bam "<<endl;
	cerr<<"e.g.,"<<programName<<" --max-hits 3 --in in.bam --out out.bam"<<endl;
	cerr<<"preprocessor options:"<<endl;
	outArgsHelp("--@import-args filename","load arguments from tab delimited file");
	cerr<<"options:"<<endl;
	outArgsHelp("--max-hits hits","specify the maximum number of hits to retain the read");
	outArgsHelp("--add-NH","Add NH:i:<numHits> to the aux fields if not exists");
	outArgsHelp("--use-NH-flag","Use NH flag to filter which is way faster if it is present also if NH flag is consistent such that Left max hits == right max hits if both >0");
	outArgsHelp("--print-NH-stat-to","print NH stat to a file. either qname<tab>NH or qname<tab>firstAlgNH<tab>secondAlgNH");
	
}

class CountStruct{
	public:
		unsigned int firstAlignmentCount;
		unsigned int secondAlignmentCount;

		
		CountStruct(unsigned int _firstAlignmentCount=0,unsigned int _secondAlignmentCount=0):firstAlignmentCount(_firstAlignmentCount),secondAlignmentCount(_secondAlignmentCount){
			
		}
		
		inline CountStruct& operator += (const CountStruct& right){
			this->firstAlignmentCount+=right.firstAlignmentCount;
			this->secondAlignmentCount+=right.secondAlignmentCount;

		}
		
		/*inline CountStruct& withFirstFragQual(unsigned char qual){
			this->firstFragBestQual=qual;	
			return *this;
		}
		
		inline CountStruct& withSecondFragQual(unsigned char qual){
			this->secondFragBestQual=qual;
			
		}*/
		
		
};







int runGetUniqReads_twoPass(OptionStruct& opts){
	
	//CountPair firstAlignmentAdder(1,0);
	//CountPair secondAlignmentAdder(0,1);
	
	char qnameBuffer[256];
	
	samfile_t* bf=samopen(opts.bamfile.c_str(),"rb",0);



	if(!bf){
		cerr<<"bam file "<<opts.bamfile<<" cannot be open for counting"<<endl;
		return 1;
	}

	
	map<string,CountStruct> readHitsMap;
	
	
	
	bam1_t *bamInfo=bam_init1();

	//first pass
	unsigned int total;
	total=0;
	
	
	
	map<string,BamReader::BamAuxStruct> auxfields;
	
	while(samread(bf,bamInfo)>=0){
		total++;
		if(total%1000000==1){
			cerr<<"first pass: passing through read "<<total<<endl;
		}
		
		string qname=BamReader::getQName(bamInfo);
		//unsigned char qual=BamReader::getMappingQual(bamInfo);
		
				
		CountStruct Adder;
		
		if(BamReader::hasMultipleFragments(bamInfo))
		{
			if(BamReader::isFirstFragment(bamInfo))
			{
				Adder.firstAlignmentCount=1;
				Adder.secondAlignmentCount=0;			
			}
			else
			{
				Adder.firstAlignmentCount=0;
				Adder.secondAlignmentCount=1;	
			}
		}
		else
		{
			Adder.firstAlignmentCount=1;
			Adder.secondAlignmentCount=0;
		}
		
		map<string,CountStruct>::iterator i=readHitsMap.find(qname);
		if(i==readHitsMap.end()){
			readHitsMap.insert(map<string,CountStruct>::value_type(qname,Adder));
		}
		else
		{
			i->second+=Adder;
		}
	}

	
	samclose(bf);
	
	cerr<<"inReads\t"<<total<<endl;
	
	
	//check if left hits == right hits
	
	
	if(opts.printStatFile!=""){
		ofstream statOutFile(opts.printStatFile.c_str());
		
		unsigned int totalFirstAlignmentCount=0;
		unsigned int totalSecondAlignmentCount=0;
		
		unsigned int totalUniqFirstAlignmentCount=0;
		unsigned int totalUniqSecondAlignmentCount=0;
		
		
		statOutFile<<"QName\tFirstAlignmentCount\tSecondAlignmentCount"<<endl;
		
		
		for(map<string,CountStruct>::iterator i=readHitsMap.begin();i!=readHitsMap.end();i++)
		{
			statOutFile<<i->first<<"\t"<<i->second.firstAlignmentCount<<"\t"<<i->second.secondAlignmentCount<<endl;
			totalFirstAlignmentCount+=i->second.firstAlignmentCount;
			totalSecondAlignmentCount+=i->second.secondAlignmentCount;
			
			if(i->second.firstAlignmentCount>0){
				totalUniqFirstAlignmentCount++;
			}
			
			if(i->second.secondAlignmentCount>0){
				totalUniqSecondAlignmentCount++;	
			}
				
		}
		
		statOutFile<<"TotalAlignments"<<"\t"<<totalFirstAlignmentCount<<"\t"<<totalSecondAlignmentCount<<endl;
		
		statOutFile<<"TotalMapped"<<"\t"<<totalUniqFirstAlignmentCount<<"\t"<<totalUniqSecondAlignmentCount<<endl;
		
		statOutFile<<"TotalMappedUnion"<<"\t"<<readHitsMap.size()<<"\t"<<readHitsMap.size()<<endl;
		
		statOutFile.close();
		
		cerr<<"TotalAlignments"<<"\tread1="<<totalFirstAlignmentCount<<"\tread2="<<totalSecondAlignmentCount<<endl;
		cerr<<"TotalMapped"<<"\tread1="<<totalUniqFirstAlignmentCount<<"\tread2="<<totalUniqSecondAlignmentCount<<endl;
		cerr<<"TotalMappedUnion"<<"\t"<<readHitsMap.size()<<endl;
	}
	//second pass
	
	
	
	if(opts.outfile!="")
	{
		
		int outTotal=0;
		
		bf = samopen(opts.bamfile.c_str(),"rb",0);
		
		samfile_t* out= samopen(opts.outfile.c_str(), "wb", bf->header);
		
		if(!bf){
			cerr<<"bam file "<<opts.bamfile<<" cannot be open"<<endl;
			return 1;
		}
		total=0;
		
		
		
		
		while(samread(bf,bamInfo)>=0){
			total++;
			if(total%1000000==1){
				cerr<<"second pass: passing through read "<<total<<endl;
			}
			
			string qname=BamReader::getQName(bamInfo);
			//unsigned char qual=BamReader::getMappingQual(bamInfo);
			
			map<string,CountStruct>::iterator i=readHitsMap.find(qname);
			if(i==readHitsMap.end()){
				cerr<<"cannot found read "<<qname<<" in second pass?";
			}
			else{
				
	
				if(i->second.firstAlignmentCount<=opts.maxHits && i->second.secondAlignmentCount<=opts.maxHits){
					//only output these:
					//TODO: do we need the best?
					//can we directly write to a bam file?
					
					if(opts.addNH && !BamReader::hasAuxField(bamInfo,"NH")){
						int NH;
						if(!BamReader::hasMultipleFragments(bamInfo) || BamReader::isFirstFragment(bamInfo)){
							//first fragment
							NH=i->second.firstAlignmentCount;
						}
						else{
							NH=i->second.secondAlignmentCount;
						}
						
						BamReader::BamAuxStruct bas;
						bas.type=BAMAUX_INTVALUE;
						bas.intValue=NH;
						
						BamReader::appendAuxField(bamInfo,"NH",&bas);
						
					}
					
					bam_write1(out->x.bam,bamInfo);
					
					outTotal++;
				
				}
			}
		}
		
		samclose(out);
		samclose(bf);
		
		cerr<<"outReads\t"<<outTotal<<endl;
	}
		

	
	bam_destroy1(bamInfo);
	
	cerr<<"<Done>"<<endl;
	return 0;
}


int runGetUniqReads_useNHFlag(OptionStruct& opts){


	
	
	int total;
	
	samfile_t*  bf = samopen(opts.bamfile.c_str(),"rb",0);
	bam1_t *bamInfo=bam_init1();
	
	
	samfile_t* out=NULL;
	
	ofstream *statOutFile=NULL;
	set<string>* qnamesRecord=NULL; //remember those outputed ones
	
	if(opts.printStatFile!=""){
		statOutFile=new ofstream(opts.printStatFile.c_str());
		qnamesRecord=new set<string>;
	}
	
	
	if(opts.outfile!="")
		out = samopen(opts.outfile.c_str(), "wb", bf->header);
	
	if(!bf){
		cerr<<"bam file "<<opts.bamfile<<" cannot be open"<<endl;
		return 1;
	}
	total=0;
	
	
	
	
	while(samread(bf,bamInfo)>=0){
		total++;
		if(total%1000000==1){
			cerr<<"Passing through read "<<total<<endl;
		}
		
		
		//unsigned char qual=BamReader::getMappingQual(bamInfo);
		
		int numHits=BamReader::getNumHits(bamInfo,-1);
		string qname=BamReader::getQName(bamInfo);
		
		if(numHits==-1){
			//string qname=BamReader::getQName(bamInfo);
			cerr<<"Read at line "<<total<<" with name "<<qname<<" has no NH flag. abort";
			break;
		}
		
		if(statOutFile){
			
			if(qnamesRecord->find(qname)==qnamesRecord->end()){
				(*statOutFile)<<qname<<"\t"<<numHits<<endl;
				qnamesRecord->insert(qname);	
			}
		}
		
	
		if(numHits<=opts.maxHits){	
			
			if(out)
				bam_write1(out->x.bam,bamInfo);
				
		}
		
	}
	
	if(out)
		samclose(out);
	
	if(statOutFile){
		statOutFile->close();
		delete statOutFile;
	}
	
	if(qnamesRecord){
		delete qnamesRecord;	
	}
	
	samclose(bf);
	bam_destroy1(bamInfo);
	
	return 0;
}

int runGetUniqReads(OptionStruct& opts){
	//if(opts.useNHFlag){
	//	return runGetUniqReads_useNHFlag(opts);
	//}else{
		return runGetUniqReads_twoPass(opts);
	//}
}

int main(int argc,char*argv[])
{
	
	vector<string> long_options;
	
	long_options.push_back("max-hits=");
	long_options.push_back("in=");
	long_options.push_back("out=");
	long_options.push_back("add-NH");
	//long_options.push_back("use-NH-flag"); //cancel use NH flag
	long_options.push_back("print-NH-stat-to=");
	
	//long_options.push_bacl("out-best-qual");
	
	
	OptionStruct opts;
	map<string,string> optmap;
	
	
	
	EasyAdvGetOptOut argsFinal=easyAdvGetOpt(argc,argv,"",&long_options);
	
	
	if(argsFinal.success){
		argsFinal.print(cerr);
		
		
		parseOptsIntoMap(argsFinal.opts,optmap);
		
	}
	else{
		printUsage(argsFinal.programName);
		return 1;
	}
	
	
	opts.maxHits=atoi(getOptValue(optmap,"--max-hits","1").c_str());
	if(!hasOpt(optmap,"--in")){
		cerr<<"input bam file not specified. abort"<<endl;
		printUsage(argsFinal.programName);
		return 1;	
	}
	
	/*if(!hasOpt(optmap,"--out")){
		cerr<<"output bam file not specified. abort"<<endl;
		printUsage(argsFinal.programName);
		return 1;
	}*/
	
	opts.bamfile=getOptValue(optmap,"--in");
	opts.outfile=getOptValue(optmap,"--out","");
	opts.addNH=hasOpt(optmap,"--add-NH");
	//opts.useNHFlag=hasOpt(optmap,"--use-NH-flag");
	opts.printStatFile=getOptValue(optmap,"--print-NH-stat-to","");
	//opts.bestQual=hasOpt(optmap,"--out-best-qual");
		
	
	
	return runGetUniqReads(opts);
	

}