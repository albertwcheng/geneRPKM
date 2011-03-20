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
	vector<BamReader*> bamfiles;
	vector<Annotation*> annotations;
	vector<string> bamfilenames;
	vector<string> bedfilenames;
	int totalNumOfReads;
	double constitutiveThresholdFrac;
	int constitutiveThresholdNum;
	bool flexmaxThresholding;
	int flexmaxThreshold; //numbp
	string regionBedOut;
	OptionStruct():totalNumOfReads(0),constitutiveThresholdFrac(0.0),constitutiveThresholdNum(0),flexmaxThresholding(false),flexmaxThreshold(0){}
};


/* output to stdout
 
1) GeneName
2) Chrom
3) Gene Start (1-based)
4) Gene End (1-based)
5) Strand
6) Length of probed region
7) Read Counts
8) RPKM
9) log2(1+RPKM)
10) MinConsUsedFrac
11) MinConsUsedNum
12) TotalNumberOfReads
 
//not implemented
11) TotalNumberOfReadsExonic
12) TotalNumberOfReadsIntronic
13) TotalNumberOfReadsIntegenic

 
 */


/* region bed out
 
 1) chrom
 2) chromStart (0-based)
 3) chromEnd (1-based)
 4) geneName
 5) score = 0
 6) strand
 7) thickStart=chromStart
 8) thickEnd=chromEnd
 9) itemRgb=0,0,0
 10) blockCount
 11) blockSizes
 12) blockStarts
 
 
 */

void printUsage(string programName){
	
	cerr<<"Usage: "<<programName<<" [options]"<<endl;
	cerr<<"preprocessor options:"<<endl;
	outArgsHelp("--@import-args filename","load arguments from tab delimited file");
	cerr<<"options:"<<endl;
	outArgsHelp("--bamfile bamfile","specify the bam file");
	outArgsHelp("--bedfile befile","specify the bed file for the gene. Can repeat this option for multiple gene annotations");
	outArgsHelp("--total-num-reads numReads","Specify total number of mapped reads. Instead of counting from the bam file");
	outArgsHelp("--constitutive-threshold-frac fraction","Specify the min fraction [0.0,1.0] of transcripts covering a region to be used in counting. Default: 1.-");
	outArgsHelp("--constitutive-threshold-num num","Specify the min number of transcripts covering a region tobe used in counting. Default: 0");
	outArgsHelp("--flexmax-thresholding bp","use a flexible max threshold such that all genes can be counted with at least bp number of basepairs");
	cerr<<"\tNote that both thresholds are applied unless --flexmax-thresholding is specified"<<endl;
	outArgsHelp("--region-bed-out bedfile","output a regions used to calculate to a bedfile");
	//outArgsHelp("--use-coding-region-only","whether to use only coding region (for genes that have coding regions");
	
}

int runGeneRPKM(OptionStruct& opts){
	
	//total number of reads not specified. count from bam files.
	if(opts.totalNumOfReads==0){
		for(vector<string>::iterator i=opts.bamfilenames.begin();i!=opts.bamfilenames.end();i++){
			opts.totalNumOfReads+=BamReader::countTotalNumOfReadsInBam(*i);
		}
		
		cerr<<"total number of reads is "<<opts.totalNumOfReads<<endl;
	}
	
	
	
	//now go to each genes and count
	for(vector<Annotation*>::iterator annoI=opts.annotations.begin();annoI!=opts.annotations.end();annoI++){
		Annotation* pannot=*annoI;
		for(Annotation::NameGeneMapI nameGeneI=pannot->name_genes_begin();nameGeneI!=pannot->name_genes_end();nameGeneI++){
			cerr<<"processing gene "<<nameGeneI->first<<endl;
		}
	}
		
	return 1;
}

int main(int argc,char*argv[])
{
	
	vector<string> long_options;
	//required:
	long_options.push_back("bamfile=");
	long_options.push_back("bedfile=");
	long_options.push_back("total-num-reads=");
	long_options.push_back("constitutive-threshold-frac=");
	long_options.push_back("constitutive-threshold-num=");
	long_options.push_back("flexmax-thresholding=");
	long_options.push_back("region-bed-out=");

	
	OptionStruct opts;
	multimap<string,string> optmap;
	
	

	EasyAdvGetOptOut argsFinal=easyAdvGetOpt(argc,argv,"",&long_options);
	
	
	if(argsFinal.success){
		argsFinal.print(cerr);
		
		
		parseOptsIntoMultiMap(argsFinal.opts,optmap);
		
	}
	else{
		printUsage(argsFinal.programName);
		return 1;
	}
	
	
	vector<string> tmpNum;
	
	getOptValues(tmpNum,optmap,"--total-num-reads");
	
	for(vector<string>::iterator i=tmpNum.begin();i!=tmpNum.end();i++){
		opts.totalNumOfReads+=atoi(i->c_str());
	}
	
	
	getOptValues(opts.bamfilenames,optmap,"--bamfile");
	getOptValues(opts.bedfilenames,optmap,"--bedfile");
	
	opts.constitutiveThresholdFrac=atof(getOptValue(optmap,"--constitutive-threshold-frac","1.0").c_str());
	opts.constitutiveThresholdNum=atoi(getOptValue(optmap,"--constitutive-threshold-num","0").c_str());
	
	if(hasOpt(optmap,"--flexmax-thresholding"))
	{	
		opts.flexmaxThresholding=true;
		opts.flexmaxThreshold=atoi(getOptValue(optmap,"--flexmax-thresholding").c_str());
	}
	
	
	opts.regionBedOut=getOptValue(optmap,"--region-bed-out","");
	
	
	
	if(opts.bamfilenames.size()==0){
		cerr<<"no bam file specified"<<endl;
		printUsage(argsFinal.programName);
		return 1;
	}
	
	if(opts.bedfilenames.size()==0){
		cerr<<"no bed file specified"<<endl;
		printUsage(argsFinal.programName);
		return 1;
	}
	
	for(vector<string>::iterator i=opts.bamfilenames.begin();i!=opts.bamfilenames.end();i++){
		opts.bamfiles.push_back(new BamReader(*i));
	}
	
	for(vector<string>::iterator i=opts.bedfilenames.begin();i!=opts.bedfilenames.end();i++){
		Annotation* annot=new Annotation;
		annot->readBedFile(*i);
		opts.annotations.push_back(annot);
	}
	
	int success_status=runGeneRPKM(opts);
	
	//now clean up
	for(vector<BamReader*>::iterator i=opts.bamfiles.begin();i!=opts.bamfiles.end();i++)
	{
		(*i)->close();
		delete *i;
	}
	
	for(vector<Annotation*>::iterator i=opts.annotations.begin();i!=opts.annotations.end();i++){
		delete *i;
	}
	
	return success_status;

}