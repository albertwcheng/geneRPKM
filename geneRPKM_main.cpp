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
#include <math.h>
using namespace std;
using namespace Gff;

#include <sam.h>

#define TABS "\n  "

void outArgsHelp(string argname,string help){
	cerr<<argname<<TABS<<help<<endl;
}

#define EXPRESSIONMODE_RPKM     1
#define EXPRESSIONMODE_RPKM_DIVHITS 2
#define EXPRESSIONMODE_FPKM  3
#define EXPRESSIONMODE_FPKM_DIVHITS 4

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
	bool forceFlexMaxBasepairPolicy;
	
	ofstream* regionBedOutStream;
	string regionBedOut;
	
	string noBlockBedOut;
	ofstream *noBlockBedOutStream;
	
	
	int expressionMode; //EXPRESSIONMODE_*
	int maxHits;
	string itemRgb;
	string fillNA;
	string prefixDataLabel;
	
	OptionStruct():totalNumOfReads(0),constitutiveThresholdFrac(0.0),constitutiveThresholdNum(0),flexmaxThresholding(false),flexmaxThreshold(0),maxHits(0),regionBedOutStream(NULL),noBlockBedOutStream(NULL),itemRgb("0,0,0"){}
	~OptionStruct(){
		if(regionBedOutStream){
			regionBedOutStream->close();
			delete regionBedOutStream;
		}
		
		if(noBlockBedOutStream){
			noBlockBedOutStream->close();
			delete noBlockBedOutStream;
		}
	}
};


/* output to stdout
 
1) GeneName
2) Chrom
3) Gene Start (1-based)
4) Gene End (1-based)
5) Strand
6) Length of probed region
7) Read Counts
8) RPKM or RPKM
9) log2(RPKM) or log2(FPKM)
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
	outArgsHelp("--bamfile bamfile","specify the bam file. Repeat this option for multiple bam files");
	outArgsHelp("--bedfile befile","specify the bed file for the gene. Can repeat this option for multiple gene annotations");
	outArgsHelp("--total-num-reads numReads","Specify total number of mapped reads. Instead of counting from the bam file");
	outArgsHelp("--constitutive-threshold-frac fraction","Specify the min fraction [0.0,1.0] of transcripts covering a region to be used in counting. Default: 1.-");
	outArgsHelp("--constitutive-threshold-num num","Specify the min number of transcripts covering a region to be used in counting. Default: 1");
	outArgsHelp("--flexmax-thresholding bp","use a flexible max threshold such that all genes can be counted with at least bp number of basepairs");
	outArgsHelp("--force-flexmax-bp-policy","enforce flexmax bp policy, if going down to 1 blockFreq and no bp width of region available, then fail the gene expression detection for that gene");
	outArgsHelp("--fpkm-divhits","[default] count reads and expression as FPKM after normalizing read numbers to number of hits (using NH flag) for each read");
	outArgsHelp("--fpkm","count reads and expression as FPKM");
	outArgsHelp("--rpkm-divhits","count reads and expression as RPKM after normalizing read numbers to the number of hits (using NH flag) for each read");
	outArgsHelp("--rpkm","count reads and expression as RPKM");
	outArgsHelp("--max-hits hits","discard reads with more than a certain number of hits");
	outArgsHelp("--label-prefix str","prefix label of RPKM/FPKM, etc by str");
	outArgsHelp("--fill-NA-with str","fill NA data with str");
	cerr<<"\tNote that both thresholds are applied unless --flexmax-thresholding is specified"<<endl;
	outArgsHelp("--region-bed-out bedfile","output a regions used to calculate to a bedfile");
	outArgsHelp("--region-bed-itemRgb rgb","set the itemRgb output for the bed out");
	outArgsHelp("--no-block-bed-out bedfile","output a bed file consisting of genes with no available blocks for gene expression estimation according to the current settings");
	//outArgsHelp("--use-coding-region-only","whether to use only coding region (for genes that have coding regions");
	
}

int runGeneRPKM(OptionStruct& opts){
	
	//total number of reads not specified. count from bam files.
	if(opts.totalNumOfReads==0){
		double totalNumOfReadsT=0.0;
		for(vector<string>::iterator i=opts.bamfilenames.begin();i!=opts.bamfilenames.end();i++){
			switch(opts.expressionMode){
				case EXPRESSIONMODE_RPKM:
					totalNumOfReadsT+=BamReader::countTotalNumOfReadsInBam(*i);
					break;
				case EXPRESSIONMODE_RPKM_DIVHITS:
					totalNumOfReadsT+=BamReader::countTotalNumOfReadsInBamDivHits(*i);
					break;
				case EXPRESSIONMODE_FPKM:
					totalNumOfReadsT+=BamReader::countTotalNumOfFragmentsInBam(*i);
					break;
				case EXPRESSIONMODE_FPKM_DIVHITS:
					totalNumOfReadsT+=BamReader::countTotalNumOfFragmentsInBamDivHits(*i);
					break;
				default:
					break;
					
			}
		}
		
		opts.totalNumOfReads=ceil(totalNumOfReadsT);
		
		cerr<<"total number of reads is "<<opts.totalNumOfReads<<endl;
	}
	
	//return 0;
	
	double countD=0.0;
	int countI=0;
	
	BamReader::AdvFragmentSetCounter afsc;
	
	/*
	1) GeneName
	2) Chrom
	3) Gene Start (1-based)
	4) Gene End (1-based)
	5) Strand
	6) Length of probed region
	7) Read Counts
	8) RPKM or RPKM
	 9) log2(RPKM) or log2(FPKM)
	 11) MinConsUsedFrac
	 12) MinConsUsedNum
	 13) TotalNumberOfReads
	 */
	//write header
	cout<<"GeneName"<<"\t";
	cout<<"Chrom"<<"\t";
	cout<<"GeneStart"<<"\t";
	cout<<"GeneEnd"<<"\t";
	cout<<"Strand"<<"\t";
	cout<<"LengthProbed"<<"\t";
	cout<<opts.prefixDataLabel<<"ReadCounts"<<"\t";
	switch (opts.expressionMode) {
		case EXPRESSIONMODE_FPKM:case EXPRESSIONMODE_FPKM_DIVHITS:
			cout<<opts.prefixDataLabel<<"FPKM"<<"\t";
			cout<<opts.prefixDataLabel<<"log2(FPKM)"<<"\t";
			//cout<<opts.prefixDataLabel<<"log2(1+FPKM)"<<"\t";
			
			break;
		case EXPRESSIONMODE_RPKM:case EXPRESSIONMODE_RPKM_DIVHITS:
			cout<<opts.prefixDataLabel<<"RPKM"<<"\t";
			cout<<opts.prefixDataLabel<<"log2(RPKM)"<<"\t";
			//cout<<opts.prefixDataLabel<<"log2(1+RPKM)"<<"\t";
			break;
		default:
			break;
	}
	cout<<"MinConsUsedFrac"<<"\t";
	cout<<"MinConsUsedNum"<<"\t";
	cout<<"TotalNumberOfReads"<<endl;
	
	//now go to each genes and count
	for(vector<Annotation*>::iterator annoI=opts.annotations.begin();annoI!=opts.annotations.end();annoI++){
		Annotation* pannot=*annoI;
		for(Annotation::NameGeneMapI nameGeneI=pannot->name_genes_begin();nameGeneI!=pannot->name_genes_end();nameGeneI++){
			//cerr<<"processing gene "<<nameGeneI->first<<endl;
			SmartPtr<Gene> gene=nameGeneI->second;
			
			double countD=0.0;
			

			pair<int,float> minUsed;
			
			
			set<SortPair<Coord,Coord> > blocks;
			///TODO
			if(opts.flexmaxThresholding){
				///add here
				minUsed=gene->getFlexMaxConstitutiveBlocks(blocks,opts.flexmaxThreshold,opts.forceFlexMaxBasepairPolicy);
			}else{
				minUsed=gene->getConstitutiveBlocks(blocks,opts.constitutiveThresholdFrac,opts.constitutiveThresholdNum);
			}
			
			//if we have regionBedOutStream
			if(blocks.size()>0){
				if(opts.regionBedOutStream){
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
					(*opts.regionBedOutStream)<<gene->chrom()<<"\t";
					int regionChromStart0=blocks.begin()->k1;
					int regionChromEnd1=blocks.rbegin()->k2;
					(*opts.regionBedOutStream)<<regionChromStart0<<"\t";
					(*opts.regionBedOutStream)<<regionChromEnd1<<"\t";
					(*opts.regionBedOutStream)<<gene->name()<<"\t";
					(*opts.regionBedOutStream)<<"0"<<"\t";
					(*opts.regionBedOutStream)<<gene->strand<<"\t";
					(*opts.regionBedOutStream)<<regionChromStart0<<"\t";
					(*opts.regionBedOutStream)<<regionChromEnd1<<"\t";
					(*opts.regionBedOutStream)<<opts.itemRgb<<"\t";
					(*opts.regionBedOutStream)<<blocks.size()<<"\t";
				
					
					set<SortPair<Coord,Coord> >::iterator i=blocks.begin();
					string blockSizes=StringUtil::str(i->k2-i->k1);
					string blockStarts=StringUtil::str(i->k1-regionChromStart0);
					i++;
					for(;i!=blocks.end();i++){
						blockSizes+=","+StringUtil::str(i->k2-i->k1);
						blockStarts+=","+StringUtil::str(i->k1-regionChromStart0);
					}
					(*opts.regionBedOutStream)<<blockSizes<<"\t";
					(*opts.regionBedOutStream)<<blockStarts<<endl;
					
					
				}
			}else{
				if(opts.noBlockBedOutStream){
					(*opts.noBlockBedOutStream)<<gene->chrom()<<"\t";
					(*opts.noBlockBedOutStream)<<gene->start0()<<"\t";
					(*opts.noBlockBedOutStream)<<gene->end1()<<"\t";
					(*opts.noBlockBedOutStream)<<gene->name()<<"\t";
					(*opts.noBlockBedOutStream)<<"0"<<"\t";
					(*opts.noBlockBedOutStream)<<gene->strand<<endl;
				}
			}
			
			//Now we have the blocks for expression calculation
			//go to each block, get number of reads or fragments.		
			
			int lengthProbed=0;
			
			bool hasChromInAnyBams=false;
			
			for(vector<BamReader*>::iterator bfi=opts.bamfiles.begin();bfi!=opts.bamfiles.end();bfi++){
				
				BamReader* curBam=*bfi;
				
				if(!curBam->hasChromInBam(gene->chrom())){
					continue;
				}
				
				hasChromInAnyBams=true;
				
				afsc.resetCount();
				lengthProbed=0;
				for(set<SortPair<Coord,Coord> >::iterator i=blocks.begin();i!=blocks.end();i++){
					
					int blockStart0=i->k1;
					int blockEnd1=i->k2;
					lengthProbed+=blockEnd1-blockStart0;
					
					switch (opts.expressionMode) {
						case EXPRESSIONMODE_FPKM:
							curBam->fetchFragmentCountOverlappingRegion(afsc,gene->chrom(),blockStart0,blockEnd1,true,opts.maxHits);
							break;
						case EXPRESSIONMODE_FPKM_DIVHITS:
							curBam->fetchFragmentCountOverlappingRegionDivideByNumHits(afsc,gene->chrom(),blockStart0,blockEnd1,true,opts.maxHits);
							break;
						case EXPRESSIONMODE_RPKM:
							countD+=curBam->fetchCountOverlappingRegion(gene->chrom(),blockStart0,blockEnd1,true,opts.maxHits);
							break;
						case EXPRESSIONMODE_RPKM_DIVHITS:
							countD+=curBam->fetchCountOverlappingRegionDivideByNumHits(gene->chrom(),blockStart0,blockEnd1,true,opts.maxHits);
							break;
						default:
							break;
					}
				}
				
				
				//transmit read info from afsc to countI or countD
				if(opts.expressionMode==EXPRESSIONMODE_FPKM || opts.expressionMode==EXPRESSIONMODE_FPKM_DIVHITS){
					countD+=afsc;
				}
				
				
				
				
				
			}
			
			//end of gene: output!
			/*
			1) GeneName
			2) Chrom
			3) Gene Start (1-based)
			4) Gene End (1-based)
			5) Strand
			6) Length of probed region
			7) Read Counts
			8) RPKM or RPKM
			9) log2(RPKM) or log2(FPKM)
			10) MinConsUsedFrac
			11) MinConsUsedNum
			12) TotalNumberOfReads			
			*/
			
			cout<<gene->name()<<"\t";
			cout<<gene->chrom()<<"\t";
			cout<<(gene->start0()+1)<<"\t";
			cout<<gene->end1()<<"\t";
			cout<<gene->strand<<"\t";
			cout<<lengthProbed<<"\t";
			if(lengthProbed==0 || !hasChromInAnyBams){
				cout<<opts.fillNA<<"\t";
				cout<<opts.fillNA<<"\t";
				cout<<opts.fillNA<<"\t";
				//cout<<opts.fillNA<<"\t";
			}else{
				cout<<countD<<"\t";
				double RPKM=countD/(float(opts.totalNumOfReads)/1e6)/(float(lengthProbed)/1e3);
				
				cout<<RPKM<<"\t";
				if(RPKM==0.0){
					cout<<opts.fillNA<<"\t";
				}else {
					cout<<(log(RPKM)/log(2))<<"\t";
				}

				
				//cout<<(log(RPKM+1.0)/log(2))<<"\t";
			}
			cout<<minUsed.second<<"\t";
			cout<<minUsed.first<<"\t";
			cout<<opts.totalNumOfReads<<endl;
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
	long_options.push_back("region-bed-itemRgb=");
	long_options.push_back("max-hits=");
	long_options.push_back("force-flexmax-bp-policy");
	long_options.push_back("rpkm");
	long_options.push_back("fpkm");
	long_options.push_back("rpkm-divhits");
	long_options.push_back("fpkm-divhits");
	long_options.push_back("label-prefix=");
	long_options.push_back("fill-NA-with=");
	long_options.push_back("no-block-bed-out=");
	
	
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
	opts.constitutiveThresholdNum=atoi(getOptValue(optmap,"--constitutive-threshold-num","1").c_str());
	
	if(hasOpt(optmap,"--flexmax-thresholding"))
	{	
		opts.flexmaxThresholding=true;
		opts.flexmaxThreshold=atoi(getOptValue(optmap,"--flexmax-thresholding").c_str());
	}
	
	//cerr<<opts.flexmaxThresholding<<endl;
	//cerr<<opts.constitutiveThresholdFrac<<endl;
	//cerr<<opts.constitutiveThresholdNum<<endl;
	
	opts.regionBedOut=getOptValue(optmap,"--region-bed-out","");
	
	opts.maxHits=atoi(getOptValue(optmap,"--max-hits","0").c_str());
	
	opts.expressionMode=EXPRESSIONMODE_FPKM_DIVHITS; //default
	
	opts.regionBedOut=getOptValue(optmap,"--region-bed-out");
	opts.noBlockBedOut=getOptValue(optmap,"--no-block-bed-out");
	
	opts.itemRgb=getOptValue(optmap,"--region-bed-itemRgb","0,0,0");
	
	opts.prefixDataLabel=getOptValue(optmap,"--label-prefix","");
	opts.fillNA=getOptValue(optmap,"--fill-NA-with","NA");
	
	if(opts.regionBedOut.length()>0){
		opts.regionBedOutStream=new ofstream(opts.regionBedOut.c_str());
	}
	
	if(opts.noBlockBedOut.length()>0){
		opts.noBlockBedOutStream=new ofstream(opts.noBlockBedOut.c_str());
	}
	
	if(hasOpt(optmap,"--rpkm")){
		opts.expressionMode=EXPRESSIONMODE_RPKM;
	}else if(hasOpt(optmap,"--rpkm-divhits")){
		opts.expressionMode=EXPRESSIONMODE_RPKM_DIVHITS;
	}else if(hasOpt(optmap,"--fpkm")){
		opts.expressionMode=EXPRESSIONMODE_FPKM;
	}else if(hasOpt(optmap,"--fpkm-divhits")){
		//this is default anyway
		opts.expressionMode=EXPRESSIONMODE_FPKM_DIVHITS;
	}
	
	opts.forceFlexMaxBasepairPolicy=hasOpt(optmap,"--force-flexmax-bp-policy");
	
	
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