if [[ $CPPUTILCLASSES == "" ]]; then
	echo "\$CPPUTILCLASSES not specified"
	exit
fi

if [[ $CPPBIOCLASSES == "" ]]; then
	echo "\$CPPBIOCLASSES not specified"
	exit
fi

if [[ $SAMTOOLPATH == "" ]]; then
	echo "\$SAMTOOLPATH not specified"
	exit
fi

g++ -o geneRPKM -I$SAMTOOLPATH -I$CPPUTILCLASSES -I$CPPBIOCLASSES -L$SAMTOOLPATH -lbam -lz -lm geneRPKM_main.cpp AdvGetOptCpp/AdvGetOpt.cpp $SAMTOOLPATH/libbam.a 