#include <stdio.h>
#include <stdlib.h>
#include <Correct.h>
#include <htab.h>
#include <iostream>
#include "CommandLines.h"
#include "ketopt.h"
int aaamain(int argc, char *argv[])
{

    
    int i, ret;
    yak_reset_realtime();
    int c;
    init_opt1(&asm_opt);
    std::cout<<"chenggong"<<std::endl;
    // if (!CommandLine_process(argc, argv, &asm_opt)) return 0;
    
    // //std::cout<<"success"<<std::endl;
    // //ret = ha_assemble();
    // std::cout<<"success"<<std::endl;
    // ret = ha_assemble();
    // destory_opt(&asm_opt);
	// fprintf(stderr, "[M::%s] Version: %s\n", __func__, HA_VERSION);
	// fprintf(stderr, "[M::%s] CMD:", __func__);


	// for (i = 0; i < argc; ++i)
	// 	fprintf(stderr, " %s", argv[i]);
	// fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, yak_realtime(), yak_cputime(), yak_peakrss_in_gb());
    // return ret;
    // return 0;
}