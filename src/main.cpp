#ifndef MAIN_CPP
#define MAIN_CPP
#include "utils.hpp"
#include <dechat_1.hpp>
#include <cstdlib>
#include "CommandLines.h"
#include "CommandLines1.h"
#include "Correct.h"
#include "ketopt.h"
#include "correct_round1.h"
chat_opt_t chat_opt;
int main(int argc, char **argv)
{

	
	init_opt(&chat_opt);
	//PRINT_LINE_FUNC();
	if (!Dechat_command(argc, argv, &chat_opt, &asm_opt))
		return 0;
	//PRINT_LINE_FUNC();
	// //.......................................Round 1.......................................
	std::cout << "asm_opt->output_file_name:" << asm_opt.output_file_name << std::endl;
	std::cout << "chat_opt.fast:" << chat_opt.fast << std::endl;

	correct_round1(&chat_opt);
	

	correct_round2(&chat_opt, &asm_opt);
	std::cout << "chat_opt.max_ov_diff_final:" << chat_opt.max_ov_diff_final << std::endl;
	std::cout << "chat_opt.abundance_min:" << chat_opt.abundance_min << std::endl;
	std::cout << "chat_opt.second_number_of_round:" << chat_opt.second_number_of_round << std::endl;
	destory_opt(&asm_opt);
	fprintf(stderr, "[M::%s] Version: %s\n", __func__, HA_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (int i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, yak_realtime(), yak_cputime(), yak_peakrss_in_gb());

	return EXIT_SUCCESS;
}

#endif // MAIN_CPP