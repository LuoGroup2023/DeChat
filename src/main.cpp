#ifndef MAIN_CPP
#define MAIN_CPP
#include "utils.hpp"
#include <dechat_1.hpp>
#include <cstdlib>
#include <Aligner.h>
#include "CommandLines.h"
#include "CommandLines1.h"
#include "Correct.h"
#include "ketopt.h"
#include "correct_round1.h"
chat_opt_t chat_opt;
int main(int argc, char **argv)
{

	// try
	// {
	// 	dechat_1().run(argc, argv);
	// }
	// catch (Exception &e)
	// {
	// 	std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
	// 	return EXIT_FAILURE;
	// }
	// std::cout << "End dechat" << std::endl;
	//.......................................DBG default 参数.......................................
	// params.minimizerSeedDensity = 5;
	// params.minimizerLength = 19;
	// 			params.minimizerWindowSize = 30;
	// 			params.maxClusterExtend = 5;
	// 			params.alignmentBandwidth = 5;
	// 			params.maxCellsPerSlice = 10000;

	//.......................................TODO.......................................

	// params.graphFile = "graph.gfa";
	// 	std::vector<std::string> Reads;
	// 	Reads.push_back("reads.fa");
	// 	params.fastqFiles = Reads;
	// 	params.outputCorrectedFile = "reads.corrected.fa";
	// 	params.numThreads = 24;

	//........................................Test........................................

	// void Print_help(hifiasm_opt_t* asm_opt);
	// void Print_version();
	init_opt(&chat_opt);
	PRINT_LINE_FUNC();
	if (!Dechat_command(argc, argv, &chat_opt, &asm_opt))
		return 0;
	PRINT_LINE_FUNC();
	// //.......................................Round 1.......................................
	std::cout << "asm_opt->output_file_name:" << asm_opt.output_file_name << std::endl;
	std::cout << "chat_opt.fast:" << chat_opt.fast << std::endl;
	if (chat_opt.fast)
	{
		try
		{
			char *myArgv[] = {"dechat", "-in", chat_opt.read_file_names, "-kmer-size", strdup(std::to_string(chat_opt.k_mer_length).c_str()), "-abundance-min", strdup(std::to_string(chat_opt.abundance_min).c_str()), "-nb-cores", strdup(std::to_string(chat_opt.thread_num).c_str())};
			int myArgc = sizeof(myArgv) / sizeof(char *);
			dechat_1().run(myArgc, myArgv);

			// std::string inputFilename = "reads.unitigs.fa";
			// std::ifstream inputFile(inputFilename);
			// std::string outputFilename = "graph.gfa";
			// std::ofstream outputFile(outputFilename);

			// std::string line;
			// std::string name, segment;
			// std::vector<std::string> optional, links;
			// bool firstLine = false;
			// while (std::getline(inputFile, line)){

			// }

			std::string command = "python /home/work/yichen/tools/bcalm-1/scripts/convertToGFA.py reads.unitigs.fa graph.gfa " + std::to_string(chat_opt.k_mer_length) + " -s";
			int status = system(command.c_str());
			// //.......................GraphAligner......................................................//
			AlignerParams params;
			std::vector<std::string> Reads;
			Reads.push_back("reads.fa");
			params.fastqFiles = Reads;
			params.graphFile = "graph.gfa";
			params.outputGAMFile = "";
			params.outputJSONFile = "";
			params.outputGAFFile = "";
			params.outputCorrectedFile = "recorrected.fa";
			params.outputCorrectedClippedFile = "";
			params.numThreads = chat_opt.thread_num;
			params.alignmentBandwidth = 5;
			params.dynamicRowStart = false;
			params.maxCellsPerSlice = 10000;
			params.verboseMode = false;
			params.mxmLength = 20;
			params.mumCount = 0;
			params.memCount = 0;
			params.seederCachePrefix = "";
			params.selectionECutoff = -1;
			params.compressCorrected = false;
			params.compressClipped = false;
			params.minimizerSeedDensity = 5;
			params.minimizerLength = 19;
			params.minimizerWindowSize = 30;
			params.seedClusterMinSize = 1;
			params.minimizerDiscardMostNumerousFraction = 0.0002;
			params.maxClusterExtend = 5;
			params.preciseClippingIdentityCutoff = 0.66;
			params.Xdropcutoff = 50;
			params.DPRestartStride = 0;
			params.multimapScoreFraction = 0.9;
			params.cigarMatchMismatchMerge = false;
			params.minAlignmentScore = 0;
			params.hpcCollapse = false;
			params.includeCigar = true;
			params.clipAmbiguousEnds = -1;
			params.maxTraceCount = 10;
			params.overlapIncompatibleCutoff = 0.3;
			params.realignFile = "";
			params.uniqueMemBonusFactor = 1.0;
			params.lowMemoryMEMIndexConstruction = false;
			params.MEMindexUsesWaveletTree = true;
			params.MEMwindowsize = 0;
			params.useDiploidHeuristic = false;
			params.diploidHeuristicCacheFile = "";
			Aamain(argc, argv);
			alignReads(params);
			chat_opt.output_dir_ec = dirname(std::string(chat_opt.read_file_names)) + "recorrected.fa";
		}
		catch (Exception &e)
		{
			std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
			return EXIT_FAILURE;
		}
		std::cout << "End dechat" << std::endl;
	}
	else
	{
		correct_round1(&chat_opt);
	}

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