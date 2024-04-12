#ifndef Aligner_h
#define Aligner_h

#include <string>
#include <vector>
#include "AlignmentGraph.h"
#include "vg.pb.h"
#include "AlignmentSelection.h"

struct AlignerParams
{
	std::string graphFile;
	std::vector<std::string> fastqFiles;
	size_t numThreads;
	size_t alignmentBandwidth;
	bool dynamicRowStart;
	size_t maxCellsPerSlice;
	std::vector<std::string> seedFiles;
	std::string outputGAMFile;
	std::string outputJSONFile;
	std::string outputGAFFile;
	std::string outputCorrectedFile;
	std::string outputCorrectedClippedFile;
	bool verboseMode;
	size_t mxmLength;
	size_t mumCount;
	size_t memCount;
	std::string seederCachePrefix;
	double selectionECutoff;
	bool compressCorrected;
	bool compressClipped;
	size_t minimizerLength;
	size_t minimizerWindowSize;
	double minimizerSeedDensity;
	size_t seedClusterMinSize;
	double minimizerDiscardMostNumerousFraction;
	size_t maxClusterExtend;
	double preciseClippingIdentityCutoff;
	int Xdropcutoff;
	size_t DPRestartStride;
	double multimapScoreFraction;
	bool cigarMatchMismatchMerge;
	double minAlignmentScore;
	bool hpcCollapse;
	bool includeCigar;
	int clipAmbiguousEnds;
	double overlapIncompatibleCutoff;
	size_t maxTraceCount;
	std::string realignFile;
	double uniqueMemBonusFactor;
	bool lowMemoryMEMIndexConstruction;
	bool MEMindexUsesWaveletTree;
	size_t MEMwindowsize;
	bool useDiploidHeuristic;
	std::vector<size_t> diploidHeuristicK;
	std::string diploidHeuristicCacheFile;
};

void alignReads(AlignerParams params);
void replaceDigraphNodeIdsWithOriginalNodeIds(vg::Alignment& alignment, const AlignmentGraph& graph);
int Aamain(int argc, char** argv);
#endif
