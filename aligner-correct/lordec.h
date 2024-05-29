#ifndef LORDEC_H
#define LORDEC_H
#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <stack>
#include <getopt.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <gatb/gatb_core.hpp>
#define GATB_V141

#define MIN_KMER_LEN 4
#define MIN_SOLID_THR 1

#define MIN_ERROR_RATE 0.0
#define MAX_ERROR_RATE 0.5

#define MIN_TRIALS 1
#define MAX_TRIALS 100

#define MIN_NB_BRANCH 1
#define MAX_NB_BRANCH 10000

#define MIN_THREADS 0
#define MAX_THREADS 64


// other constants
//#define MAX_READ_LEN 500000
#define MAX_PATH_LEN 5000

// default values for correction
#define DEF_ERROR_RATE 0.40
#define DEF_MAX_BRANCH 200
#define DEF_TRIALS 5
#define DEF_THREADS 0
#define DEF_MAX_ABUNDANCE 2147483647

// alignment scores
#define ALIGN_MATCH 1
#define ALIGN_MISMATCH -3
#define ALIGN_INDEL -2

// Constants for trimming and splitting
#define DEF_MIN_REGION_LG 100

// Constants for statistics
#define STAT_FOUND 0
#define STAT_FOUND_LEN1 1
#define STAT_TOOLONG 2
#define STAT_EXPLOSION 3
#define STAT_NOPATH 4

#define STAT_TAIL "TAIL"
#define STAT_END2END "END2END"
#define STAT_GAPEXTEND "GAPEXTEND"

// File extensions
#define HDF5_EXT ".h5"
#define FASTA_EXT ".fa"
#define FASTQ_EXT ".fq"

// Utilities 

// min macro
#define MIN(a,b) ((a) < (b) ? (a) : (b))

// max macro
#define MAX(a,b) ((a) > (b) ? (a) : (b))





#define TO_STR2(x) #x
#define STR(x) TO_STR2(x)
#define LORDECSTRVERSION STR(LORDECVERSION)
#define GATBSTRVERSION STR(GATBVERSION)

// FLAGS
bool DEBUG_l = true;

bool DEBUG2 = false;

// Global variables for correction parameters
float max_error_rate = DEF_ERROR_RATE;
int max_branch = DEF_MAX_BRANCH;
int max_trials = DEF_TRIALS;
int threads = DEF_THREADS;

int strict_mode = 0;

#ifndef GATB_V110
static const int KSIZE_4 = KMER_SPAN(3);
#endif


static const int progressBarWidth = 70;

// usage
std::string usageString = "\n\n-i|--long_reads <long read FASTA/Q file>\n-2|--short_reads <short read FASTA/Q file(s)>\n-k|--kmer_len <k-mer size>\n-o|--corrected_read_file <output reads file>\n-s|solid_threshold <solid k-mer abundance threshold>\n[-t|--trials <number of paths to try from a k-mer>]\n[-b|--branch <maximum number of branches to explore>]\n[-e|--errorrate <maximum error rate>]\n[-T|--threads <number of threads>]\n[-S|--stat_file <out statistics file>]\n[-c|--complete_search]\n[-a|--abundance-max <abundance max threshold for k-mers>]\n[-O|--out-tmp <GATB graph creation temporary files directory>]\n[-p|--progress]\n[-g|--graph_named_like_output]\n";

void usage(std::string prog) {
  std::cerr << "Usage :" << std::endl << std::endl;
  std::cerr << prog << usageString << std::endl;
}

std::string dirname(std::string source) {
  source.erase(std::find(source.rbegin(), source.rend(), '/').base(), source.end());
  return source;
}

std::string basename(std::string pathname) {
  std::string dirn = dirname(pathname);
  size_t pos = pathname.find(dirn);
  size_t len = dirn.length();
  const std::string empty = "";
  return pathname.replace(pos, len, empty);
}

//////////////////////////////////////////////////
// complement reverse a dna seq
void reverse(char *target, char *source, int len) {
  for(int i = 0; i < len; i++) {
    switch(source[len-i-1]) {
      case 'a':
        target[i] = 't';
        break;
      case 'c':
        target[i] = 'g';
        break;
      case 'g':
        target[i] = 'c';
        break;
      case 't':
        target[i] = 'a';
        break;
      case 'n':
        target[i] = 'n';
        break;
      case 'A':
        target[i] = 'T';
        break;
      case 'C':
        target[i] = 'G';
        break;
      case 'G':
        target[i] = 'C';
        break;
      case 'T':
        target[i] = 'A';
        break;
      case 'N':
        target[i] = 'N';
        break;
    }
  }
  target[len] = '\0';
}
//////////////////////////////////////////////////
// string manipulation lowercase - uppercase
void copy_lower_case(char *target, char *source, int len) {
  for(int i = 0; i < len; i++) {
    target[i] = tolower(source[i]);
  }
}

void copy_upper_case(char *target, char *source, int len) {
  for(int i = 0; i < len; i++) {
    target[i] = toupper(source[i]);
  }
}

//////////////////////////////////////////////////
// STACK class
class stack_element {
  public:

    Node node;
    int pos;
    char nt;

    stack_element(Node n, int p, char c);
};


stack_element::stack_element(Node n, int p, char c) {
  node=n;pos=p;nt=c;
}

//////////////////////////////////////////////////
// Attempt to correct the head/tail of a read. Find a correction for
// the longest part of the head/tail that can be corrected given the
// maximum error threshold. Of the alternative possible paths return
// the one with the smallest edit distance as compared to the original
// head/tail of the read.
//
// Return the length of the found path, -1 if no path was found
// int extend(Graph graph, int kmer_len,  Node begin, char *extended_path, int *max_b,
// 	   char *read_part, int part_len,
// 	   int dp[][MAX_PATH_LEN+1], int *best_ed, int *best_part_len) {
int extend(Graph graph, int kmer_len,  Node begin, char *extended_path, int *max_b,
    char *read_part, int part_len,
    int *dp, int *best_ed, int *best_part_len) {

  std::stack <stack_element *> nodes;
  char path[MAX_PATH_LEN+1];
  // char *path = extended_path;
  int best_len = -1;
  *best_ed = max_error_rate*(part_len+kmer_len)+1;
  *best_part_len=0;

  // store the neighbors of begin node in the stack

#if defined(GATB_V110)
  Graph::Vector<Edge> neighbors = graph.successors<Edge>(begin);
#elif defined(GATB_V130) || defined(GATB_V140) || defined(GATB_V141)
  GraphVector<Edge> neighbors = graph.successorsEdge(begin);
#else
  Graph::Vector<Edge> neighbors = graph.successorsEdge(begin);
#endif

  //  if (neighbors.size() == 1)
  for (int i = neighbors.size()-1; i>= 0; i--) {
    nodes.push(new stack_element(neighbors[i].to, 1, ascii(neighbors[i].nt)));
  }

  while(nodes.size() > 0 && *max_b > 0) {
    stack_element *current = nodes.top();
    nodes.pop();
    Node cnode = current->node;
    int pos = current->pos;
    path[pos-1] = current->nt;
    delete current;

    // compute the current row of the DP matrix and save its minimum
    dp[pos*(MAX_PATH_LEN+1)+0] = pos;
    // dp[pos][0] = pos;
    int min = MAX_PATH_LEN;
    int mink = -1;
    for(int k=1; k <= part_len; k++) {
      // if (path[pos-1] == read_part[k-1]) {
      // 	dp[pos][k] = MIN(dp[pos-1][k-1], MIN(dp[pos-1][k]+1, dp[pos][k-1]+1));
      // } else {
      // 	dp[pos][k] = MIN(dp[pos-1][k-1]+1, MIN(dp[pos-1][k]+1, dp[pos][k-1]+1));
      // }
      // if (dp[pos][k] < min) {
      // 	min = dp[pos][k];
      // 	mink = k;
      // }

      if (path[pos-1] == read_part[k-1]) {
        dp[pos*(MAX_PATH_LEN+1)+k] = MIN(dp[(pos-1)*(MAX_PATH_LEN+1)+k-1], MIN(dp[(pos-1)*(MAX_PATH_LEN+1)+k]+1, dp[pos*(MAX_PATH_LEN+1)+k-1]+1));
      }
      else {
        dp[pos*(MAX_PATH_LEN+1)+k] = MIN(dp[(pos-1)*(MAX_PATH_LEN+1)+k-1]+1, MIN(dp[(pos-1)*(MAX_PATH_LEN+1)+k]+1, dp[pos*(MAX_PATH_LEN+1)+k-1]+1));
      }
      if (dp[pos*(MAX_PATH_LEN+1)+k] < min) {
        min = dp[pos*(MAX_PATH_LEN+1)+k];
        mink = k;
      }
    }

    // Check if this path can eventually be under the error threshold
    if (min <= max_error_rate*(part_len+kmer_len)+1) {
      // Is it currently under the error threshold?
      if ((mink+kmer_len)*max_error_rate >= min) {
        // Is it better than what we have found before?
        if (mink > *best_part_len || (mink == *best_part_len && min < *best_ed)) {
          best_len = pos;
          *best_part_len = mink;
          *best_ed = min;
          memcpy(extended_path, path, pos);
        }
      }
#if defined(GATB_V110)
      Graph::Vector<Edge> neighbors = graph.successors<Edge>(cnode);
#elif defined(GATB_V130) || defined(GATB_V140) || defined(GATB_V141)
      GraphVector<Edge> neighbors = graph.successorsEdge(cnode);
#else
      Graph::Vector<Edge> neighbors = graph.successorsEdge(cnode);
#endif

      //      if (neighbors.size() == 1)
      for (int i = neighbors.size()-1; i >= 0; i--) {
        Edge e = neighbors[i];
        nodes.push(new stack_element(e.to, pos+1, ascii(e.nt)));
      }
      if (neighbors.size() == 0) {
        (*max_b)--;
      }
    }
    else {
      (*max_b)--;
    }
  }

  // Empty the stack and free memory
  while (nodes.size() > 0) {
    stack_element *current = nodes.top();
    nodes.pop();

    delete current;
  }

  // Check for strict mode
  if (strict_mode && *max_b <= 0) {
    return -1;
  }

  // Compute the best correcting extension
  int align_dp1[MAX_PATH_LEN];
  int align_dp2[MAX_PATH_LEN];

  int max=0, maxi=0, maxj=0;

  // Initialize the first row
  for (int i = 0; i <= best_len; i++) {
    align_dp1[i] = i*ALIGN_INDEL;
  }

  for (int j = 1; j <= *best_part_len; j++) {
    align_dp2[0] = j*ALIGN_INDEL;
    for(int i = 1; i <= best_len; i++) {
      if (read_part[j-1] == extended_path[i-1]) { // match
        align_dp2[i] = MAX(align_dp1[i-1] + ALIGN_MATCH, MAX(align_dp1[i] + ALIGN_INDEL, align_dp2[i-1] + ALIGN_INDEL));
      }
      else { // mismatch
        align_dp2[i] = MAX(align_dp1[i-1] + ALIGN_MISMATCH, MAX(align_dp1[i] + ALIGN_INDEL, align_dp2[i-1] + ALIGN_INDEL));
      }
      if (align_dp2[i] >= max) {
        max = align_dp2[i];
        maxi = i;
        maxj = j;
      }
    }

    if (j+1 <= *best_part_len) {
      j++;
      align_dp1[0] = j*ALIGN_INDEL;
      for(int i = 1; i <= best_len; i++) {
        if (read_part[j-1] == extended_path[i-1]) { // match
          align_dp1[i] = MAX(align_dp2[i-1] + ALIGN_MATCH, MAX(align_dp2[i] + ALIGN_INDEL, align_dp1[i-1] + ALIGN_INDEL));
        }
        else { // mismatch
          align_dp1[i] = MAX(align_dp2[i-1] + ALIGN_MISMATCH, MAX(align_dp2[i] + ALIGN_INDEL, align_dp1[i-1] + ALIGN_INDEL));
        }
        if (align_dp1[i] >= max) {
          max = align_dp1[i];
          maxi = i;
          maxj = j;
        }
      }
    }
  }

  best_len = maxi;
  *best_part_len = maxj;

  extended_path[best_len] = '\0';
  return best_len;
}

  //////////////////////////////////////////////////

  /* Find the best path (smallest edit distance as compared against the
     read) between begin and end allowing at most max_branch branches to
     be explored.

     Return the length of the path, -1 if the path was not found */
  // int best_path(Graph graph, Node begin, Node end, char *best_path, int *max_b,
  // 	      char *read_part, int part_len, int dp[][MAX_PATH_LEN+1], int *best_ed) {
int best_path(Graph graph, Node begin, Node end, char *best_path, int *max_b,
    char *read_part, int part_len, int *dp, int *best_ed) {

  std::stack <stack_element *> nodes;
  char path[MAX_PATH_LEN+1];
  *best_ed = max_error_rate*part_len+1;
  int best_len = -1;
  int i;

  // store the neighbors of begin node in the stack
#if defined(GATB_V110)
  Graph::Vector<Edge> neighbors = graph.successors<Edge>(begin);
#elif defined(GATB_V130) || defined(GATB_V140) || defined(GATB_V141)
  GraphVector<Edge> neighbors = graph.successorsEdge(begin);
#else
  Graph::Vector<Edge> neighbors = graph.successorsEdge(begin);
#endif
  for (i = neighbors.size()-1; i>= 0; i--) {
    nodes.push(new stack_element(neighbors[i].to, 1, ascii(neighbors[i].nt)));
  }

  while (nodes.size() > 0 && *max_b > 0) {
    stack_element *current = nodes.top();
    nodes.pop();
    Node cnode = current->node;
    int pos = current->pos;
    path[pos-1] = current->nt;

    delete current;

    if (pos >= MAX_PATH_LEN) {
      std::cout << "Too long path searched for" << std::endl;
      exit(EXIT_FAILURE);
    }

    // compute the current row of the DP matrix and save its minimum
    // dp[pos][0] = pos;
    // int min = MAX_PATH_LEN;
    // for(int k=1; k <=part_len; k++) {
    //   if (path[pos-1] == read_part[k-1]) {
    // 	dp[pos][k] = MIN(dp[pos-1][k-1], MIN(dp[pos-1][k]+1, dp[pos][k-1]+1));
    //   } else {
    // 	dp[pos][k] = MIN(dp[pos-1][k-1]+1, MIN(dp[pos-1][k]+1, dp[pos][k-1]+1));
    //   }
    //   if (dp[pos][k] < min)
    // 	min = dp[pos][k];
    // }

    dp[pos*(MAX_PATH_LEN+1)+0] = pos;
    int min = MAX_PATH_LEN;
    for(int k=1; k <=part_len; k++) {
      if (path[pos-1] == read_part[k-1]) {
        dp[pos*(MAX_PATH_LEN+1)+k] = MIN(dp[(pos-1)*(MAX_PATH_LEN+1)+k-1], MIN(dp[(pos-1)*(MAX_PATH_LEN+1)+k]+1, dp[pos*(MAX_PATH_LEN+1)+k-1]+1));
      }
      else {
        dp[pos*(MAX_PATH_LEN+1)+k] = MIN(dp[(pos-1)*(MAX_PATH_LEN+1)+k-1]+1, MIN(dp[(pos-1)*(MAX_PATH_LEN+1)+k]+1, dp[pos*(MAX_PATH_LEN+1)+k-1]+1));
      }
      if (dp[pos*(MAX_PATH_LEN+1)+k] < min)
        min = dp[pos*(MAX_PATH_LEN+1)+k];
    }

    // Can the current path eventually be extended to be better that
    // the best we have found so far?
    if (min < *best_ed) {
      // Have we reached the end node?
      if (cnode == end) {
        (*max_b)--;
        // int ed = dp[pos][part_len];
        int ed = dp[pos*(MAX_PATH_LEN+1)+part_len];
        // Is this path better than the previous one we have found?
        if (ed < *best_ed) {
          *best_ed = ed;
          memcpy(best_path, path, pos);
          best_len = pos;
        }
      }
      else {

#if defined(GATB_V110)
        Graph::Vector<Edge> neighbors = graph.successors<Edge>(cnode);
#elif defined(GATB_V130) || defined(GATB_V140) || defined(GATB_V141)
        GraphVector<Edge> neighbors = graph.successorsEdge(cnode);
#else
        Graph::Vector<Edge> neighbors = graph.successorsEdge(cnode);
#endif
        for(i=neighbors.size()-1; i>=0; i--) {
          Edge e = neighbors[i];
          nodes.push(new stack_element(e.to, pos+1, ascii(e.nt)));
        }
        if (neighbors.size() == 0) {
          (*max_b)--;
        }
      }
    }
    else {
      (*max_b)--;
    }
  }

  // Empty the stack and free memory
  while (nodes.size() > 0) {
    stack_element *current = nodes.top();
    nodes.pop();

    delete current;
  }

  if (strict_mode && *max_b <= 0) {
    return -1;
  }
  else {
    return best_len;
  }
}
  //////////////////////////////////////////////////

#define DEADEND 0
#define OUTBRANCHING -1
#define TOOLONG -3
#define OTHER -5

// A function to explore simple path from a starting node, up to a
// maximal path length, without paying attention to in-branching nodes
// along the path. Direction always forward; use reverse node to go in
// the other direction
//
// return 0 for a dead-end -1 for out branching ;-3 path too long; -4 cycle ???
int traverseForwardUniSimplePath (Graph graph, Node & node, int limitLg, char *path, int *pathLg) {
  *pathLg = 0;
  path[0] = '\0';

  while (true) {
#if defined(GATB_V110)
    Graph::Vector<Edge> neighbors = graph.successors<Edge>(node);
#elif defined(GATB_V130) || defined(GATB_V140) || defined(GATB_V141)
    GraphVector<Edge> neighbors = graph.successorsEdge(node);
#else
    Graph::Vector<Edge> neighbors = graph.successorsEdge(node);
#endif
    if (neighbors.size() == 0) // dead end
      return DEADEND;
    if ( neighbors.size() > 1 ) // out branching node
      return OUTBRANCHING;
    // invariant neighbors contains a single neighbor
    if ( *pathLg + 1 > limitLg )
      return TOOLONG;
    else {
      // the 0 element of the vector contains the edge to the only out going neighbor
      path[*pathLg] = ascii(neighbors[0].nt);
      path[*pathLg+1] = '\0';
      (*pathLg)++;
      node = neighbors[0].to;
    }
  }

  // We exit the loop by a return so we'll never end up here.
  return OTHER;
}

//////////////////////////////////////////

// Some statistics
long long int path_found=0;
long long int path_len1=0;
long long int path_nopath=0;
long long int path_explosion=0;
long long int path_toolong=0;

//////////////////////////////////////////////////

int correct_one_read(Sequence seq, char *corrected, Graph graph, IFile *statFile, ISynchronizer *syncStat, int kmer_len, int max_read_len) {
  int stats[10];
  char path[MAX_PATH_LEN+1];
  char path2[MAX_PATH_LEN+1];
  int corrected_len = 0;
  char buffer[MAX_PATH_LEN+1];
  char *read = seq.getDataBuffer();
  int read_len = seq.getDataSize();
  int *dp;
  // int dp[MAX_PATH_LEN+1][MAX_PATH_LEN+1];

  // Allocate memory for dp table
  // std::cout<<"Allocate memory for dp table"<<std::endl;
  dp = new int[(MAX_PATH_LEN+1)*(MAX_PATH_LEN+1)];

  read[read_len] = '\0';

#ifdef OLD_GATB
  Kmer<>::Model model(kmer_len);
  Kmer<>::Model::Iterator itKmer(model);
#else
  Kmer<KSIZE_4>::ModelDirect model(kmer_len);
  Kmer<KSIZE_4>::ModelDirect::Iterator itKmer(model);
#endif
  // Model<LargeInt<1> > model(kmer_len);
  // Model<LargeInt<1> >::Iterator itKmer(model);

  // Initialize DP matrix
  // std::cout<<"Initialize DP matrix"<<std::endl;
  for(int k=0; k <= (1+max_error_rate)*read_len+1; k++) {
    // dp[0][k] = k;
    dp[0*(MAX_PATH_LEN+1)+k] = k;
  }

  // We set the data from which we want to extract kmers.
  itKmer.setData(seq.getData());

  // We iterate the kmers.
  // First enumerate all positions with solid k-mers
  int pos = 0;
  int num_solid = 0;
  int *solid = new int[max_read_len];
  // std::cout<<"enumerate all positions with solid k-mers"<<std::endl;
  for (itKmer.first(); !itKmer.isDone(); itKmer.next())   {
#ifdef OLD_GATB
    Node node = graph.buildNode(Data((char *)(model.toString(*itKmer).c_str())));
#else
    Node node = graph.buildNode(Data((char *)(model.toString(itKmer->value()).c_str())));
#endif
    // The graph has false positives -> filtering out nodes without incoming or outgoing edges helps
    if (graph.contains(node) && graph.indegree(node) >= 1 && graph.outdegree(node)>=1) {
      solid[num_solid++] = pos;
    }
    pos++;
  }
  // std::cout<<"read_len"<<read_len<<std::endl;
  // Do we need this? (May be something strange with seq.getDataSize()?)
  read_len=strlen(read);
  // std::cout<<"num_solid"<<num_solid<<std::endl;
  if (num_solid == 0) {
    // No solid k-mers -> just copy the read as it is
    copy_lower_case(&corrected[corrected_len], &read[0], read_len);
    corrected_len += read_len;
    delete [] solid;
    delete [] dp;
    return corrected_len;
  }
  // std::cout<<"copy_lower_case"<<std::endl;
  if (solid[0] > 0) {
    // There is a head to correct
    if ((1.0+max_error_rate)*solid[0] < MAX_PATH_LEN) {

      // Attempt to correct the head of the read
      Node n = graph.buildNode(Data(&read[solid[0]]));
      n = graph.reverse(n);
      reverse(buffer, &read[0], solid[0]);

#ifdef DEBUG2

        std::cout << "Searching for path from " << graph.toString(n)
        << " with max length " << (1.0 + max_error_rate)*(solid[0]) <<  std::endl;
#endif

      int ed;
      int best_part_len;

      // for(int k=0; k <= solid[0]; k++) {
      //   dp[0][k] = k;
      // }


      int max_b = max_branch;
      int len = extend(graph, kmer_len, n, path, &max_b, buffer, solid[0], dp, &ed, &best_part_len);

      if (len >= 0) {
        path[len] = '\0';
#ifdef DEBUG2
        std::cout << "Found path of length " << len << " best_part_len " << best_part_len << std::endl;
        std::cout << path << std::endl;
        std::cout << buffer << std::endl;
#endif
        copy_lower_case(&corrected[corrected_len], &read[0], solid[0]-best_part_len);
        corrected_len += solid[0]-best_part_len;

        reverse(buffer, path, len);
        copy_upper_case(&corrected[corrected_len], buffer, len);
        copy_upper_case(&corrected[corrected_len+len], &read[solid[0]], kmer_len);
        corrected_len += len+kmer_len;

        if (statFile != NULL) {
          LocalSynchronizer local(syncStat);
          stats[0] = solid[0];
          stats[1] = STAT_FOUND;
          stats[2] = len;
          stats[3] = ed;
          statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_TAIL);
        }

      }
      else {
        // No correction of the head was found
        copy_lower_case(&corrected[corrected_len], &read[0], solid[0]);
        copy_upper_case(&corrected[corrected_len+solid[0]], &read[solid[0]], kmer_len);
        corrected_len += (solid[0]+kmer_len);

        if (statFile != NULL) {
          if (max_b == 0) { // branch limit was exhausted before finding a path
            LocalSynchronizer local(syncStat);
            stats[0] = solid[0];
            stats[1] = STAT_EXPLOSION;
            stats[2] = 0;
            stats[3] = 0;
            statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_TAIL);
          }
          else { // There is really no path
            LocalSynchronizer local(syncStat);
            stats[0] = solid[0];
            stats[1] = STAT_NOPATH;
            stats[2] = 0;
            stats[3] = 0;
            statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_TAIL);
          }
        }
      }
    }
    else {
      // The head is too long: no correction was attempted (memory constraints)
      copy_lower_case(&corrected[corrected_len], &read[0], solid[0]);
      copy_upper_case(&corrected[corrected_len+solid[0]], &read[solid[0]], kmer_len);
      corrected_len += (solid[0]+kmer_len);
      if (statFile != NULL) {
        LocalSynchronizer local(syncStat);
        stats[0] = solid[0];
        stats[1] = STAT_TOOLONG;
        stats[2] = 0;
        stats[3] = 0;
        statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_TAIL);
      }
    }
  }
  else {
    // The first k-mer is solid
    copy_upper_case(&corrected[corrected_len], &read[solid[0]], kmer_len);
    corrected_len += kmer_len;
  }

  ////////////////////////////////////////////////////////
  // INNER REGION CORRECTION

  if (num_solid >= 2) {
    // There is an intermediate part to correct

    // Data structure for a correction for the region between two solid k-mers
    struct Path {
      int ed;      // Edit distance between the region and the corrected region
      int len;     // Length of the corrected region
      char *str;   // Corrected region
    };

    // We use the graph from boost library
    // Path graph: describes the corrections that we have found between the solid k-mers
    //   - Nodes: solid k-mers: identified by their position in the solid-array
    //   - Edges: Found paths between the solid k-mers in dbg

    // Edge in the path graph
    typedef std::pair<int, int> Edge;
    // Path graph: no properties attached to the nodes, Path struct (above) gives the properties attached to the edges
    typedef boost::adjacency_list < boost::listS, boost::vecS, boost::directedS, boost::no_property, Path > path_graph_t;
    typedef boost::graph_traits < path_graph_t >::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits < path_graph_t >::edge_descriptor edge_descriptor;

    // Reserve memory for edges and edge properties
    Edge *paths;
    struct Path *path_prop;
    paths = (Edge *)MemoryAllocatorStdlib::singleton().malloc(max_read_len*max_trials*sizeof(Edge));
    path_prop = (struct Path *)MemoryAllocatorStdlib::singleton().malloc(max_read_len*max_trials*sizeof(struct Path));

    int num_paths = 0;

    for (int i = 0; i < num_solid-1; i++) {
      int found = 0;
      if (statFile != NULL) {
        statFile->print("solid[i]=%d\n", solid[i]);
      }
      // Number of correction paths we have attempted to find for this k-mer
      int trials = 0;
      for(int j = i+1; j < num_solid && trials < max_trials; j++) {
        if (solid[j] == solid[i] + (j-i) ) { // run of solid k-mers (all adjacent)
          trials++;
          if (num_paths >= max_read_len*max_trials) {
            std::cout << "Not enough memory allocated" << std::endl;
            exit(EXIT_FAILURE);
          }
          paths[num_paths] = std::make_pair(i,j);
          path_prop[num_paths].ed = 0;
          path_prop[num_paths].len = solid[j]-solid[i];
          path_prop[num_paths].str = (char *)MemoryAllocatorStdlib::singleton().malloc((solid[j]-solid[i]+1)*sizeof(char));
          copy_upper_case(path_prop[num_paths].str, &read[solid[i]+kmer_len], solid[j]-solid[i]);
          num_paths++;
          found++;
          if (statFile != NULL) {
            LocalSynchronizer local(syncStat);
            path_len1++;
            stats[0] = solid[j]-solid[i];
            stats[1] = STAT_FOUND_LEN1;
            stats[2] = solid[j]-solid[i];
            stats[3] = 0;
            statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_END2END);
          }
        }
        else if (solid[j] - solid[i] < (1.0-max_error_rate)*kmer_len) { // solid k-mers are within the same k-mer: too near from eachother
          continue;
        }
        else if ((1.0+max_error_rate)*(solid[j] - solid[i]) >= MAX_PATH_LEN) { // solid k-mers are too far from eachother: not enough memory for the DP matrix
          if (statFile) {
            LocalSynchronizer local(syncStat);
            path_toolong++;
            stats[0] = solid[j] - solid[i];
            stats[1] = STAT_TOOLONG;
            stats[2] = 0;
            stats[3] = 0;
            statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_END2END);
          }
          break;
        }
        else { // solid k-mers are neither too near, nor to far: try find a path
          trials++;

          Node ni = graph.buildNode(Data(&read[solid[i]]));
          Node nj = graph.buildNode(Data(&read[solid[j]]));

          strncpy(buffer, &read[solid[i]+kmer_len], solid[j]-solid[i]);
          buffer[solid[j]-solid[i]] = '\0';

#ifdef DEBUG2
          std::cout << "Searching for path " << graph.toString(ni) << " -> " << graph.toString(nj)
            << " with max length " << (1.0 + max_error_rate)*(solid[j]-solid[i]) <<  std::endl;
#endif
          int max_b = max_branch;
          int ed = MAX_PATH_LEN;                       // computed edit distance value


          int len = best_path(graph, ni, nj, path, &max_b, buffer, solid[j]-solid[i], dp, &ed);
          if (len >= 0) {
            path[len] = '\0';
#ifdef DEBUG2
            std::cout << "Found path of length " << len << std::endl;
            std::cout << path << std::endl;
            std::cout << buffer << std::endl;
#endif
            if (num_paths >= max_read_len*max_trials) {
              std::cout << "Not enough memory allocated" << std::endl;
              exit(EXIT_FAILURE);
            }
            paths[num_paths] = std::make_pair(i,j); // insert a new edge in the path graph
            // record the prop of the found path
            path_prop[num_paths].ed= ed;
            path_prop[num_paths].len = len;
            path_prop[num_paths].str = (char *)MemoryAllocatorStdlib::singleton().malloc((len+1)*sizeof(char));
            copy_upper_case(path_prop[num_paths].str, path, len); // copy the path to be saved with the path
            num_paths++;
            found++;

            if (statFile != NULL) {
              LocalSynchronizer local(syncStat);
              stats[0] = solid[j]-solid[i];
              stats[1] = STAT_FOUND;
              stats[2] = len;
              stats[3] = ed;
              statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_END2END);
              path_found++;
            }
          }
          else { // no found path
#ifdef DEBUG2
            std::cout << "None found: " << max_b << std::endl;
#endif

            if (statFile != NULL) {
              LocalSynchronizer local(syncStat);
              if (max_b == 0) { // branch limit was exhausted before finding a path
                stats[0] = solid[j] - solid[i];
                stats[1] = STAT_EXPLOSION;
                stats[2] = 0;
                stats[3] = 0;
                statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_END2END);
                path_explosion++;
              }
              else {              // there was really no path
                stats[0] = solid[j] - solid[i];
                stats[1] = STAT_NOPATH;
                stats[2] = 0;
                stats[3] = 0;
                statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_END2END);
                path_nopath++;
              }
            }
          }
        }
      } // End of trials for k-mer at position i
      if (found == 0) {
        if (solid[i+1]-solid[i] > kmer_len) {
          // Attempt to extend from the left side up to half the gap
          Node n = graph.buildNode(Data(&read[solid[i]]));
          int l = (solid[i+1]-solid[i]-kmer_len)/2;

          if ((1.0+max_error_rate)*(l+kmer_len) >= MAX_PATH_LEN-1) {
            l = (int)(MAX_PATH_LEN/(1.0+max_error_rate))-2-kmer_len;
          }

#ifdef DEBUG2
          std::cout << "Extending from left of gap from " << graph.toString(n)
            << " with max length " << (1.0 + max_error_rate)*l <<  std::endl;
#endif
          strncpy(buffer, &read[solid[i]+kmer_len], l);
          buffer[l] = '\0';

          int ed = 0;
          int best_part_len = 0;
          int max_b = max_branch;
          int len = extend(graph, kmer_len, n, path, &max_b, buffer, l, dp, &ed, &best_part_len);

#ifdef DEBUG2
          if (len >= 0) {
            std::cout << "Found left extension of length " << len << " best_part_len " << best_part_len << " edit distance " << ed  << std::endl;
            std::cout << path << std::endl;
            std::cout << buffer << std::endl;
          }
#endif

          if (statFile != NULL) {
            if (len >= 0) {
              {
                LocalSynchronizer local(syncStat);
                stats[0] = l;
                stats[1] = STAT_FOUND;
                stats[2] = len;
                stats[3] = ed;
                statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_GAPEXTEND);
              }
            }
            else {
              if (max_b == 0) {  // Branching limit exhausted
                {
                  LocalSynchronizer local(syncStat);
                  stats[0] = l;
                  stats[1] = STAT_EXPLOSION;
                  stats[2] = 0;
                  stats[3] = 0;
                  statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_GAPEXTEND);
                }
              }
              else { // Really no path
                {
                  LocalSynchronizer local(syncStat);
                  stats[0] = l;
                  stats[1] = STAT_NOPATH;
                  stats[2] = 0;
                  stats[3] = 0;
                  statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_GAPEXTEND);
                }
              }
            }
          }

          // Attempt to extend the gap from the right up to half the gap
          Node n2 = graph.reverse(graph.buildNode(Data(&read[solid[i+1]])));

          // If the gap is of odd length add the extra base to this extension
          if ((solid[i+1]-solid[i]-kmer_len) % 2 == 1)
            l++;

#ifdef DEBUG2
          std::cout << "Extending from right of gap from " << graph.toString(n2)
            << " with max length " << (1.0 + max_error_rate)*l <<  std::endl;
#endif

          reverse(buffer, &read[solid[i+1]-l], l);
          buffer[l] = '\0';

          int ed2 = 0;
          int best_part_len2 = 0;
          int max_b2 = max_branch;
          int len2 = extend(graph, kmer_len, n2, path2, &max_b2, buffer, l, dp, &ed2, &best_part_len2);

#ifdef DEBUG2
          if (len2 >= 0) {
            std::cout << "Found right extension of length " << len2 << " best_part_len " << best_part_len2 << " edit distance " << ed2 << std::endl;
            std::cout << path2 << std::endl;
            std::cout << buffer << std::endl;
          }
#endif

          if (statFile != NULL) {
            if (len2 >= 0) {
              {
                LocalSynchronizer local(syncStat);
                stats[0] = l;
                stats[1] = STAT_FOUND;
                stats[2] = len2;
                stats[3] = ed2;
                statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_GAPEXTEND);
              }
            }
            else {
              if (max_b2 == 0) {  // Branching limit exhausted
                {
                  LocalSynchronizer local(syncStat);
                  stats[0] = l;
                  stats[1] = STAT_EXPLOSION;
                  stats[2] = 0;
                  stats[3] = 0;
                  statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_GAPEXTEND);
                }
              }
              else { // Really no path
                {
                  LocalSynchronizer local(syncStat);
                  stats[0] = l;
                  stats[1] = STAT_NOPATH;
                  stats[2] = 0;
                  stats[3] = 0;
                  statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_GAPEXTEND);
                }
              }
            }
          }


          if (num_paths >= max_read_len*max_trials) {
            std::cout << "Not enough memory allocated" << std::endl;
            exit(EXIT_FAILURE);
          }

          if (len < 0) {
            ed = 0;
            len = 0;
            best_part_len = 0;
          }
          if (len2 < 0) {
            ed2 = 0;
            len2 = 0;
            best_part_len2 = 0;
          }

          int uncorrected_len = solid[i+1]-solid[i]-kmer_len - best_part_len - best_part_len2;
          if (uncorrected_len < 0) {
            std::cout << "Overlapping gap extensions!" << std::endl;
            exit(EXIT_FAILURE);
          }

          paths[num_paths] = std::make_pair(i,i+1);
          path_prop[num_paths].ed = ed+ed2+uncorrected_len;
          path_prop[num_paths].len = len+uncorrected_len+len2+kmer_len;
          path_prop[num_paths].str = (char *)MemoryAllocatorStdlib::singleton().malloc((path_prop[num_paths].len+1)*sizeof(char));

          if (len > 0) {
            copy_upper_case(path_prop[num_paths].str, path, len); // Copy the left extension
          }
          if (uncorrected_len > 0) {
            copy_lower_case(&path_prop[num_paths].str[len], &read[solid[i] + kmer_len + best_part_len], uncorrected_len); // Copy the uncorrected part
          }
          if (len2 > 0) {
            reverse(buffer, path2, len2);
            copy_upper_case(&path_prop[num_paths].str[len+uncorrected_len], buffer, len2); // Copy the right extension
          }
          copy_upper_case(&path_prop[num_paths].str[len+uncorrected_len+len2], &read[solid[i+1]], kmer_len); // Copy the right k-mer

          num_paths++;

        }
        else {
          // Add a dummy edge if the kmers overlap so that at least
          // one path from 1st to last solid k-mer is found
          if (num_paths >= max_read_len*max_trials) {
            std::cout << "Not enough memory allocated" << std::endl;
            exit(EXIT_FAILURE);
          }
          paths[num_paths] = std::make_pair(i,i+1);
          path_prop[num_paths].ed = solid[i+1]-solid[i];
          path_prop[num_paths].len = solid[i+1]-solid[i];
          path_prop[num_paths].str = (char *)MemoryAllocatorStdlib::singleton().malloc((solid[i+1]-solid[i]+1)*sizeof(char));
          if (solid[i+1]-solid[i]-kmer_len > 0) {
            copy_lower_case(path_prop[num_paths].str, &read[solid[i]+kmer_len], solid[i+1]-solid[i]-kmer_len);
            copy_upper_case(&path_prop[num_paths].str[solid[i+1]-solid[i]-kmer_len], &read[solid[i+1]], kmer_len);
          }
          else {
            copy_upper_case(path_prop[num_paths].str, &read[solid[i]+kmer_len], solid[i+1]-solid[i]);
          }
          num_paths++;
        }
      }
    }

    path_graph_t path_graph(num_solid);

    for (int i = 0; i < num_paths; i++) {
      boost::add_edge(paths[i].first, paths[i].second, path_prop[i], path_graph);
    }

    std::vector<vertex_descriptor> p(num_vertices(path_graph));
    std::vector<int> d(num_vertices(path_graph));
    vertex_descriptor s = vertex(0, path_graph);

    // Find the shortest paths from the first solid k-mer
    boost::dijkstra_shortest_paths(path_graph, s,
        boost::weight_map(get(&Path::ed, path_graph)).
        predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, path_graph))).
        distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, path_graph))));

    // Backtrack the shortest path from the last solid k-mer to the first
    int *successor = new int[max_read_len];
    int curr = num_solid-1;

    while (curr != 0) {
      if (curr < 0 || curr >= num_solid) {
        std::cout << "Curr: " << curr << std::endl;
        exit(EXIT_FAILURE);
      }
      successor[p[curr]] = curr;
      curr = p[curr];
    }

    // Traverse the shortest path and correct the read
    while(curr != num_solid-1) {
      boost::graph_traits < path_graph_t>::edge_descriptor  e = boost::edge(curr, successor[curr], path_graph).first;
      if (path_graph[e].len < 0 || corrected_len + path_graph[e].len >= max_read_len) {
        std::cout << "Invalid path length in the path graph" << std::endl;
        exit(EXIT_FAILURE);
      }
      strncpy(&corrected[corrected_len], path_graph[e].str, path_graph[e].len);
      corrected_len += path_graph[e].len;
      curr = successor[curr];
    }

#ifdef DEBUG2
    std::cout << "Found paths:" << std::endl;

    boost::graph_traits < path_graph_t >::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(path_graph); ei != ei_end; ++ei) {
      boost::graph_traits < path_graph_t >::edge_descriptor e = *ei;
      boost::graph_traits < path_graph_t >::vertex_descriptor
        u = boost::source(e, path_graph), v = boost::target(e, path_graph);
      std::cout << u << "->" << v << " " << path_graph[e].ed << std::endl;
    }

    std::cout << "distances and parents:" << std::endl;
    boost::graph_traits < path_graph_t >::vertex_iterator vi, vend;
    for (boost::tie(vi, vend) = vertices(path_graph); vi != vend; ++vi) {
      std::cout << "distance(" << *vi << ") = " << d[*vi] << ", ";
      std::cout << "parent(" << *vi << ") = " << p[*vi] << std::
        endl;
    }
    std::cout << std::endl;
#endif

    // Free the memory of the path graph
    for(int i = 0; i < num_paths; i++) {
      MemoryAllocatorStdlib::singleton().free(path_prop[i].str);
    }

    MemoryAllocatorStdlib::singleton().free(paths);
    MemoryAllocatorStdlib::singleton().free(path_prop);
    delete [] successor;
  }

  // Attempt to correct the tail of the read
  if (solid[num_solid-1] < read_len-kmer_len) {
    if ((1.0+max_error_rate)*(read_len-solid[num_solid-1]-kmer_len) < MAX_PATH_LEN) {
      Node n = graph.buildNode(Data(&read[solid[num_solid-1]]));

#ifdef DEBUG2
      std::cout << "Searching for path from " << graph.toString(n)
        << " with max length " << (1.0 + max_error_rate)*(read_len-solid[num_solid-1]-kmer_len) <<  std::endl;
#endif
      strncpy(buffer, &read[solid[num_solid-1]+kmer_len], read_len-solid[num_solid-1]-kmer_len);
      buffer[read_len-solid[num_solid-1]-kmer_len] = '\0';

      int ed;
      int best_part_len;

      // for(int k=0; k <= read_len-solid[num_solid-1]-kmer_len; k++) {
      //   dp[0][k] = k;
      // }

      int max_b = max_branch;
      int len = extend(graph, kmer_len, n, path, &max_b, buffer, read_len-solid[num_solid-1]-kmer_len, dp, &ed, &best_part_len);
      if (len >= 0) {
        path[len] = '\0';
#ifdef DEBUG2
        std::cout << "Found path of length " << len << " best_part_len " << best_part_len << std::endl;
        std::cout << path << std::endl;
        std::cout << buffer << std::endl;
#endif
        copy_upper_case(&corrected[corrected_len], path, len);
        corrected_len += len;
        if (read_len - solid[num_solid-1] - kmer_len > best_part_len) {
          copy_lower_case(&corrected[corrected_len], &read[solid[num_solid-1]+kmer_len+best_part_len],
              read_len - solid[num_solid-1]-kmer_len-best_part_len);
          corrected_len += read_len - solid[num_solid-1]-kmer_len-best_part_len;
        }

        if (statFile != NULL) {
          LocalSynchronizer local(syncStat);
          stats[0] = read_len-solid[num_solid-1]-kmer_len;
          stats[1] = STAT_FOUND;
          stats[2] = len;
          stats[3] = ed;
          statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_TAIL);
        }
      }
      else {
        copy_lower_case(&corrected[corrected_len], &read[solid[num_solid-1]+kmer_len], read_len - solid[num_solid-1]-kmer_len);
        corrected_len += (read_len - solid[num_solid-1]-kmer_len);

        if (statFile != NULL) {
          if (max_b == 0) { // branch limit was exhausted before finding a path
            LocalSynchronizer local(syncStat);
            stats[0] = read_len-solid[num_solid-1]-kmer_len;
            stats[1] = STAT_EXPLOSION;
            stats[2] = 0;
            stats[3] = 0;
            statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_TAIL);
          }
          else { // There is really no path
            LocalSynchronizer local(syncStat);
            stats[0] = read_len-solid[num_solid-1]-kmer_len;
            stats[1] = STAT_NOPATH;
            stats[2] = 0;
            stats[3] = 0;
            statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_TAIL);
          }
        }
      }
    }
    else {
      copy_lower_case(&corrected[corrected_len], &read[solid[num_solid-1]+kmer_len], read_len - solid[num_solid-1]-kmer_len);
      corrected_len += (read_len - solid[num_solid-1]-kmer_len);
      if (statFile != NULL) {
        LocalSynchronizer local(syncStat);
        stats[0] = read_len-solid[num_solid-1]-kmer_len;
        stats[1] = STAT_TOOLONG;
        stats[2] = 0;
        stats[3] = 0;
        statFile->print("%d %d %d %d %s\n", stats[0], stats[1], stats[2], stats[3], STAT_TAIL);
      }
    }
    // else
    //   copy_lower_case(&corrected[corrected_len], &read[solid[num_solid-1]+kmer_len], read_len - solid[num_solid-1]-kmer_len);
    //   corrected_len += (read_len - solid[num_solid-1]-kmer_len);
  }

  delete [] solid;
  delete [] dp;

  return corrected_len;
}



//////////////////////////////////////////////////
// check wheter a file is present and readable 
bool is_readable( const std::string & file ) 
{ 
    std::ifstream fichier( file.c_str() ); 
    return !fichier.fail(); 
} 




int AA_main(int argc, char* argv[]);

#endif