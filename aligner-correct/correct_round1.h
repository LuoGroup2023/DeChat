#ifndef _TOOL_CORRECT_ROUND1_H_
#define _TOOL_CORRECT_ROUND1_H_
#include "utils.hpp"
#include <cstdlib>
#include "lordec.h"
#include "CommandLines.h"
#include "CommandLines1.h"
#include "Correct.h"
#include "ketopt.h"
#include <iostream>

void correct_round2(chat_opt_t *chat_opt, hifiasm_opt_t *asm_opt);
void correct_round1(chat_opt_t *chat_opt)
{
    std::cout << "correct_round1 thread:" << chat_opt->thread_num << std::endl;
    // PRINT_LINE_FUNC();
    // if (chat_opt->dBGFile != NULL)
    //     std::cout << "dBG use " << chat_opt->dBGFile << std::endl;
    std::string pacbioFile = chat_opt->read_file_names;
    BankFasta bsize(pacbioFile);
    BankFasta::Iterator itSeqSize(bsize);
    size_t max_read_len = 0;
    size_t seqSize;
    long long nbSeq = 0;
    for (itSeqSize.first(); !itSeqSize.isDone(); itSeqSize.next())
    {
        seqSize = itSeqSize->getDataSize();
        if (seqSize > max_read_len)
        {
            max_read_len = seqSize;
        }
        nbSeq++;
    }
    // PRINT_LINE_FUNC();
    max_read_len = max_read_len * 2;
    // std::cout << "max_read_len:" << max_read_len << std::endl;
    std::string ONTGraph;
    int comaPosition = std::string::npos;
    // PRINT_LINE_FUNC();
    if (chat_opt->dBGFile == NULL)
    {
        std::cout << "dBG use ONT reads" << chat_opt->read_file_names << std::endl;
        ONTGraph = chat_opt->read_file_names;
        std::cerr << "creating the graph from file(s): " << chat_opt->read_file_names << std::endl;
    }
    else
    {
        std::cout << "dBG use " << chat_opt->dBGFile << std::endl;
        std::cerr << "creating the graph from file(s): " << chat_opt->dBGFile << std::endl;
        ONTGraph = chat_opt->dBGFile;
    }

    comaPosition = ONTGraph.find(",");

    if (comaPosition != std::string::npos)
    {
        ONTGraph = ONTGraph.substr(0, comaPosition) + "_multi";
    }
    ONTGraph = "reads_k" + std::to_string(chat_opt->k_mer_length) + "_s" + "1" + ".h5";
    std::cout << "ONTGraph:" << ONTGraph << std::endl;
    std::string dirn = dirname(ONTGraph);
    int writeGraphPossible;
    if (dirn == "")
    {
        writeGraphPossible = access("./", W_OK);
    }
    else
    {
        writeGraphPossible = access(dirn.c_str(), W_OK);
    }
    // if it's not possible : we put/read it in the output directory
    if (writeGraphPossible != 0)
    {
        std::cout << "Impossible to write in " << dirname(ONTGraph) << "\n";
        ONTGraph = dirname(chat_opt->read_file_names) + basename(ONTGraph);
        std::cout << "Graph will be written in output directory : " << ONTGraph << "\n";
    }

    Graph graph;
    // PRINT_LINE_FUNC();
    if (DEBUG_l)
    {
    }
    try
    {
        // v106: open IBank from 1/ list of filenames 2/ a file of filenames
        IBank *b;
        if (chat_opt->dBGFile == NULL)
        {
            std::cout << "dBG use ONT reads" << std::endl;
            b = Bank::open(chat_opt->read_file_names);
        }
        else
        {
            std::cout << "dBG use " << chat_opt->dBGFile << std::endl;
            b = Bank::open(chat_opt->dBGFile);
        }

        std::string outTmpPath = "";
        // outTmpPath = outTmpPath + dirname(std::string(chat_opt->read_file_names)) + "recorrected.fa";
        // PRINT_LINE_FUNC();
        std::cout << outTmpPath.c_str() << "   ONTGraph:" << ONTGraph << std::endl;
        graph = Graph::create(b, (const char *)"%s -out %s -kmer-size %d -abundance-min %d -bloom cache -debloom original -debloom-impl basic -nb-cores %d -abundance-max 2147483647", outTmpPath.c_str(), ONTGraph.c_str(), chat_opt->k_mer_length, chat_opt->abundance_min, chat_opt->thread_num);
        // PRINT_LINE_FUNC();
        if (is_readable(ONTGraph))
        {
            std::cerr << "!!! file present : " << ONTGraph << std::endl;
        }
        else
        {
            std::cerr << "!!! file NOT present : " << ONTGraph << std::endl;
        }
    }
    catch (Exception &e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
        return;
    }

    if (DEBUG_l)
    {
        std::cerr << "graph created" << std::endl;
    }

    std::cout << graph.getInfo() << std::endl;

    IBank *ptrBankPB = NULL;

    try
    {

        ptrBankPB = Bank::open(pacbioFile);
        // std::cout << "ptrBankPB = Bank::open(pacbioFile) done" << std::endl;
    }
    catch (gatb::core::system::Exception &bankPBExc)
    {
        std::cout << "Error message PacBio bank " << bankPBExc.getMessage() << std::endl;
        return;
    }
    Iterator<Sequence> *itSeq = ptrBankPB->iterator();
    chat_opt->output_dir_ec = dirname(std::string(chat_opt->read_file_names)) + "recorrected.fa";
    BankFasta output(chat_opt->output_dir_ec); // argv[5]
    cout << "Found " << nbSeq << " reads.\n";
    cout << "Correcting reads...\n";
    IFile *statFile = NULL;
    ISynchronizer *syncStat = NULL;

    std::cout << "max_read_len:" << max_read_len << std::endl;
    // progress
    long long nbSeqProcessed = 0;
    ProgressManager *pmCorrect = new ProgressManager(nbSeq, "reads");

    // Access to the output file must be synchronized
    ISynchronizer *sync = System::thread().newSynchronizer();
    ISynchronizer *syncCount = System::thread().newSynchronizer();

    // Initialize a counter
    long long processedCount = 0;

    IDispatcher::Status status = Dispatcher(threads).iterate(itSeq, [&](const Sequence &seq)
                                                             {
    {
        // Increment the counter
        LocalSynchronizer local(syncCount);
        processedCount++;
    }

    // Check if we have processed 50000 reads to enable DEBUG2
    // if (processedCount == 30000)
    // // 498100
    // {
    //     // Enable DEBUG2
    //     DEBUG2 = true;
    // }
    // DEBUG2 = true;
    if (seq.getDataSize() >= chat_opt->k_mer_length)
    {
        char *read = new char[max_read_len];
        int read_len = 0;
        char *buffer = new char[max_read_len];
        if (seq.getDataSize() > max_read_len)
        {
            std::cout << "Too long read" << std::endl;
            exit(EXIT_FAILURE);
        }

        // Correct the read backward
        copy_upper_case(buffer, seq.getDataBuffer(), seq.getDataSize());
        buffer[seq.getDataSize()] = '\0';
        reverse(read, buffer, strlen(buffer));
        Sequence seq1(read);
        seq1._comment = seq.getComment();
        read_len = correct_one_read(seq1, buffer, graph, statFile, syncStat, chat_opt->k_mer_length, max_read_len);
        // Correct the read forward
        copy_upper_case(read, buffer, read_len);
        buffer[read_len] = '\0';
        reverse(buffer, read, read_len);
        Sequence seq2(buffer);
        seq2._comment = seq.getComment();
        read_len = correct_one_read(seq2, read, graph, statFile, syncStat, chat_opt->k_mer_length, max_read_len);

        read[read_len] = '\0';
        Sequence s(read);
        s._comment = seq.getComment();
        {
            LocalSynchronizer local(sync);
            output.insert(s);
        }

        // Conditional DEBUG2 output
        if (DEBUG2)
        {
            std::cout << "Found path of length " << read_len << std::endl;
            std::cout << "Buffer: " << buffer << std::endl;
            std::cout << "Read: " << read << std::endl;
        }

        delete[] read;
        delete[] buffer;
    } });

    std::cout << std::endl;
    delete pmCorrect;

    output.flush();

    delete ptrBankPB;
    delete sync;
    delete syncCount;

    if (statFile != NULL)
    {
        statFile->flush();
        delete statFile;

        delete syncStat;

        std::cout << "Path statistics:" << std::endl;
        std::cout << "Path found: " << path_found << std::endl;
        std::cout << "Path (of length 1) found: " << path_len1 << std::endl;
        std::cout << "No path found: " << path_nopath << std::endl;
        std::cout << "Combinatorial explosion: " << path_explosion << std::endl;
        std::cout << "K-mers too distant: " << path_toolong << std::endl;

        std::cout << std::endl
                  << "Total: " << (path_found + path_len1 + path_nopath + path_explosion + path_toolong) << std::endl;
    }
}

void correct_round2(chat_opt_t *chat_opt, hifiasm_opt_t *asm_opt)
{
    int ret;
    yak_reset_realtime();
    int c;
    init_opt1(asm_opt);
    std::cout << "successful" << std::endl;
    //.................传入参数..................
    asm_opt->num_reads = 1;
    asm_opt->read_file_names = new char *[asm_opt->num_reads];
    asm_opt->thread_num = chat_opt->thread_num;
    asm_opt->read_file_names[0] = strdup(chat_opt->output_dir_ec.c_str());

    //
    std::cout << asm_opt->read_file_names[0] << std::endl;
    std::string outputname = chat_opt->outReadFile;
    std::cout << "chat_opt->outReadFile:" << chat_opt->outReadFile << std::endl;
    asm_opt->output_file_name = chat_opt->outReadFile;
    std::cout << "asm_opt->output_file_name:" << asm_opt->output_file_name << std::endl;

    asm_opt->max_ov_diff_ec = chat_opt->max_ov_diff_ec;
    asm_opt->number_of_round = chat_opt->second_number_of_round;
    std::cout << "asm_opt->max_ov_diff_ec:" << asm_opt->max_ov_diff_ec << std::endl;
    std::cout << "asm_opt->number_of_round:" << asm_opt->number_of_round << std::endl;
    check_option1(asm_opt);

    ret = ha_assemble();

    return;
}

#endif