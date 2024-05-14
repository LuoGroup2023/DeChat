#define __STDC_LIMIT_MACROS
#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include "CommandLines1.h"
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#define ko_no_argument 0
#define ko_required_argument 1
#define ko_optional_argument 2

static ko_longopt_t long_options1[] = {
    {"version", ko_required_argument, 58},
    {"fast", ko_no_argument, 100},
    {"dbg-gfa", ko_no_argument, 301},
    {"write-paf", ko_no_argument, 302},
    {"write-ec", ko_no_argument, 303},
    {"skip-triobin", ko_no_argument, 304},
    {"max-od-ec", ko_no_argument, 305},
    {"max-od-final", ko_no_argument, 306},
    {"ex-list", ko_required_argument, 307},
    {"ex-iter", ko_required_argument, 308},
    {"hom-cov", ko_required_argument, 309},
    {"pri-range", ko_required_argument, 310},
    {"lowQ", ko_required_argument, 312},
    {"min-hist-cnt", ko_required_argument, 313},
    {"h1", ko_required_argument, 314},
    {"h2", ko_required_argument, 315},
    {"enzyme", ko_required_argument, 316},
    {"b-cov", ko_required_argument, 317},
    {"h-cov", ko_required_argument, 318},
    {"m-rate", ko_required_argument, 319},
    {"primary", ko_no_argument, 320},
    {"t-occ", ko_required_argument, 321},
    {"seed", ko_required_argument, 322},
    {"n-perturb", ko_required_argument, 323},
    {"f-perturb", ko_required_argument, 324},
    {"n-hap", ko_required_argument, 325},
    {"n-weight", ko_required_argument, 326},
    {"l-msjoin", ko_required_argument, 327},
    {"purge-max", ko_required_argument, 328},
    {"dp-er", ko_required_argument, 330},
    {"max-kocc", ko_required_argument, 331},
    {"hg-size", ko_required_argument, 332},
    {"ul", ko_required_argument, 333},
    {"unskew", ko_no_argument, 334},
    {"kpt-rate", ko_required_argument, 335},
    {"ul-rate", ko_required_argument, 336},
    {"dbg-het-cnt", ko_no_argument, 337},
    {"ul-tip", ko_required_argument, 338},
    {"low-het", ko_no_argument, 339},
    {"s-base", ko_required_argument, 340},
    {"bin-only", ko_no_argument, 341},
    {"ul-round", ko_required_argument, 342},
    {"prt-raw", ko_no_argument, 343},
    {"integer-correct", ko_required_argument, 344},
    {"dbg-ovec", ko_no_argument, 345},
    {"path-max", ko_required_argument, 346},
    {"path-min", ko_required_argument, 347},
    {"trio-dual", ko_no_argument, 348},
    {"ul-cut", ko_required_argument, 349},
    {"dual-scaf", ko_no_argument, 350},
    {"scaf-gap", ko_required_argument, 351},
    {"sec-in", ko_required_argument, 352},
    {"somatic-cov", ko_required_argument, 353},
    // { "path-round",     ko_required_argument, 348},
    {0, 0, 0}};

void Print_help(chat_opt_t *chat_opt)
{
    fprintf(stderr, "Repeat and haplotype aware error correction in nanopore sequencing reads with DeChat\n");
    fprintf(stderr, "Usage: dechat [options] -o <output> -t <thread>  -i <reads> <...>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  Input/Output:\n");
    fprintf(stderr, "       -o STR       prefix of output files [%s]\n", chat_opt->outReadFile);
    fprintf(stderr,  "                   The output for the stage 1 of correction is \"recorrected.fa\", \n");
    fprintf(stderr, "                    The final corrected file is \"file name\".ec.fa;\n");
    fprintf(stderr, "       -t INT       number of threads [%d]\n", chat_opt->thread_num);
    fprintf(stderr, "       -h           show help information\n");
    fprintf(stderr, "       -v --version show version number\n");
    fprintf(stderr, "       -i           input reads file\n");
    fprintf(stderr, "       -k INT       k-mer length (must be <64) [%d]\n", chat_opt->k_mer_length);
    fprintf(stderr, "  Error correction stage 1 (dBG):\n");
    fprintf(stderr, "       -r1           set the maximal abundance threshold for a k-mer in dBG [%d]\n",chat_opt->abundance_min);
    fprintf(stderr, "  Error correction stage 2 (MSA):\n");
    fprintf(stderr, "       -r            round of correction in alignment [%d]\n",chat_opt->second_number_of_round);
    fprintf(stderr, "       -e            maximum allowed error rate used for filtering overlaps [%.2f]\n",chat_opt->max_ov_diff_ec);
}
// 初始化默认参数
void init_opt(chat_opt_t *chat_opt)
{
    // memset(chat_opt, 0, sizeof(hifichat_opt_t));
    /// chat_opt->flag = 0;
    // chat_opt->flag = HA_F_PARTITION;
    chat_opt->abundance_min = 2;
    chat_opt->coverage = -1;
    chat_opt->num_reads = 0;
    chat_opt->fast = 0;
    chat_opt->read_file_names = NULL;
    // chat_opt->outReadFile = (char*)(DEFAULT_OUTPUT);
    chat_opt->required_read_name = NULL;
    chat_opt->fn_bin_poy = NULL;
    // chat_opt->ar = NULL;
    chat_opt->thread_num = 1;
    chat_opt->k_mer_length = 21;
    chat_opt->k_mer_length_str = "21";
    chat_opt->ul_mer_length = 19;
    chat_opt->trans_mer_length = 31;
    chat_opt->mz_win = 51;
    chat_opt->ul_mz_win = 19;
    chat_opt->trans_win = 31;
    chat_opt->mz_rewin = 1000;
    chat_opt->ul_mz_rewin = 360;
    chat_opt->mz_sample_dist = 500;
    chat_opt->bf_shift = 37;
    chat_opt->max_kmer_cnt = 2000;
    chat_opt->high_factor = 5.0;
    chat_opt->max_ov_diff_ec = 0.04; // 最大错误率默认设置为4%，最低设置为2%
    chat_opt->max_ov_diff_final = 0.03;
    chat_opt->hom_cov = 20;
    chat_opt->het_cov = -1024;
    // chat_opt->max_n_chain = MIN_N_CHAIN;
    chat_opt->min_hist_kmer_cnt = 5;
    chat_opt->load_index_from_disk = 1;
    chat_opt->write_index_to_disk = 1;
    chat_opt->second_number_of_round = 3;
    chat_opt->adapterLen = 0;
    chat_opt->clean_round = 4;
    /// chat_opt->small_pop_bubble_size = 100000;
    chat_opt->small_pop_bubble_size = 0;
    chat_opt->large_pop_bubble_size = 10000000;
    chat_opt->min_drop_rate = 0.2;
    chat_opt->max_drop_rate = 0.8;
    chat_opt->max_hang_Len = 1000;
    chat_opt->max_hang_rate = 0.8;
    chat_opt->gap_fuzz = 1000;
    chat_opt->min_overlap_Len = 50;
    chat_opt->min_overlap_coverage = 0;
    chat_opt->max_short_tip = 3;
    chat_opt->max_short_ul_tip = 6;
    chat_opt->min_cnt = 2;
    chat_opt->mid_cnt = 5;
    chat_opt->purge_level_primary = 3;
    chat_opt->purge_level_trio = 0;
    chat_opt->purge_simi_rate_l2 = 0.75;
    chat_opt->purge_simi_rate_l3 = 0.55;
    chat_opt->trans_base_rate = 0.93;
    chat_opt->trans_base_rate_sec = 0.5;
    chat_opt->purge_overlap_len = 1;
    /// chat_opt->purge_overlap_len_hic = 50;
    chat_opt->recover_atg_cov_min = -1024;
    chat_opt->recover_atg_cov_max = INT_MAX;
    chat_opt->hom_global_coverage = -1;
    chat_opt->hom_global_coverage_set = 0;
    chat_opt->pur_global_coverage = -1;
    chat_opt->bed_inconsist_rate = 70;
    /// chat_opt->bub_mer_length = 3;
    chat_opt->bub_mer_length = 1000000;
    chat_opt->b_low_cov = 0;
    chat_opt->b_high_cov = -1;
    chat_opt->m_rate = 0.75;
    chat_opt->hap_occ = 1;
    chat_opt->polyploidy = 2;
    chat_opt->trio_flag_occ_thres = 60;
    chat_opt->seed = 11;
    chat_opt->n_perturb = 10000;
    chat_opt->f_perturb = 0.1;
    chat_opt->n_weight = 3;
    chat_opt->is_alt = 0;
    chat_opt->misjoin_len = 500000;
    chat_opt->scffold = 0;
    chat_opt->dp_min_len = 2000;
    chat_opt->dp_e = 0.0025;
    chat_opt->hg_size = -1;
    chat_opt->kpt_rate = -1;
    chat_opt->infor_cov = 3;
    chat_opt->s_hap_cov = 3;
    chat_opt->ul_error_rate = 0.2 /**0.15**/;
    chat_opt->ul_error_rate_low = 0.1;
    chat_opt->ul_error_rate_hpc = 0.2;
    chat_opt->ul_ec_round = 3;
    chat_opt->is_dbg_het_cnt = 0;
    chat_opt->is_low_het_ul = 0;
    chat_opt->is_base_trans = 1;
    chat_opt->is_read_trans = 1;
    chat_opt->is_topo_trans = 1;
    chat_opt->is_bub_trans = 1;
    chat_opt->bin_only = 0;
    chat_opt->ul_clean_round = 1;
    chat_opt->prt_dbg_gfa = 0;
    chat_opt->integer_correct_round = 0;
    chat_opt->dbg_ovec_cal = 0;
    chat_opt->min_path_drop_rate = 0.2;
    chat_opt->max_path_drop_rate = 0.6;
    chat_opt->hifi_pst_join = 1;
    chat_opt->ul_pst_join = 1;
    chat_opt->trio_cov_het_ovlp = -1;
    chat_opt->ul_min_base = 0;
    chat_opt->self_scaf = 0;
    chat_opt->self_scaf_min = 250000;
    chat_opt->self_scaf_reliable_min = 5000000;
    chat_opt->self_scaf_gap_max = 3000000;
    // chat_opt->sec_in = NULL;
    chat_opt->somatic_cov = -1;
}

int check_option(chat_opt_t *chat_opt)
{
    if (chat_opt->read_file_names == nullptr)
    {
        fprintf(stderr, "[ERROR] missing input: please specify a read file\n");
        return 0;
    }

    if (chat_opt->outReadFile == nullptr)
    {
        fprintf(stderr, "[ERROR] missing output: please specify the output name (-o)\n");
        return 0;
    }
    if (chat_opt->thread_num < 1)
    {
        fprintf(stderr, "[ERROR] the number of threads must be > 0 (-t)\n");
        return 0;
    }
    std::cout << chat_opt->second_number_of_round << std::endl;
    if (chat_opt->second_number_of_round < 1)
    {
        fprintf(stderr, "[ERROR] the number of rounds for correction must be > 0 (-r)\n");
        return 0;
    }
    return 1;
}

// 参数读入
int Dechat_command(int argc, char *argv[], chat_opt_t *chat_opt, hifiasm_opt_t *asm_opt)
{
    ketopt_t opt = KETOPT_INIT;
    int c;
    int option_index = 0;
    // PRINT_LINE_FUNC();
    while ((c = ketopt(&opt, argc, argv, 1, "hvt:i:k:o:r:r1:e:hifi", long_options1)) >= 0)
    {
        if (c == 'h')
        {
            Print_help(chat_opt);
            return 0;
        }
        else if (c == 'v'|| c == 58)
        {
            std::ifstream file("../version.txt");
            std::string content;
            std::getline(file, content);
            std::cout<<"DeChat version:" << content << std::endl;
            file.close();
            return 0;
        }
        else if (c == 't')
            chat_opt->thread_num = atoi(opt.arg);
        else if (c == 'k')
            chat_opt->k_mer_length = atoi(opt.arg);
        else if (c == 'i')
            chat_opt->read_file_names = opt.arg;
        else if (c == 'o')
            chat_opt->outReadFile = opt.arg, asm_opt->output_file_name = opt.arg;
        else if (c == 'r')
            chat_opt->second_number_of_round = atoi(opt.arg);
        else if (c == 'r1')
            chat_opt->abundance_min = atoi(opt.arg);
        else if (c == 'e')
            chat_opt->max_ov_diff_ec = strtod(opt.arg,NULL);
        else if (c== 'hifi')
            chat_opt->hifi = atoi(opt.arg);
    }
    // PRINT_LINE_FUNC();
    if (argc == 1)
    {
        Print_help(chat_opt);
        return 0;
    }
    // std::cout << chat_opt->second_number_of_round << std::endl;
    // std::cout << "asm_opt->output_file_name:" << asm_opt->output_file_name << std::endl;
    // chat_opt->fast = 1;
    return check_option(chat_opt);
}