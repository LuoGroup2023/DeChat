#ifndef __COMMAND1_LINE_PARSER__
#define __COMMAND1_LINE_PARSER__

#define __STDC_LIMIT_MACROS
#include <pthread.h>
#include <stdint.h>
#include <string>
#include "ketopt.h"
#include "CommandLines.h"
// #define DEFAULT_OUTPUT "hifiasm.corrected"
#define HA_VERSION "0.19.8-r603"
// #define PRINT_LINE_FUNC() printf("Line %d in function %s\n", __LINE__, __func__)
typedef struct {
    //...................................Test.......................................
    char* read_file_names;
    int thread_num;
    int k_mer_length;
    char* k_mer_length_str;
    char* outReadFile;
    int second_number_of_round;
    int fast;
    int abundance_min;
    int hifi;
    //...................................hifiasm....................................
    int flag;
    int num_reads;
    // char** read_file_names;
    // char* output_file_name;
    char* required_read_name;
	char *fn_bin_yak[2];
	char *fn_bin_list[2];
    char *fn_bin_poy;
	char *extract_list;
    int extract_iter;
    int ul_mer_length;
    int trans_mer_length;
    int bub_mer_length;
	int mz_win;
    int ul_mz_win;
    int trans_win;
    int mz_rewin;
    int ul_mz_rewin;
	int mz_sample_dist;
	int bf_shift;
	int max_kmer_cnt;
	double high_factor; // coverage cutoff set to high_factor*hom_cov
	double max_ov_diff_ec;
	double max_ov_diff_final;
	int hom_cov;
    int het_cov;
    int b_low_cov;
    int b_high_cov;
    double m_rate;
	int max_n_chain; // fall-back max number of chains to consider
	int min_hist_kmer_cnt;
    int load_index_from_disk;
    int write_index_to_disk;
    //int number_of_round;
    int adapterLen;
    int clean_round;
    int roundID;
    int max_hang_Len;
    int gap_fuzz;
    int min_overlap_Len;
    int min_overlap_coverage;
    int max_short_tip;
    int max_short_ul_tip;
    int min_cnt;
    int mid_cnt;
    int purge_level_primary;
    int purge_level_trio;
    int purge_overlap_len;
    ///int purge_overlap_len_hic;
    int recover_atg_cov_min;
    int recover_atg_cov_max;
    int hom_global_coverage;
    int hom_global_coverage_set;
    int pur_global_coverage;
    int bed_inconsist_rate;

    float max_hang_rate;
    float min_drop_rate;
    float max_drop_rate;
    float purge_simi_rate_l2;
    float purge_simi_rate_l3;
    float purge_simi_thres;
    float trans_base_rate;
    float trans_base_rate_sec;
    float min_path_drop_rate;
    float max_path_drop_rate;
    // uint64_t path_clean_round;

    ///float purge_simi_rate_hic;

    long long small_pop_bubble_size;
    long long large_pop_bubble_size;
    long long num_bases;
    long long num_corrected_bases;
    long long num_recorrected_bases;
	long long mem_buf;
    long long coverage;
    int hap_occ;
    int polyploidy;
    int trio_flag_occ_thres;
    uint64_t seed;
    int32_t n_perturb;
    double f_perturb;
    int32_t n_weight;
    uint32_t is_alt;
    uint64_t misjoin_len;
    uint64_t scffold;
    int32_t dp_min_len;
    float dp_e;
    int64_t hg_size;
    float kpt_rate;
    int64_t infor_cov, s_hap_cov, trio_cov_het_ovlp;
    double ul_error_rate, ul_error_rate_low, ul_error_rate_hpc;
    int32_t ul_ec_round;
    uint8_t is_dbg_het_cnt;
    uint8_t is_low_het_ul;
    uint8_t is_base_trans;
    uint8_t is_read_trans;
    uint8_t is_topo_trans;
    uint8_t is_bub_trans;
    uint8_t bin_only;
    int32_t ul_clean_round;
    int32_t prt_dbg_gfa;
    int32_t integer_correct_round;
    uint8_t dbg_ovec_cal;
    uint8_t hifi_pst_join, ul_pst_join;
    uint32_t ul_min_base;
    uint8_t self_scaf;
    uint64_t self_scaf_min;
    uint64_t self_scaf_reliable_min;
    int64_t self_scaf_gap_max;
    int64_t somatic_cov;
    std::string output_dir_ec;
}chat_opt_t;

extern chat_opt_t chat_opt;






void init_opt(chat_opt_t *chat_opt);
int check_option(chat_opt_t *chat_opt);
int Dechat_command(int argc, char *argv[],chat_opt_t* chat_opt,hifiasm_opt_t* asm_opt);
#endif