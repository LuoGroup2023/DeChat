#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include "Process_Read.h"
#include "htab.h"
#include "Correct.h"
#include "kalloc.h"
#include <assert.h>


#define UL_FLANK 512
uint8_t seq_nt6_table[256] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 0, 5, 1,  5, 5, 5, 2,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  3, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 0, 5, 1,  5, 5, 5, 2,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  3, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

char bit_t_seq_table[256][4] = {{0}};
char bit_t_seq_table_rc[256][4] = {{0}};
char s_H[5] = {'A', 'C', 'G', 'T', 'N'};
char rc_Table[6] = {'T', 'G', 'C', 'A', 'N', 'N'};


void init_All_reads(All_reads* r)
{
	memset(r, 0, sizeof(All_reads));
	r->index_size = READ_INIT_NUMBER;
	r->read_length = (uint64_t*)malloc(sizeof(uint64_t)*r->index_size);
	r->name_index_size = READ_INIT_NUMBER;
	r->name_index = (uint64_t*)malloc(sizeof(uint64_t)*r->name_index_size);
	r->name_index[0] = 0;
}

void recover_UC_Read_RC(UC_Read* r, All_reads* R_INF, uint64_t ID)
{
	r->length = Get_READ_LENGTH((*R_INF), ID);
	uint8_t* src = Get_READ((*R_INF), ID);

	if (r->length + 4 > r->size) {
		r->size = r->length + 4;
		r->seq = (char*)realloc(r->seq,sizeof(char)*(r->size));
	}

	long long last_chr = r->length % 4;
	long long i = r->length / 4 - 1 + (last_chr != 0);
	long long index = 0;

	if(src) {
		if(last_chr!=0) {
			memcpy(r->seq + index, bit_t_seq_table_rc[src[i]] + 4 - last_chr, last_chr);
			index = last_chr;
			i--;
		}

		while (i >= 0) {
			memcpy(r->seq + index, bit_t_seq_table_rc[src[i]], 4);
			i--;
			index = index + 4;
		}

		if (R_INF->N_site[ID]) {
			for (i = 1; i <= (long long)R_INF->N_site[ID][0]; i++) {
				r->seq[r->length - R_INF->N_site[ID][i] - 1] = 'N';
			}
		}
	} else {///N
		memset(r->seq, 'N', r->length);
	}

	
}

void reverse_complement(char* pattern, uint64_t length)
{
	uint64_t i = 0;
	uint64_t end = length / 2;
	char k;
	uint64_t index;

	for (i = 0; i < end; i++)
	{
		
		index = length - i - 1;
		k = pattern[index];
		pattern[index] = RC_CHAR(pattern[i]);
		pattern[i] = RC_CHAR(k);
	}

	if(length&(uint64_t)1)
	{
		pattern[end] = RC_CHAR(pattern[end]);
	}
}


void recover_UC_Read(UC_Read* r, const All_reads *R_INF, uint64_t ID)
{
	r->length = Get_READ_LENGTH((*R_INF), ID);
	uint8_t* src = Get_READ((*R_INF), ID);

	if (r->length + 4 > r->size) {
		r->size = r->length + 4;
		r->seq = (char*)realloc(r->seq,sizeof(char)*(r->size));
	}

	uint64_t i = 0;

	if(src) {
		while ((long long)i < r->length) {
			memcpy(r->seq+i, bit_t_seq_table[src[i>>2]], 4);
			i = i + 4;
		}

		if (R_INF->N_site[ID]) {
			for (i = 1; i <= R_INF->N_site[ID][0]; i++) r->seq[R_INF->N_site[ID][i]] = 'N';
		}
	} else {///N
		memset(r->seq, 'N', r->length);
	}

	r->RID = ID;	
}


void malloc_All_reads(All_reads* r)
{
	r->read_size = (uint64_t*)malloc(sizeof(uint64_t)*r->total_reads);
	memcpy(r->read_size, r->read_length, sizeof(uint64_t)*r->total_reads);

	r->read_sperate = (uint8_t**)malloc(sizeof(uint8_t*)*r->total_reads);
	long long i = 0;
	for (i = 0; i < (long long)r->total_reads; i++)
	{
		r->read_sperate[i] = (uint8_t*)malloc(sizeof(uint8_t)*(r->read_length[i]/4+1));
	}

	r->cigars = (Compressed_Cigar_record*)malloc(sizeof(Compressed_Cigar_record)*r->total_reads);
	r->second_round_cigar = (Compressed_Cigar_record*)malloc(sizeof(Compressed_Cigar_record)*r->total_reads);
	r->paf = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*r->total_reads);
	r->reverse_paf = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*r->total_reads);
	///r->pb_regions = (kvec_t_u64_warp*)malloc(r->total_reads*sizeof(kvec_t_u64_warp));

	for (i = 0; i < (long long)r->total_reads; i++)
	{
		r->second_round_cigar[i].size = r->cigars[i].size = 0;
		r->second_round_cigar[i].length = r->cigars[i].length = 0;
		r->second_round_cigar[i].record = r->cigars[i].record = NULL;

		r->second_round_cigar[i].lost_base_size = r->cigars[i].lost_base_size = 0;
		r->second_round_cigar[i].lost_base_length = r->cigars[i].lost_base_length = 0;
		r->second_round_cigar[i].lost_base = r->cigars[i].lost_base = NULL;
		init_ma_hit_t_alloc(&(r->paf[i]));
		init_ma_hit_t_alloc(&(r->reverse_paf[i]));
		///kv_init(r->pb_regions[i].a);
	}

	r->name = (char*)malloc(sizeof(char)*r->total_name_length);
	r->N_site = (uint64_t**)calloc(r->total_reads, sizeof(uint64_t*));
	r->trio_flag = (uint8_t*)malloc(r->total_reads*sizeof(uint8_t));
	memset(r->trio_flag, AMBIGU, r->total_reads*sizeof(uint8_t));
}

void init_UC_Read(UC_Read* r)
{
	r->length = 0;
	r->size = 0;
	r->seq = NULL;
	if (bit_t_seq_table[0][0] == 0)
	{
		uint64_t i = 0;

		for (i = 0; i < 256; i++)
		{
			bit_t_seq_table[i][0] = s_H[((i >> 6)&(uint64_t)3)];
			bit_t_seq_table[i][1] = s_H[((i >> 4)&(uint64_t)3)];
			bit_t_seq_table[i][2] = s_H[((i >> 2)&(uint64_t)3)];
			bit_t_seq_table[i][3] = s_H[(i&(uint64_t)3)];

			bit_t_seq_table_rc[i][0] = RC_CHAR(bit_t_seq_table[i][3]);
			bit_t_seq_table_rc[i][1] = RC_CHAR(bit_t_seq_table[i][2]);
			bit_t_seq_table_rc[i][2] = RC_CHAR(bit_t_seq_table[i][1]);
			bit_t_seq_table_rc[i][3] = RC_CHAR(bit_t_seq_table[i][0]);
		}
	}
}
void destory_UC_Read(UC_Read* r)
{
	free(r->seq);
}



#define COMPRESS_BASE {c = seq_nt6_table[(uint8_t)src[i]];\
		if (c >= 4)\
		{\
			c = 0;\
			(*N_site_lis)[N_site_i] = i;\
			N_site_i++;\
		}\
		i++;}\

void ha_compress_base(uint8_t* dest, char* src, uint64_t src_l, uint64_t** N_site_lis, uint64_t N_site_occ)
{
	///N_site_lis saves the pos of all Ns in this read
	///N_site_lis[0] is the number of Ns
	free((*N_site_lis));
	if (N_site_occ)
	{
		(*N_site_lis) = (uint64_t*)malloc(sizeof(uint64_t)*(N_site_occ + 1));
		(*N_site_lis)[0] = N_site_occ;
	}
	else
	{
		(*N_site_lis) = NULL;
	}

	uint64_t i = 0;
	uint64_t N_site_i = 1;
	uint64_t dest_i = 0;
	uint8_t tmp = 0;
	uint8_t c = 0;

	while (i + 4 <= src_l)
	{
		tmp = 0;

		COMPRESS_BASE;
		tmp = tmp | (c<<6);

		COMPRESS_BASE;
		tmp = tmp | (c<<4);

		COMPRESS_BASE;
		tmp = tmp | (c<<2);

		COMPRESS_BASE;
		tmp = tmp | c;

		dest[dest_i] = tmp;

		dest_i++;
	}

	//at most 3 bases here
	uint64_t shift = 6;
	if (i < src_l)
	{
		tmp = 0;

		while (i < src_l)
		{
			COMPRESS_BASE;
			tmp = tmp | (c << shift);
			shift = shift -2;
		}
		
		dest[dest_i] = tmp;
		dest_i++;
	}
}



void ha_insert_read_len(All_reads *r, int read_len, int name_len)
{
	r->total_reads++;
	r->total_reads_bases += (uint64_t)read_len;
	r->total_name_length += (uint64_t)name_len;

	// must +1
	if (r->index_size < r->total_reads + 2) {
		r->index_size = r->index_size * 2 + 2;
		r->read_length = (uint64_t*)realloc(r->read_length, sizeof(uint64_t) * r->index_size);
		r->name_index_size = r->name_index_size * 2 + 2;
		r->name_index = (uint64_t*)realloc(r->name_index, sizeof(uint64_t) * r->name_index_size);
	}

	r->read_length[r->total_reads - 1] = read_len;
	r->name_index[r->total_reads] = r->name_index[r->total_reads - 1] + name_len;
}


//从原始DNA序列中提取子区间
void recover_UC_Read_sub_region(char* r, int64_t start_pos, int64_t length, uint8_t strand, All_reads* R_INF, int64_t ID)
{
	int64_t readLen = Get_READ_LENGTH((*R_INF), ID);
	int64_t end_pos = start_pos + length - 1, begLen, tailLen, offset, mn, src_i, des_i, i;
	uint8_t* src = Get_READ((*R_INF), ID);

	if(strand == 0) {
		offset = start_pos&3;
		begLen = 4-offset;
		if(begLen > length) begLen = length;
		tailLen = (length-begLen)&3;
		mn = (length - begLen - tailLen)>>2;
		src_i = start_pos; des_i = 0; i = 0;

		if(begLen > 0) {
			memcpy(r+des_i, bit_t_seq_table[src[src_i>>2]]+offset, begLen);
			des_i += begLen; src_i += begLen;
		}

		for (i = 0; i < mn; i++) {
			memcpy(r+des_i, bit_t_seq_table[src[src_i>>2]], 4);
			des_i += 4; src_i += 4;
		}

		if(tailLen > 0) {
			memcpy(r+des_i, bit_t_seq_table[src[src_i>>2]], tailLen);
			des_i += tailLen; src_i += tailLen;
		}

		if (R_INF->N_site[ID]) {
			for (i = 1; i <= (long long)R_INF->N_site[ID][0]; i++) {
				if ((long long)R_INF->N_site[ID][i] >= start_pos && (long long)R_INF->N_site[ID][i] <= end_pos) {
					r[R_INF->N_site[ID][i] - start_pos] = 'N';
				}
				else if((long long)R_INF->N_site[ID][i] > end_pos) {
					break;
				}
			}
		}
	}
	else {
		start_pos = readLen - start_pos - 1;
		end_pos = readLen - end_pos - 1;

		begLen = (start_pos+1)&3;
		offset = 4 - begLen;
		if(begLen > length) begLen = length;
		tailLen = (length-begLen)&3;
		mn = (length - begLen - tailLen)>>2;
		src_i = start_pos; des_i = 0; i = 0;

		if(begLen > 0) {
			memcpy(r+des_i, bit_t_seq_table_rc[src[src_i>>2]]+offset, begLen);
			des_i += begLen; src_i -= begLen;
		}

		for (i = 0; i < mn; i++) {
			memcpy(r+des_i, bit_t_seq_table_rc[src[src_i>>2]], 4);
			des_i += 4; src_i -= 4;
		}

		if(tailLen > 0) {
			memcpy(r+des_i, bit_t_seq_table_rc[src[src_i>>2]], tailLen);
			des_i += tailLen; src_i -= tailLen;
		}

		/**
		if (R_INF->N_site[ID]) {
			offset = readLen - start_pos - 1;
			for (i = 1; i <= (long long)R_INF->N_site[ID][0]; i++) {
				if ((long long)R_INF->N_site[ID][i] >= end_pos && (long long)R_INF->N_site[ID][i] <= start_pos) {
					r[readLen - R_INF->N_site[ID][i] - 1 - offset] = 'N';
				}
				else if((long long)R_INF->N_site[ID][i] > start_pos) {
					break;
				}
			}
		}
		**/
		if (R_INF->N_site[ID]) {
			start_pos = readLen - start_pos - 1; end_pos = readLen - end_pos - 1;
			for (i = 1; i <= (long long)R_INF->N_site[ID][0]; i++) {
				offset = readLen - R_INF->N_site[ID][i] - 1;
				if(offset >= start_pos && offset <= end_pos) {
					r[offset - start_pos] = 'N';
				} else if(offset < start_pos) {
					break;
				}
			}
		}
	}
}

uint32_t retrieve_u_cov(const ul_idx_t *ul, uint64_t id, uint8_t strand, uint64_t pos, uint8_t dir, int64_t *pi)
{
	uint64_t *a = ul->cc->interval.a + ul->cc->idx[id], cc = 0, ff = 0;
	int64_t a_n = ul->cc->idx[id+1]-ul->cc->idx[id], k = 0, cc_i = pi? *pi:0;
	if(a_n == 0) return 0;
	if(cc_i + 1 >= a_n || cc_i < 0) cc_i = 0;
	if(strand) pos = ul->ug->u.a[id].len - pos - 1;
	if(dir == 0) {
		for (k = cc_i; k + 1 < a_n; k++) {
			if(pos>=(a[k]>>32) && pos<(a[k+1]>>32)) {
				cc = (uint32_t)a[k];
				ff = 1;
				break;
			}
		}

		if(ff == 0) {
			for (k = 0; k < cc_i; k++) {
					if(pos>=(a[k]>>32) && pos<(a[k+1]>>32)) {
					cc = (uint32_t)a[k];
					ff = 1;
					break;
				}
			}
		}
	} else {
		for (k = cc_i; k >= 0; k--) {
			if(pos>=(a[k]>>32) && pos<(a[k+1]>>32)) {
				cc = (uint32_t)a[k];
				ff = 1;
				break;
			}
		}

		if(ff == 0) {
			for (k = cc_i+1; k + 1 < a_n; k++) {
				if(pos>=(a[k]>>32) && pos<(a[k+1]>>32)) {
					cc = (uint32_t)a[k];
					ff = 1;
					break;
				}
			}
		}
	}
	
	if(pi) *pi = ff?k:0;
	return cc;
}

void retrieve_u_seq(UC_Read* i_r, char* i_s, ma_utg_t *u, uint8_t strand, int64_t s, int64_t l, void *km)
{
	if(u->m == 0 || u->n == 0) return;
	if(l < 0) l = u->len;
	char *r = NULL, *a = NULL;
	int64_t e = s + l, ssp, sep, rs, re, des_i;
	uint64_t k, rId, ori, r_l;
	if(i_r) {
        i_r->length = l; i_r->RID = 0;
        if(i_r->length > i_r->size) {
            i_r->size = i_r->length;
			if(!km) REALLOC(i_r->seq, i_r->size);
			else KREALLOC(km, i_r->seq, i_r->size);
            // i_r->seq = (char*)realloc(i_r->seq,sizeof(char)*(i_r->size));
        }
        r = i_r->seq;
    }
    if(i_s) r = i_s;

	if(strand == 1) {
		sep = u->len - s; 
        ssp = u->len - e;
        s = ssp; e = sep; 
	}
	for (k = l = des_i = 0; k < u->n; k++) {
		rId = u->a[k]>>33;
		ori = u->a[k]>>32&1;
		r_l = (uint32_t)u->a[k];
		if(r_l == 0) continue;
		ssp = l; sep = l + r_l;
		l += r_l;
		if(sep <= s) continue;
		if(ssp >= e) break;
		rs = MAX(ssp, s); re = MIN(sep, e); 
		a = r + des_i; des_i += re - rs;
		recover_UC_Read_sub_region(a, rs-ssp, re-rs, ori, &R_INF, rId);
	}
	if(strand == 1) {
		char t;
		re = (e - s);
		l = re>>1; 
		for (k = 0; k < (uint64_t)l; k++) {
			des_i = re - k - 1;
			t = r[des_i];
			r[des_i] = RC_CHAR(r[k]);
			r[k] = RC_CHAR(t);
		}
		if(re&1) r[l] = RC_CHAR(r[l]);
	}
}


void write_All_reads(All_reads* r, char* read_file_name)
{
    fprintf(stderr, "Writing reads to disk... \n");
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "w");
	fwrite(&asm_opt.adapterLen, sizeof(asm_opt.adapterLen), 1, fp);
    fwrite(&r->index_size, sizeof(r->index_size), 1, fp);
	fwrite(&r->name_index_size, sizeof(r->name_index_size), 1, fp);
	fwrite(&r->total_reads, sizeof(r->total_reads), 1, fp);
	fwrite(&r->total_reads_bases, sizeof(r->total_reads_bases), 1, fp);
	fwrite(&r->total_name_length, sizeof(r->total_name_length), 1, fp);

	uint64_t i = 0;
	uint64_t zero = 0;
	for (i = 0; i < r->total_reads; i++)
	{
		if (r->N_site[i] != NULL)
		{
			///number of Ns
			fwrite(&r->N_site[i][0], sizeof(r->N_site[i][0]), 1, fp);
			if (r->N_site[i][0])
			{
				fwrite(r->N_site[i]+1, sizeof(r->N_site[i][0]), r->N_site[i][0], fp);
			}
		}
		else
		{
			fwrite(&zero, sizeof(zero), 1, fp);
		}
	}

	fwrite(r->read_length, sizeof(uint64_t), r->total_reads, fp);
	for (i = 0; i < r->total_reads; i++)
	{
		fwrite(r->read_sperate[i], sizeof(uint8_t), r->read_length[i]/4+1, fp);
	}
	
	fwrite(r->name, sizeof(char), r->total_name_length, fp);
	fwrite(r->name_index, sizeof(uint64_t), r->name_index_size, fp);
	fwrite(r->trio_flag, sizeof(uint8_t), r->total_reads, fp);
	fwrite(&(asm_opt.hom_cov), sizeof(asm_opt.hom_cov), 1, fp);
	fwrite(&(asm_opt.het_cov), sizeof(asm_opt.het_cov), 1, fp);

    free(index_name);    
	fflush(fp);
    fclose(fp);
    fprintf(stderr, "Reads has been written.\n");
}