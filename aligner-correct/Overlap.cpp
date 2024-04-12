#include "Overlaps.h"
#include "htab.h"
uint32_t debug_purge_dup = 0;

KDQ_INIT(uint64_t)
KDQ_INIT(uint32_t)

#define ma_hit_key_tn(a) ((a).tn)
KRADIX_SORT_INIT(hit_tn, ma_hit_t, ma_hit_key_tn, member_size(ma_hit_t, tn))
void init_ma_hit_t_alloc(ma_hit_t_alloc* x)
{
    x->size = 0;
    x->buffer = NULL;
    x->length = 0;
    x->is_fully_corrected = 0;
	x->is_abnormal = 0;
}
void ma_hit_sort_tn(ma_hit_t *a, long long n)
{
	radix_sort_hit_tn(a, a + n);
}


void resize_ma_hit_t_alloc(ma_hit_t_alloc* x, uint32_t size)
{
	if (size > x->size) {
		x->size = size;
		kroundup32(x->size);
		REALLOC(x->buffer, x->size);
	}
}

void clear_ma_hit_t_alloc(ma_hit_t_alloc* x)
{
    x->length = 0;
}

void add_ma_hit_t_alloc(ma_hit_t_alloc* x, ma_hit_t* element)
{
	if (x->length + 1 > x->size) {
		x->size = x->length + 1;
		kroundup32(x->size);
		REALLOC(x->buffer, x->size);
	}
	x->buffer[x->length++] = *element;
}