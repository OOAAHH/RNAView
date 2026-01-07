#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "rnaview_profile.h"

RnaViewProfile RNAVIEW_PROFILE = {0};

static void copy_str(char *dst, size_t cap, const char *src)
{
    size_t n;
    if (!dst || cap == 0)
        return;
    if (!src)
    {
        dst[0] = '\0';
        return;
    }
    n = strlen(src);
    if (n >= cap)
        n = cap - 1;
    memcpy(dst, src, n);
    dst[n] = '\0';
}

long long rnaview_profile_now_ns(void)
{
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0)
        return 0;
    return (long long)ts.tv_sec * 1000000000LL + (long long)ts.tv_nsec;
}

int rnaview_profile_is_enabled(void)
{
    return RNAVIEW_PROFILE.enabled != 0;
}

void rnaview_profile_begin(const char *input_path, long num_residue)
{
    const char *path = getenv("RNAVIEW_PROFILE_JSON");
    if (!path || !path[0])
    {
        RNAVIEW_PROFILE.enabled = 0;
        return;
    }

    memset(&RNAVIEW_PROFILE, 0, sizeof(RNAVIEW_PROFILE));
    RNAVIEW_PROFILE.enabled = 1;
    RNAVIEW_PROFILE.num_residue = num_residue;
    copy_str(RNAVIEW_PROFILE.input_path, sizeof(RNAVIEW_PROFILE.input_path), input_path);
    copy_str(RNAVIEW_PROFILE.json_path, sizeof(RNAVIEW_PROFILE.json_path), path);
    RNAVIEW_PROFILE.begin_ns = rnaview_profile_now_ns();
}

static void json_kv_ll(FILE *fp, const char *k, long long v, int comma)
{
    fprintf(fp, "    \"%s\": %lld%s\n", k, v, comma ? "," : "");
}

static void json_kv_long(FILE *fp, const char *k, long v, int comma)
{
    fprintf(fp, "    \"%s\": %ld%s\n", k, v, comma ? "," : "");
}

void rnaview_profile_dump(void)
{
    FILE *fp;

    if (!RNAVIEW_PROFILE.enabled || !RNAVIEW_PROFILE.json_path[0])
        return;

    RNAVIEW_PROFILE.end_ns = rnaview_profile_now_ns();

    fp = fopen(RNAVIEW_PROFILE.json_path, "w");
    if (!fp)
        return;

    fprintf(fp, "{\n");
    fprintf(fp, "  \"schema_version\": 1,\n");
    fprintf(fp, "  \"input\": \"%s\",\n", RNAVIEW_PROFILE.input_path);
    fprintf(fp, "  \"num_residue\": %ld,\n", RNAVIEW_PROFILE.num_residue);

    fprintf(fp, "  \"counts\": {\n");
    json_kv_long(fp, "cand_pairs", RNAVIEW_PROFILE.cand_pairs, 1);
    json_kv_long(fp, "all_pairs_check_pairs_calls", RNAVIEW_PROFILE.all_pairs_check_pairs_calls, 1);
    json_kv_long(fp, "all_pairs_base_stack_calls", RNAVIEW_PROFILE.all_pairs_base_stack_calls, 1);
    json_kv_long(fp, "all_pairs_hbond_pair_calls", RNAVIEW_PROFILE.all_pairs_hbond_pair_calls, 1);
    json_kv_long(fp, "all_pairs_hbond_pair_h_catalog_calls", RNAVIEW_PROFILE.all_pairs_hbond_pair_h_catalog_calls, 1);
    json_kv_long(fp, "all_pairs_lw_pair_type_calls", RNAVIEW_PROFILE.all_pairs_lw_pair_type_calls, 1);
    json_kv_long(fp, "all_pairs_lw_get_hbond_ij_calls", RNAVIEW_PROFILE.all_pairs_lw_get_hbond_ij_calls, 1);
    json_kv_long(fp, "best_pair_check_pairs_calls", RNAVIEW_PROFILE.best_pair_check_pairs_calls, 0);
    fprintf(fp, "  },\n");

    fprintf(fp, "  \"times_ns\": {\n");
    json_kv_ll(fp, "total", RNAVIEW_PROFILE.end_ns - RNAVIEW_PROFILE.begin_ns, 1);
    json_kv_ll(fp, "base_info", RNAVIEW_PROFILE.base_info_ns, 1);
    json_kv_ll(fp, "all_pairs_total", RNAVIEW_PROFILE.all_pairs_total_ns, 1);
    json_kv_ll(fp, "all_pairs_candidate", RNAVIEW_PROFILE.all_pairs_candidate_ns, 1);
    json_kv_ll(fp, "all_pairs_check_pairs", RNAVIEW_PROFILE.all_pairs_check_pairs_ns, 1);
    json_kv_ll(fp, "all_pairs_base_stack", RNAVIEW_PROFILE.all_pairs_base_stack_ns, 1);
    json_kv_ll(fp, "all_pairs_hbond_pair", RNAVIEW_PROFILE.all_pairs_hbond_pair_ns, 1);
    json_kv_ll(fp, "all_pairs_hbond_pair_h_catalog", RNAVIEW_PROFILE.all_pairs_hbond_pair_h_catalog_ns, 1);
    json_kv_ll(fp, "all_pairs_lw_pair_type", RNAVIEW_PROFILE.all_pairs_lw_pair_type_ns, 1);
    json_kv_ll(fp, "all_pairs_lw_get_hbond_ij", RNAVIEW_PROFILE.all_pairs_lw_get_hbond_ij_ns, 1);
    json_kv_ll(fp, "best_pair_total", RNAVIEW_PROFILE.best_pair_total_ns, 1);
    json_kv_ll(fp, "best_pair_check_pairs", RNAVIEW_PROFILE.best_pair_check_pairs_ns, 0);
    fprintf(fp, "  }\n");

    fprintf(fp, "}\n");
    fclose(fp);
}

void rnaview_profile_add_all_pairs_hbond_pair_h_catalog(long long delta_ns)
{
    if (!RNAVIEW_PROFILE.enabled)
        return;
    RNAVIEW_PROFILE.all_pairs_hbond_pair_h_catalog_ns += delta_ns;
    RNAVIEW_PROFILE.all_pairs_hbond_pair_h_catalog_calls++;
}

void rnaview_profile_add_all_pairs_lw_get_hbond_ij(long long delta_ns)
{
    if (!RNAVIEW_PROFILE.enabled)
        return;
    RNAVIEW_PROFILE.all_pairs_lw_get_hbond_ij_ns += delta_ns;
    RNAVIEW_PROFILE.all_pairs_lw_get_hbond_ij_calls++;
}
