#ifndef RNAVIEW_PROFILE_H
#define RNAVIEW_PROFILE_H

typedef struct RnaViewProfile
{
    int enabled;
    long num_residue;

    long cand_pairs;
    long all_pairs_check_pairs_calls;
    long all_pairs_base_stack_calls;
    long all_pairs_hbond_pair_calls;
    long all_pairs_hbond_pair_h_catalog_calls;
    long all_pairs_lw_pair_type_calls;
    long all_pairs_lw_get_hbond_ij_calls;
    long best_pair_check_pairs_calls;

    long long begin_ns;
    long long end_ns;

    long long base_info_ns;
    long long all_pairs_total_ns;
    long long all_pairs_candidate_ns;
    long long all_pairs_check_pairs_ns;
    long long all_pairs_base_stack_ns;
    long long all_pairs_hbond_pair_ns;
    long long all_pairs_hbond_pair_h_catalog_ns;
    long long all_pairs_lw_pair_type_ns;
    long long all_pairs_lw_get_hbond_ij_ns;
    long long best_pair_total_ns;
    long long best_pair_check_pairs_ns;

    char input_path[1024];
    char json_path[1024];
} RnaViewProfile;

extern RnaViewProfile RNAVIEW_PROFILE;

int rnaview_profile_is_enabled(void);
long long rnaview_profile_now_ns(void);
void rnaview_profile_begin(const char *input_path, long num_residue);
void rnaview_profile_add_all_pairs_hbond_pair_h_catalog(long long delta_ns);
void rnaview_profile_add_all_pairs_lw_get_hbond_ij(long long delta_ns);
void rnaview_profile_dump(void);

#endif /* RNAVIEW_PROFILE_H */
