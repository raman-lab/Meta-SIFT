##  Cython  ##
# Author: Kristopher Kieft, UW-Madison
# Feb 2020 (updated Dec 2021)

# methods for parent motif finding script

cpdef min_max_size(str seq_aa, int Xs, int min_size, int max_size):
    cdef int l_aa, l_aa_Xs

    l_aa = len(seq_aa)
    l_aa_Xs = l_aa-Xs
    if l_aa_Xs >= min_size and l_aa_Xs <= max_size:
        return True,l_aa,l_aa_Xs
    return False,l_aa,l_aa_Xs


cpdef filter_scores(list score_check, int length, long double cutoff):
    cdef long double sc_zero, sc_minus_one, sc_minus_two, combo, aa_score

    sc_zero = score_check[0]
    sc_minus_one = score_check[-1]
    sc_minus_two = score_check[-2]
    combo = sc_minus_one + sc_minus_two

    aa_score = sum(score_check)
    if aa_score/length >= cutoff and sc_zero >= 0.75 and combo/aa_score < 0.8:
        return True,aa_score
    return False,aa_score

 
cpdef get_score_max(long double score_max, long double m, int length):
    score_max += m/length
    return score_max


cpdef get_zeros(list md_values):
    cdef int zeros
    zeros = len([val for val in md_values if val > 0])
    return zeros


cpdef get_final_scores(long double score_max, int l_aa_Xs, int windows, int win_check):
    cdef long double final, e_val

    final = (score_max/(l_aa_Xs))
    e_val = windows*win_check*(2**(-1*score_max))
    return final,e_val


cpdef evalue_final(int windows, long double final, long double threshold, long double e_val, long double evalue):
    if windows > 0 and final >= threshold and e_val <= evalue:
        return True
    return False


cpdef mot_dict_math(int p, int sum_p):
    cdef long double new 

    new = p/sum_p
    return new

cpdef score_check_len(int l_sc, int length):
 
    if l_sc == length:
        return True
    return False

cpdef win_check_add(int win_check):
    win_check += 1
    return win_check

cpdef windows_add(int windows):
    windows += 1
    return windows