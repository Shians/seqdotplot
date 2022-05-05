#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

using std::string;
using std::string_view;
using std::vector;

inline const size_t matched_chars(string_view s1, string_view s2) {
    size_t n_max = s1.length();
    size_t n_match = 0;

    for (auto i = 0; i < n_max; i++) {
        if (s1[i] == s2[i]) {
            n_match++;
        }
    }

    return(n_match);
}

// [[Rcpp::export]]
DataFrame compute_dotplot_data(String r_s1, String r_s2, int wsize, int step, int n_mismatches) {
    string const s1(r_s1.get_cstring());
    string const s2(r_s2.get_cstring());

    auto const s1_len = s1.length();
    auto const s2_len = s2.length();

    int const ncol = s1_len - wsize + 1;
    int const nrow = s2_len - wsize + 1;

    vector<int> x_coord;
    vector<int> y_coord;

    string_view s1_v(s1);
    string_view s2_v(s2);

    for (auto i = 0; i < ncol; i += step) {
        auto sub_v1 = s1_v.substr(i, wsize);
        R_CheckUserInterrupt();
        for (auto j = 0; j < nrow; j += step) {
            auto sub_v2 = s2_v.substr(j, wsize);
            if (wsize - matched_chars(sub_v1, sub_v2) <= n_mismatches) {
                x_coord.push_back(i);
                y_coord.push_back(j);
            }
        }
    }

    return DataFrame::create(
        Named("x") = x_coord,
        Named("y") = y_coord
    );
}

