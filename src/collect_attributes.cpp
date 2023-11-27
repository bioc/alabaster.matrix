#include "Rcpp.h"
#include <cmath>
#include <algorithm>

//[[Rcpp::export(rng=false)]]
double lowest_double() {
    return std::numeric_limits<double>::lowest();
}

//[[Rcpp::export(rng=false)]]
double highest_double() {
    return std::numeric_limits<double>::max();
}

//[[Rcpp::export(rng=false)]]
Rcpp::List collect_double_attributes(Rcpp::NumericVector x) {
    bool has_missing = false;
    for (auto y : x) {
        if (ISNA(y)) {
            has_missing = true;
            break;
        }
    }

    bool has_nan = false;
    bool has_posinf = false;
    bool has_neginf = false;
    {
        for (auto y : x) {
            if (!ISNA(y) && std::isnan(y)) {
                has_nan = true;
                break;
            }
        }
        for (auto y : x) {
            if (!ISNA(y) && std::isinf(y) && y > 0) {
                has_posinf = true;
                break;
            }
        }
        for (auto y : x) {
            if (!ISNA(y) && std::isinf(y) && y < 0) {
                has_neginf = true;
                break;
            }
        }
    }

    bool non_integer = true;
    double minv = R_PosInf, maxv = R_NegInf;
    if (!has_nan && !has_posinf && !has_neginf) {
        non_integer = false;
        for (auto y : x) {
            if (!ISNA(y)) {
                if (std::floor(y) != y) {
                    non_integer = true;
                    break;
                }
            }
        }

        if (!non_integer) {
            for (auto y : x) {
                if (!ISNA(y)) {
                    minv = std::min(y, minv);
                    maxv = std::max(y, maxv);
                }
            }
        }
    }

    bool has_lowest = false, has_highest = false;
    {
        double lowest = lowest_double();
        for (auto y : x) {
            if (!ISNA(y) && y == lowest) {
                has_lowest = true;
                break;
            }
        }

        double highest = highest_double();
        for (auto y : x) {
            if (!ISNA(y) && y == highest) {
                has_highest = true;
                break;
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("range") = Rcpp::NumericVector::create(minv, maxv),
        Rcpp::Named("missing") = Rcpp::LogicalVector::create(has_missing),
        Rcpp::Named("non_integer") = Rcpp::LogicalVector::create(non_integer),
        Rcpp::Named("has_NaN") = Rcpp::LogicalVector::create(has_nan),
        Rcpp::Named("has_Inf") = Rcpp::LogicalVector::create(has_posinf),
        Rcpp::Named("has_nInf") = Rcpp::LogicalVector::create(has_neginf),
        Rcpp::Named("has_lowest") = Rcpp::LogicalVector::create(has_lowest),
        Rcpp::Named("has_highest") = Rcpp::LogicalVector::create(has_highest)
    );
}
