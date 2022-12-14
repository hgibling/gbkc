// count - Score reads based on comparison of k-mer count profiles to those for known alleles

#include <algorithm>
#include <cmath>
#include <getopt.h>
#include <iostream>
#include <map>
#include <math.h>
#include <numeric>
#include <omp.h>
#include <set>
#include <sstream>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include <zlib.h>

#include "check-profiles.h"
#include "count.h"
#include "kseq.h"

using std::cout;
using std::greater;
using std::istringstream;
using std::pair;
using std::map;
using std::set;
using std::sort;
using std::string;
using std::unordered_set;
using std::vector;


//
// Define functions
//

// Calculate lambda
double calculate_lambda(const double read_length, const size_t k, const double coverage, const double sequencing_error)
{
    double lambda = (read_length - k + 1) * (coverage / read_length) * pow((1 - sequencing_error), k);
    return lambda;
}

//// Original Cortex formula
double calculate_lambda_cortex(const double read_length, const size_t k, const double coverage, const double sequencing_error)
{
    double lambda = (read_length - k + 1) * (coverage / read_length) * (1 - (sequencing_error * k));
    return lambda;
}

// Log factorial (from Jared Simpson)
double log_factorial(size_t c)
{
    double result = 0;
    while(c > 0) {
        result += log(c--); //slow
    }
    return result;
}

// Poisson distribution probability mass function (from Jared Simpson)
double log_poisson_pmf(const double c, const double lambda)
{
    double f_c = log_factorial(c);
    double p = (double)c * log(lambda) - lambda - f_c;
    return p;
}

// Score reads k-mer count against allele k-mer count (individual k-mers)
double score_kmer(const size_t read_count, const size_t allele_count, const double lambda, const double lambda_error)
{
    double score = 0;
    if (allele_count == 0) {
        score = log_poisson_pmf(read_count, lambda_error);
    } else {
        score = log_poisson_pmf(read_count, (lambda * allele_count));
    }
    return score;
}

// Score read k-mer count proflile against allele k-mer count profile
double score_profile(const kmer_count_map& read_map, const kmer_count_map& allele_map, const unordered_set<string>& allele_union_kmers, const double lambda, const double lambda_error)
{
    double score = 0;
    double single_score = 0;
    for (auto iter = allele_union_kmers.begin(); iter != allele_union_kmers.end(); ++iter) {
        string kmer = *iter;
        auto read_iter = read_map.find(kmer);
        size_t kmer_count_in_read = read_iter != read_map.end() ? read_iter->second : 0;
        auto allele_iter = allele_map.find(kmer);
        size_t kmer_count_in_allele = allele_iter != allele_map.end() ? allele_iter->second : 0;
        single_score = score_kmer(kmer_count_in_read, kmer_count_in_allele, lambda, lambda_error);
        score += single_score;
    }
    return score;
}

// Calculate Pearson correlation coefficient
double calculate_pearson(const vector<size_t>& x_vector, const vector<size_t>& y_vector) {
    double x_sum = 0;
    double y_sum = 0;
    for (size_t i = 0; i < x_vector.size(); ++i) {
        x_sum += x_vector[i];
        y_sum += y_vector[i];
    }

    double x_mean = x_sum / x_vector.size();
    double y_mean = y_sum / y_vector.size();
    double sum_minus_mean_multiplied = 0;
    double sum_sq_minus_mean_x = 0;
    double sum_sq_minus_mean_y = 0;

    for (size_t i = 0; i < x_vector.size(); ++i) {
        double x_minus_mean = x_vector[i] - x_mean;
        double y_minus_mean = y_vector[i] - y_mean;
        sum_minus_mean_multiplied += x_minus_mean * y_minus_mean;
        sum_sq_minus_mean_x += x_minus_mean * x_minus_mean;
        sum_sq_minus_mean_y += y_minus_mean * y_minus_mean;
    }

    double r = sum_minus_mean_multiplied / (sqrt(sum_sq_minus_mean_x * sum_sq_minus_mean_y));
    return(r);
}

// Score with Pearson correlation coefficient
double score_pearson(const kmer_count_map& read_map, const kmer_count_map& allele_map, const unordered_set<string>& allele_union_kmers) {
    vector<size_t> read_count_vector;
    vector<size_t> allele_count_vector;

    for (auto iter = allele_union_kmers.begin(); iter != allele_union_kmers.end(); ++iter) {
        string kmer = *iter;
        auto read_iter = read_map.find(kmer);
        size_t kmer_count_in_read = read_iter != read_map.end() ? read_iter->second : 0;
        read_count_vector.push_back(kmer_count_in_read);
        auto allele_iter = allele_map.find(kmer);
        size_t kmer_count_in_allele = allele_iter != allele_map.end() ? allele_iter->second : 0;
        allele_count_vector.push_back(kmer_count_in_allele);
    }

    double pearson = calculate_pearson(read_count_vector, allele_count_vector);
    return(pearson);
}



//
// Help message
//

static const char *COUNT_USAGE_MESSAGE =
"Score reads based on comparison of k-mer count profiles to those for known alleles.\n\n"
"Usage: gbkc count [options]\n\n"
"Commands:\n"
"       -a       multi-fasta file of alleles/haplotypes of interest\n"
"       -1       multi-fasta/q file of sequencing reads to score (first in pair)\n"
"       -2       [optional] multi-fasta/q file of sequencing reads to score (second in pair)\n"
"       -d       score diploid genotypes (if flag is not used, haploid scoring is performed)\n"
"       -k       lower value of k-mer to use (default: 11)\n"
"       -K       upper value of k-mer to use\n"
"       -i       increment size between k-mer values (default: 4)\n"
"       -l       read length\n"
"       -e       sequencing error rate\n"
"       -c       sequencing coverage\n"
"       -L       error rate for lambda (default: 1)\n"
"       -m       method for calculating lambda (from coverage and error rate, or from flank k-mer counts; one of: coverage, mean, median)\n"
"       -M       manual lambda entry (overrides lambda method selection)\n"
"       -f       multi-fasta file of the two flanking sequences surrounding region of interest\n"
"       -N       print only the top N scores per k-mer (default: print all)\n"
"       -o       output file name (default: results.csv)\n"
"       -t       number of threads (default: 1)\n";


//
// Program commands
//

int countMain(int argc, char** argv) {

    if(argc <= 1) {
        cout << COUNT_USAGE_MESSAGE;
        return 0;
    };


    //
    // Read command line arguments
    //

    string input_alleles_file;
    string input_reads_file1;
    string input_reads_file2;
    bool is_diploid = false;
    size_t lower_k = 11;
    size_t upper_k = 0;
    size_t increment_k = 4;
    double read_length = -1;
    double sequencing_error = -1;
    double coverage = -1;
    double lambda_error = 1;
    string lambda_method;
    double manual_lambda = -1;       // temporary for troubleshooting
    string input_flanks_file;
    int top_N = -1;
    string output_name = "count-results.csv";
    size_t num_threads = 1;

    for (char c; (c = getopt_long(argc, argv, "a:1:2:dk:K:i:l:e:c:L:m:M:f:N:o:t:", NULL, NULL)) != -1;) {
        istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'a': arg >> input_alleles_file; break;
            case '1': arg >> input_reads_file1; break;
            case '2': arg >> input_reads_file2; break;
            case 'd': is_diploid = true; break;
            case 'k': arg >> lower_k; break;
            case 'K': arg >> upper_k; break;
            case 'i': arg >> increment_k; break;
            case 'l': arg >> read_length; break;        // TODO: calculate from input reads files
            case 'e': arg >> sequencing_error; break;
            case 'c': arg >> coverage; break;
            case 'L': arg >> lambda_error; break;
            case 'm': arg >> lambda_method; break;
            case 'M': arg >> manual_lambda; break;     // temporary for troubleshooting
            case 'f': arg >> input_flanks_file; break;
            case 'N': arg >> top_N; break;
            case 'o': arg >> output_name; break;
            case 't': arg >> num_threads; break;
            default: exit(EXIT_FAILURE);
        }
    }



    // Check for parameter arguments
    if (input_alleles_file.empty()) {
        fprintf(stderr, "No file for allele sequences. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (input_reads_file1.empty()) {
        fprintf(stderr, "No file for read sequences. One file must be specified for argument -1. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (read_length <= 0) {
        fprintf(stderr, "Read length must be greater than 0. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (lower_k == 0 || lower_k > read_length || upper_k == 0 || upper_k > read_length) {
        fprintf(stderr, "Value of k must be non-zero and less than the read length. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (lower_k > upper_k || increment_k < 0 || increment_k > upper_k) {
        fprintf(stderr, "Lower k value and k increment size must not be greater than upper k value. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (sequencing_error < 0 || sequencing_error > 1) {
        fprintf(stderr, "Sequencing error rate must be between 0 and 1. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (coverage <= 0) {
        fprintf(stderr, "Coverage must be greater than 0. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (lambda_method != "coverage" && lambda_method != "mean" && lambda_method != "median") {
        fprintf(stderr, "Method must be one of 'coverage', 'mean', or 'median'. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (input_flanks_file.empty()) {
        fprintf(stderr, "No file for flank sequences. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (num_threads <= 0) {
        fprintf(stderr, "Number of threads must be greater than 0. Check parameters.\n");
        exit(EXIT_FAILURE);
    }


    //
    // Read files
    //

    // TODO: Handle case where file names are null

    // Get allele sequences
    vector<sequence_record> alleles = read_sequences_from_file(input_alleles_file);

    // Get list of allele names
    set<string> allele_names;
    for (size_t i = 0; i < alleles.size(); ++i) {
        allele_names.insert(alleles[i].name.c_str());
    }

    // Get read sequences
    vector<sequence_record> reads1 = read_sequences_from_file(input_reads_file1);
    vector<sequence_record> reads2;
    if (!input_reads_file2.empty()) {
        reads2 = read_sequences_from_file(input_reads_file2);
    }

    // Get flank sequences
    vector<sequence_record> flanks = read_sequences_from_file(input_flanks_file);


    //
    // Print handy information
    //

    fprintf(stderr, "Number of threads used: %zu\n", num_threads);
    fprintf(stderr, "Input reads: %s", input_reads_file1.c_str());
    fprintf(stderr, " %s", input_reads_file2.c_str());
    fprintf(stderr, "\nInput alleles: %s\n", input_alleles_file.c_str());
    fprintf(stderr, "Input flank sequences: %s\n", input_flanks_file.c_str());
    fprintf(stderr, "Number of alleles: %zu\n", allele_names.size());
    if (!is_diploid) {
        fprintf(stderr, "Haploid profiles used\n");
    }
    else if (is_diploid) {
        fprintf(stderr, "Diploid profiles used\n");
    }
    fprintf(stderr, "Input coverage: %f X, sequencing error: %f %%\n", coverage, sequencing_error);
    fprintf(stderr, "Input lambda error: %f\n", lambda_error);
    fprintf(stderr, "Lower k value: %zu, upper k value: %zu, k increment: %zu\n", lower_k, upper_k, increment_k);
    if (manual_lambda != -1) {
        fprintf(stderr, "Manual lambda value selected: %f\n", manual_lambda);
    }
    else {
        fprintf(stderr, "Selected lambda method: %s\n", lambda_method.c_str());
    }


    //
    // Determine all possible genotypes if needed
    //

    set<string> genotype_names;
    vector<pair<string, string>> genotypes;

    if (is_diploid) {
        genotypes = pairwise_comparisons(allele_names, true);
        for (auto iter = genotypes.begin(); iter != genotypes.end(); ++iter) {
            string name = iter->first + "/" + iter->second;
            genotype_names.insert(name);
        }
    }


    //
    // Iterate over all values of k
    //

    // Get all values of k
    vector<size_t> k_values;

    size_t counter = lower_k;
    while (counter < upper_k) {
        k_values.push_back(counter);
        counter += increment_k;
    }

    k_values.push_back(upper_k);


    FILE * output;
    output = fopen(output_name.c_str(), "w");

    // Iterate
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (size_t k = 0; k < k_values.size(); ++k) {


        //
        // Count k-mers in alleles/genotypes
        //

        // All k-mers in each allele/genotype
        map<string, kmer_count_map> allele_kmer_counts;
        map<string, kmer_count_map> genotype_kmer_counts;


        // Union of k-mers from all alleles
        unordered_set<string> allele_kmers;

        // Iterate over each allele
        for (size_t a = 0; a < alleles.size(); ++a) {
            map<string, size_t> single_allele_kmer_counts = count_kmers(alleles[a].sequence, k_values[k]);
            allele_kmer_counts[alleles[a].name.c_str()] = single_allele_kmer_counts;
            for (auto iter = single_allele_kmer_counts.begin(); iter != single_allele_kmer_counts.end(); ++iter) {
                allele_kmers.insert(iter->first);
            }
        }

        // Combine profiles for diploid genotypes if applicable
        if (is_diploid) {
            for (auto iter1 = genotypes.begin(); iter1 != genotypes.end(); ++iter1) {
                kmer_count_map combined_alelle_map = allele_kmer_counts[iter1->first];
                for (auto iter2 = allele_kmer_counts[iter1->second].begin(); iter2 != allele_kmer_counts[iter1->second].end(); ++iter2) {
                    combined_alelle_map[iter2->first] += iter2->second;
                }
                string genotype_name = iter1->first + "/" + iter1->second;
                genotype_kmer_counts[genotype_name] = combined_alelle_map;
            }
        }


        //
        // Count k-mers in reads
        //

        // k-mers and counts for each read
        map<string, kmer_count_map> each_read_kmer_counts;

        // k-mers and counts across all reads
        kmer_count_map all_reads_kmer_counts;

        // Iterate over each read
        for (size_t r = 0; r < reads1.size(); ++r) {
            if (reads1[r].sequence.length() > k_values[k]) {
                map<string, size_t> single_read_kmer_counts = count_kmers(reads1[r].sequence, k_values[k]);
                each_read_kmer_counts[reads1[r].name.c_str()] = single_read_kmer_counts;
                for (auto iter = single_read_kmer_counts.begin(); iter != single_read_kmer_counts.end(); ++iter) {
                    all_reads_kmer_counts[iter->first] += iter->second;
                }
            }
        }
        if (!input_reads_file2.empty()) {
            for (size_t r = 0; r < reads2.size(); ++r) {
                if (reads2[r].sequence.length() > k_values[k]) {
                    map<string, size_t> single_read_kmer_counts = count_kmers(reads2[r].sequence, k_values[k]);
                    each_read_kmer_counts[reads2[r].name.c_str()] = single_read_kmer_counts;
                    for (auto iter = single_read_kmer_counts.begin(); iter != single_read_kmer_counts.end(); ++iter) {
                        all_reads_kmer_counts[iter->first] += iter->second;
                    }
                }
            }
        }


        //
        // Calculate lambda from coverage, error rate, read length, and k
        //

        double calculated_lambda = calculate_lambda_cortex(read_length, k_values[k], coverage, sequencing_error);
        // Adjust lambda for diploid calling
        if (is_diploid) {
            // each allele contributes to half of the coverage
            calculated_lambda = calculated_lambda / 2;
            manual_lambda = manual_lambda / 2;
        }


        //
        // Calculate lambda from flanking kmer counts if necessary
        //

        double estimated_lambda_mean = 0;
        double estimated_lambda_median = 0;
        if (lambda_method != "coverage") {
            // Iterate over each flank sequence and combine
            kmer_count_map combined_flanks_counts;
            for (size_t f = 0; f < flanks.size(); ++f) {
                kmer_count_map single_flank_kmer_counts = count_kmers(flanks[f].sequence, k_values[k]);
                if (f == 0) {
                    combined_flanks_counts = single_flank_kmer_counts;
                }
                else {
                    for (auto iter = single_flank_kmer_counts.begin(); iter != single_flank_kmer_counts.end(); ++iter) {
                        combined_flanks_counts[iter->first] += iter->second;
                    }
                }
            }

            // Get vector of k-mers in flank sequences
            vector<string> flank_kmers;
            for (auto iter = combined_flanks_counts.begin(); iter != combined_flanks_counts.end(); ++iter) {
                flank_kmers.push_back(iter->first);
            }

            // Get counts from k-mers unique to flank sequences
            kmer_count_map flank_unqiue_kmer_counts;
            for (auto iter = flank_kmers.begin(); iter != flank_kmers.end(); ++iter) {
                string flank_kmer = *iter;
                auto flank_iter = allele_kmers.find(flank_kmer);
                if (flank_iter == allele_kmers.end()) {
                    flank_unqiue_kmer_counts[flank_kmer] = combined_flanks_counts[flank_kmer];
                }
            }
            if (flank_unqiue_kmer_counts.size() == 0) {
                fprintf(stderr, "No k-mers are unique to the flank sequences. Skipping this value of k. Consider selecting 'coverage' for lambda calculation method instead.\n");
                continue;
            }

            // Get counts for flank-unique k-mers in reads
            kmer_count_map all_reads_flank_kmer_counts;
            for (auto iter = flank_unqiue_kmer_counts.begin(); iter != flank_unqiue_kmer_counts.end(); ++iter) {
                string flank_kmer = iter->first;
                auto read_iter = all_reads_kmer_counts.find(flank_kmer);
                size_t flank_kmer_count_in_read = read_iter != all_reads_kmer_counts.end() ? read_iter->second : 0;
                all_reads_flank_kmer_counts[flank_kmer] = flank_kmer_count_in_read;
            }

            // Calculate estimated lambda
            double estimated_lambda_sum = 0;
            vector<double> read_flank_kmers_lambda;
            for (auto iter = flank_unqiue_kmer_counts.begin(); iter != flank_unqiue_kmer_counts.end(); ++iter) {
                // k-mer count in reads / k-mer count in flank sequences
                // TODO: add dummy count? adjust after or keep?
                double est_lambda = all_reads_flank_kmer_counts[iter->first] / iter->second;
                estimated_lambda_sum += est_lambda;
                read_flank_kmers_lambda.push_back(est_lambda);
            }
            // Lambda mean
            estimated_lambda_mean = estimated_lambda_sum / read_flank_kmers_lambda.size();
            // Lambda median
            sort(read_flank_kmers_lambda.begin(), read_flank_kmers_lambda.end());
            size_t median_position = read_flank_kmers_lambda.size() / 2;
                // value is rounded down if vector size is odd -- correct position for 0-based indexing
            if ((read_flank_kmers_lambda.size() % 2) == 0) {
                estimated_lambda_median = (read_flank_kmers_lambda[median_position] + read_flank_kmers_lambda[median_position - 1]) / 2;
            }
            else {
                estimated_lambda_median = read_flank_kmers_lambda[median_position];
            }

            if ((estimated_lambda_median == 0) & (lambda_method == "median")) {
                fprintf(stderr, "The median count for flank k-mers is 0. Skipping this value of k. Consider selecting 'mean' or 'coverage' for lambda calculation method instead.\n");
                continue;
            }
            else if ((estimated_lambda_mean == 0) & (lambda_method == "mean")) {
                fprintf(stderr, "The mean count for flank k-mers is 0. Skipping this value of k. Consider selecting 'coverage' for lambda calculation method instead.\n");
                continue;
            }
            

            // Adjust lambda for diploid calling
            if (is_diploid) {
                // each allele contributes to half of the coverage
                estimated_lambda_mean = estimated_lambda_mean / 2;
                estimated_lambda_median = estimated_lambda_median / 2;
            }
        }


        //
        // Select lambda for scoring
        //

        double lambda;
        if (manual_lambda > 0) {
            lambda = manual_lambda;
        }
        else if (lambda_method == "coverage") {
            lambda = calculated_lambda;

        }
        else if (lambda_method == "mean") {
            lambda = estimated_lambda_mean;

        }
        else {      // median
            lambda = estimated_lambda_median;

        }


        //
        // Print more handy information
        //

        fprintf(stderr, "---\nInformation for k = %zu\n", k_values[k]);
        fprintf(stderr, "Lambda calculated from sequence coverage, error, read length: %f\n", calculated_lambda);
        if (lambda_method == "mean" || lambda_method == "median") {
            fprintf(stderr, "Lambda calculated from flank kmer counts: mean: %f, median: %f\n", estimated_lambda_mean, estimated_lambda_median);
        }
        if (manual_lambda != -1 || manual_lambda != -0.5) {
            fprintf(stderr, "Manual lambda provided: %f\n", manual_lambda);
        }


        //
        // Score k-mer count profiles for each allele/genotype
        //

        map<string, double> all_scores;
        if (!is_diploid) {
            for (auto iter = allele_names.begin(); iter != allele_names.end(); ++iter) {
                string a = *iter;
                all_scores[a] = score_profile(all_reads_kmer_counts, allele_kmer_counts[a], allele_kmers, lambda, lambda_error);
            }
        }
        else if (is_diploid) {
            for (auto iter = genotype_names.begin(); iter != genotype_names.end(); ++iter) {
                string g = *iter;
                all_scores[g] = score_profile(all_reads_kmer_counts, genotype_kmer_counts[g], allele_kmers, lambda, lambda_error);
            }
        }


        // Pearson scoring
        // for (auto iter = allele_names.begin(); iter != allele_names.end(); ++iter) {
        //          string a = *iter;
        //          all_scores[a] = score_pearson(all_reads_kmer_counts, allele_kmer_counts[a], allele_kmers);
        //     }


        //
        // Save output
        //

        #pragma omp critical
        {
            // Create vector to store score:genotype values, sort and print top N
            vector<pair<double, string>> all_scores_vector;
            for (auto iter = all_scores.begin(); iter != all_scores.end(); ++iter) {
                all_scores_vector.push_back(make_pair(iter->second, iter->first));
            }

            sort(all_scores_vector.begin(), all_scores_vector.end(), greater<pair<double, string>>());
            int N_printed = 0;
            for (size_t i = 0; i < all_scores_vector.size(); ++i) {
                fprintf(output, "%zu,%s,%f\n", k_values[k], all_scores_vector[i].second.c_str(), all_scores_vector[i].first);
                N_printed += 1;
                if (top_N > 0 && N_printed == top_N) {
                    break;
                }
            }

            // Re-adjust manual lambda if diploid so it can be recalculated correctly for next k
            if (is_diploid) {
                manual_lambda = manual_lambda * 2;
            }
        }
    }

    fclose(output);


    //
    // Finished
    //

    return 0;
}