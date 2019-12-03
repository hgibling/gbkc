// distance - Score reads based on comparison of k-mer pair distances to those for known alleles

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
#include "distance.h"
#include "kseq.h"


//
// Define functions
//

// Get positions for every k-mer in a sequence
kmer_position_map kmer_positions(const std::string& sequence, const size_t k)
{
    kmer_position_map out_map;
    for (size_t i = 0; i < sequence.size() - k + 1; ++i) {
        out_map.insert(std::pair<std::string, size_t>(sequence.substr(i, k), i));
    }
    return out_map;
}

// Get pairs of outermost k-mers from a read pair
std::pair<std::string, std::string> outer_kmers(const std::string& first_read, const std::string& second_read, const size_t k, const size_t p)
{
    std::pair<std::string, std::string> out_pair;
    if (p == 0) {
        out_pair = make_pair(first_read.substr(0, k), reverse_complement(second_read.substr(0, k)));
    }
    else if (p == 1) {
        out_pair = make_pair(second_read.substr(0, k), reverse_complement(first_read.substr(0, k)));
    }
    else {
        fprintf(stderr, "outer_kmer position must be 0 or 1.\n");
        exit(EXIT_FAILURE);
    }
    return out_pair;
}

// Get outer distance between two k-mers given their positions
size_t outer_distance(const size_t first_kmer_position, const size_t second_kmer_position, const size_t k, const size_t penalty)
{
    size_t distance;
    distance = (second_kmer_position >= first_kmer_position) ? ((second_kmer_position + k) - first_kmer_position) : penalty;
    return distance;
}

// Log normal distribution probability density function
double log_normal_pdf(const double distance, const double mean, const double standard_dev)
{
    double variance = pow(standard_dev, 2);
    double pi2 = M_PI * 2;
    double probability = (- log(sqrt(variance * pi2))) - (pow((distance - mean), 2)) / (2 * variance);
    return probability;
}

// Get the geometric mean of a vector of numbers
double geometric_mean(const std::vector<double>& numbers)
{
    double sum = 0;
    size_t size = numbers.size();
    for (auto iter = numbers.begin(); iter != numbers.end(); ++iter) {
        sum += log(*iter);
    }
    double geomean = exp(sum/size); 
    return geomean;
}

// Score reads outer k-mers with allele k-mer outer distance
double score_kmer_distances(const std::pair<std::string, std::string>& read_kmer_pair, const kmer_position_map& allele_distances, const double fragment_length, const double fragment_stdev, const size_t k, const size_t penalty, const std::string& method)
{
    double score = 0;
    std::vector<double> kmer_pair_scores;
    auto check1 = allele_distances.find(read_kmer_pair.first);
    auto check2 = allele_distances.find(read_kmer_pair.second);
    // If either kmer is not present in the list of allele kmer distances, use the penalty score instead
    if ((check1 == allele_distances.end() || check2 == allele_distances.end())) {
        score = log_normal_pdf(penalty, fragment_length, fragment_stdev);
    }
    else {
        auto range1 = allele_distances.equal_range(read_kmer_pair.first);
        auto range2 = allele_distances.equal_range(read_kmer_pair.second);
        for (auto iter1 = range1.first; iter1 != range1.second; ++iter1) {
            for (auto iter2 = range2.first; iter2 != range2.second; ++iter2) {
                size_t distance = outer_distance(iter1->second, iter2->second, k, penalty);
                kmer_pair_scores.push_back(log_normal_pdf(distance, fragment_length, fragment_stdev));
            }
        }
        if (method == "sum") {
            // Sum each log probability
            std::vector<double> antilog;
            for (auto iter = kmer_pair_scores.begin(); iter != kmer_pair_scores.end(); ++iter) {
                antilog.push_back(exp(*iter));
            }
            score = log(std::accumulate(antilog.begin(), antilog.end(), 0.0));
        }
        else if (method == "mean") {
            // Get the arithmetic mean of the log probabilities
            std::vector<double> antilog;
            for (auto iter = kmer_pair_scores.begin(); iter != kmer_pair_scores.end(); ++iter) {
                antilog.push_back(exp(*iter));
            }
            score = log(std::accumulate(antilog.begin(), antilog.end(), 0.0)/antilog.size());
        }
        else if (method == "geomean") {
            // Get the geometric mean of the log probabilities
            std::vector<double> antilog;
            for (auto iter = kmer_pair_scores.begin(); iter != kmer_pair_scores.end(); ++iter) {
                antilog.push_back(exp(*iter));
            }
            score = log(geometric_mean(antilog));
        }
        else if (method == "max") {
            // Take the maximum log probability
            score = *max_element(kmer_pair_scores.begin(), kmer_pair_scores.end());
        }
        else {
            exit(EXIT_FAILURE);
        }
    }
    return score;
}

// Score each read (for haploids)
double allele_score_read_kmer_pairs(const kmer_position_map& allele_positions, const kmer_pairs_map& read_pairs, const double fragment_length, const double fragment_stdev, const size_t k, const size_t penalty, const std::string& method)
{
    double score = 0;
    std::set<std::string> read_names;
    for (auto iter = read_pairs.begin(); iter != read_pairs.end(); ++iter) {
        read_names.insert(iter->first);
    }
    for (auto iter1 = read_names.begin(); iter1 != read_names.end(); ++iter1) {
        auto range = read_pairs.equal_range(*iter1);
        std::vector<double> compare_scores;
        for (auto iter2 = range.first; iter2 != range.second; ++iter2) {
            compare_scores.push_back(score_kmer_distances(iter2->second, allele_positions, fragment_length, fragment_stdev, k, penalty, method));
        }
        double max_score = *max_element(compare_scores.begin(), compare_scores.end());
        score += max_score;
    }
    return score;
}

// Score each read (for diploids)
double genotype_score_read_kmer_pairs(const kmer_position_map& allele_positions1, const kmer_position_map& allele_positions2, const kmer_pairs_map& read_pairs, const double fragment_length, const double fragment_stdev, const size_t k, const size_t penalty, const std::string& method)
{
    double score = 0;
    std::set<std::string> read_names;
    for (auto iter = read_pairs.begin(); iter != read_pairs.end(); ++iter) {
        read_names.insert(iter->first);
    }
    for (auto iter1 = read_names.begin(); iter1 != read_names.end(); ++iter1) {
        auto range = read_pairs.equal_range(*iter1);
        std::vector<double> compare_scores1;
        for (auto iter2 = range.first; iter2 != range.second; ++iter2) {
            compare_scores1.push_back(score_kmer_distances(iter2->second, allele_positions1, fragment_length, fragment_stdev, k, penalty, method));
        }
        double max_score1 = *max_element(compare_scores1.begin(), compare_scores1.end());

        std::vector<double> compare_scores2;
        for (auto iter2 = range.first; iter2 != range.second; ++iter2) {
            compare_scores2.push_back(score_kmer_distances(iter2->second, allele_positions2, fragment_length, fragment_stdev, k, penalty, method));
        }
        double max_score2 = *max_element(compare_scores2.begin(), compare_scores2.end());
        score += ((max_score1/2) + (max_score2/2));
    }
    return score;
}


//
// Help message
//

static const char *DISTANCE_USAGE_MESSAGE =
"Score reads based on comparison of k-mer pair distances to those for known alleles.\n\n"
"Usage: gbkc distance [options]\n\n"
"Commands:\n"
"       -a       multi-fasta file of alleles/haplotypes of interest\n"
"       -1       multi-fasta/q file of sequencing reads to score (first in pair)\n"
"       -2       multi-fasta/q file of sequencing reads to score (second in pair)\n"
"       -d       score diploid genotypes (if flag is not used, haploid scoring is performed)\n"
"       -k       lower value of k-mer to use (default: 11)\n"
"       -K       upper value of k-mer to use\n"
"       -r       increment size between k-mer values (default: 4)\n"
"       -l       read length\n"
"       -e       sequencing error rate (between 0 and 1)\n"
"       -c       sequencing coverage\n"
"       -f       mean fragment length\n"
"       -s       standard deviation of fragment length\n"
"       -p       penalty fragment length when k-mer pairs aren't observed in an allele (default: 10)\n"
"       -m       method for summarizing scores when kmer pairs occur more than once in an allele\n"
"       -o       output file name (default: results.csv)\n"
"       -t       number of threads (default: 1)\n";


//
// Program commands
//

int distanceMain(int argc, char** argv) {

    if(argc <= 1) {
        std::cout << DISTANCE_USAGE_MESSAGE;
        return 0;
    };


    //
    // Read command line arguments
    //

    std::string input_alleles_file;
    std::string input_reads_file1;
    std::string input_reads_file2;
    size_t lower_k = 11;
    size_t upper_k = 0;
    size_t increment_k = 4;
    double read_length = -1;
    double sequencing_error = -1;
    double coverage = -1;
    double fragment_length = -1;
    double fragment_stdev = -1;
    size_t input_penalty = 10;
    std::string method;
    std::string output_name = "distance-results.csv";
    bool is_diploid = false;
    size_t num_threads = 1;

    for (char c; (c = getopt_long(argc, argv, "a:1:2:k:K:r:l:e:c:f:s:p:m:o:dt:", NULL, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'a': arg >> input_alleles_file; break;
            case '1': arg >> input_reads_file1; break;
            case '2': arg >> input_reads_file2; break;
            case 'k': arg >> lower_k; break;
            case 'K': arg >> upper_k; break;
            case 'r': arg >> increment_k; break;
            case 'l': arg >> read_length; break;
            case 'e': arg >> sequencing_error; break;
            case 'c': arg >> coverage; break;
            case 'f': arg >> fragment_length; break;
            case 's': arg >> fragment_stdev; break;
            case 'p': arg >> input_penalty; break;
            case 'm': arg >> method; break;
            case 'o': arg >> output_name; break;
            case 'd': is_diploid = true; break;
            case 't': arg >> num_threads; break;
            default: exit(EXIT_FAILURE);
        }
    }


    // Check for correct parameter arguments
    if (input_alleles_file.empty()) {
        fprintf(stderr, "No file for allele sequences. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (input_reads_file1.empty() || input_reads_file2.empty()) {
        fprintf(stderr, "Missing file(s) for read sequences. Two files must be specified (one for -1, one for -2). Check parameters.\n");
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

    if (fragment_length <= 0 || fragment_stdev <=0) {
        fprintf(stderr, "Mean fragment length and standard deviation must be greater than 0. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (input_penalty == 0) {
        fprintf(stderr, "Penalty for when k-mer pairs don't exist in an allele must be greater than 0 (and should not be similar to the fragment size). Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (method != "sum" && method != "mean" && method != "geomean" && method != "max") {
        fprintf(stderr, "Method must be one of 'sum', 'mean', 'geomean', or 'max'. Check parameters.\n");
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
    std::vector<sequence_record> alleles = read_sequences_from_file(input_alleles_file);

    // Get list of allele names
    std::set<std::string> allele_names;
    for (size_t i = 0; i < alleles.size(); ++i) {
        allele_names.insert(alleles[i].name.c_str());
    }

    // Get read sequences
    std::vector<sequence_record> reads1 = read_sequences_from_file(input_reads_file1);
    std::vector<sequence_record> reads2 = read_sequences_from_file(input_reads_file2);
    if (reads1.size() != reads2.size()) {
        fprintf(stderr, "All reads must have a mate (need same number of reads in first and second list).\n");
        exit(EXIT_FAILURE);
    }


    //
    // Print handy information
    //

    fprintf(stderr, "Number of threads used: %zu\n", num_threads);
    fprintf(stderr, "Input reads: %s and %s\n", input_reads_file1.c_str(), input_reads_file2.c_str());
    fprintf(stderr, "Input alleles: %s\n", input_alleles_file.c_str());
    fprintf(stderr, "Number of alleles: %zu\n", allele_names.size());
    if (!is_diploid) {
        fprintf(stderr, "Haploid profiles used\n");
    }
    else if (is_diploid) {
        fprintf(stderr, "Diploid profiles used\n");
    }
    fprintf(stderr, "Lower k value: %zu, upper k value: %zu, k increment: %zu\n", lower_k, upper_k, increment_k);
    fprintf(stderr, "Input coverage: %f X, sequencing error: %f %%\n", coverage, sequencing_error);
    fprintf(stderr, "Input mean fragment length: %f, standard deviation: %f\n", fragment_length, fragment_stdev);
    fprintf(stderr, "Input penalty: %zu\n", input_penalty);
    fprintf(stderr, "Input scoring method: %s\n", method.c_str());


    //
    // Determine all possible genotypes if needed
    //

    std::vector<std::pair<std::string, std::string>> genotypes;

    if (is_diploid) {
        genotypes = pairwise_comparisons(allele_names, true);
    }


    //
    // Iterate over all values of k
    //

    // Get all values of k
    std::vector<size_t> k_values;

    size_t counter = lower_k;
    while (counter < upper_k) {
        k_values.push_back(counter);
        counter += increment_k;
    }
    
    k_values.push_back(upper_k);

    // Set up output
    FILE * output;
    output = fopen(output_name.c_str(), "w");


    // Iterate
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (size_t k = 0; k < k_values.size(); ++k) {


        //
        // Get k-mer positions for alleles
        //

        std::map<std::string, kmer_position_map> allele_positions;
        for (size_t a = 0; a < alleles.size(); ++a) {
            kmer_position_map single_allele_positions = kmer_positions(alleles[a].sequence, k_values[k]);
            allele_positions[alleles[a].name] = single_allele_positions;
        }


        //
        // Get outer k-mer pairs for each read pair
        //

        // Assuming reads are in paired order in both lists
        kmer_pairs_map read_pairs;
        for (size_t r = 0; r < reads1.size(); ++r) {

            // Get both first read kmer and reverse complement second read kmer, and reverse complement first read kmer
            // and second read kmer, since we don't know orientation of the reads
            std::pair<std::string, std::string> first_rc = outer_kmers(reads1[r].sequence, reads2[r].sequence, k_values[k], 0);
            std::pair<std::string, std::string> second_rc = outer_kmers(reads1[r].sequence, reads2[r].sequence, k_values[k], 1);

            read_pairs.insert({reads1[r].name, first_rc});
            read_pairs.insert({reads1[r].name, second_rc});
        }


        //
        // Calculate distance scores for each allele or genotype
        //

        std::map<std::string, double> all_scores;

        if (!is_diploid) {
            for (auto iter = allele_names.begin(); iter != allele_names.end(); ++iter) {
                std::string a = *iter;
                all_scores[a] = allele_score_read_kmer_pairs(allele_positions[a], read_pairs, fragment_length, fragment_stdev, k_values[k], input_penalty, method);
            }
        }
        else if (is_diploid) {
            for (auto iter = genotypes.begin(); iter != genotypes.end(); ++iter) {
                std::string genotype_name = iter->first + "/" + iter->second;
                all_scores[genotype_name] = genotype_score_read_kmer_pairs(allele_positions[iter->first], allele_positions[iter->second], read_pairs, fragment_length, fragment_stdev, k_values[k], input_penalty, method);
            }
        }


        //
        // Save output
        //

        #pragma omp critical
        {
            for (auto iter = all_scores.begin(); iter != all_scores.end(); ++iter) {
                fprintf(output, "%zu,%s,%f\n", k_values[k], iter->first.c_str(), iter->second);
            }
        }
    }

    fclose(output);


    //
    // Finished
    //

    return 0;
}