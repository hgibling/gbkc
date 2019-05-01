// distance - Score reads based on comparison of k-mer pair distances to those for known alleles

#include <cmath>
#include <getopt.h>
#include <iostream>
#include <map>
#include <math.h>
#include <numeric>
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
std::pair<std::string, std::string> outer_kmers(const std::string& first_kmer, const std::string& second_kmer, const size_t k, const size_t p)
{
    std::pair<std::string, std::string> out_pair;
    if (p == 0) {
        out_pair = make_pair(reverse_complement(first_kmer).substr(0, k), first_kmer.substr(0, k));
    }
    else if (p == 1) {
        out_pair = make_pair(first_kmer.substr(0, k), reverse_complement(first_kmer).substr(0, k));
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

// Score reads outer k-mers with allele k-mer outer distance
double score_kmer_distances(const std::pair<std::string, std::string> read_kmer_pair, const kmer_position_map allele_distances, const double fragment_length, const double fragment_stdev, const size_t k, const size_t penalty)
{
    double final_score = 0;
    std::vector<double> kmer_pair_scores;
    auto check1 = allele_distances.find(read_kmer_pair.first);
    auto check2 = allele_distances.find(read_kmer_pair.second);
    if ((check1 == allele_distances.end() || check2 == allele_distances.end())) {
        final_score = log_normal_pdf(penalty, fragment_length, fragment_stdev);
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
        // Sum each log probability
        final_score = std::accumulate(kmer_pair_scores.begin(), kmer_pair_scores.end(), 0);
    }
    return final_score;
}


// auto read_iter = read_map.find(kmer);
//         size_t kmer_count_in_read = read_iter != read_map.end() ? read_iter->second : 0;


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
"       -k       size of k-mers to use\n"
"       -l       read length\n"
"       -e       sequencing error rate (between 0 and 1)\n"
"       -c       sequencing coverage\n"
"       -f       mean fragment length\n"
"       -s       standard deviation of fragment length\n"
"       -p       penalty fragment length when k-mer pairs aren't observed in an allele (default: 10000000)\n"
"       -o       output file name (default: results.csv)\n";


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
    size_t input_k = 0;
    double read_length = -1;
    double sequencing_error = -1;
    double coverage = -1;
    double fragment_length = -1;
    double fragment_stdev = -1;
    size_t input_penalty = 10000000;
    std::string output_name = "distance-results.csv";

    for (char c; (c = getopt_long(argc, argv, "a:1:2:k:l:e:c:f:s:o:p:", NULL, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'a': arg >> input_alleles_file; break;
            case '1': arg >> input_reads_file1; break;
            case '2': arg >> input_reads_file2; break;
            case 'k': arg >> input_k; break;
            case 'l': arg >> read_length; break;
            case 'e': arg >> sequencing_error; break;
            case 'c': arg >> coverage; break;
            case 'f': arg >> fragment_length; break;
            case 's': arg >> fragment_stdev; break;
            case 'p': arg >> input_penalty; break;
            case 'o': arg >> output_name; break;
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

    if (input_k == 0 || input_k > read_length) {
        fprintf(stderr, "Value of k must be non-zero and less than the read length. Check parameters.\n");
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
        fprintf(stderr, "Penalty for when k-mer pairs don't exist in an allele must be greater than 0 (and should be a very large number, e.g. 1000000. Check parameters.\n");
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

    fprintf(stderr, "input reads: %s and %s\n", input_reads_file1.c_str(), input_reads_file2.c_str());
    fprintf(stderr, "input alleles: %s\n", input_alleles_file.c_str());
    fprintf(stderr, "number of alleles: %zu\n", allele_names.size());
    fprintf(stderr, "input k value: %zu\n", input_k);
    fprintf(stderr, "input coverage: %f X, sequencing error: %f %%\n", coverage, sequencing_error);
    fprintf(stderr, "input mean fragment length: %f, standard deviation: %f\n", fragment_length, fragment_stdev);
    fprintf(stderr, "input penalty: %zu\n", input_penalty);


    //
    // Get k-mer positions for alleles
    //

    std::map<std::string, kmer_position_map> allele_positions;
    for (size_t a = 0; a < alleles.size(); ++a) {
        kmer_position_map single_allele_positions = kmer_positions(alleles[a].sequence, input_k);
        allele_positions[alleles[a].name] = single_allele_positions;
    }


    //
    // Get outer k-mer pairs for each read pair
    //

    // Assuming reads are in paired order in both lists
    std::multimap<std::string, std::pair<std::string, std::string>> read_pairs;
    for (size_t r = 0; r < reads1.size(); ++r) {

        // Get both first read kmer and reverse complement second read kmer, and reverse complement first read kmer
        // and second read kmer, since we don't know orientation of the reads
        std::pair<std::string, std::string> first_rc = outer_kmers(reads1[r].sequence, reads2[r].sequence, input_k, 0);
        std::pair<std::string, std::string> second_rc = outer_kmers(reads1[r].sequence, reads2[r].sequence, input_k, 1);

        read_pairs.insert({reads1[r].name, first_rc});
        read_pairs.insert({reads1[r].name, second_rc});
    }


    //
    // Calculate distance scores for each allele
    //

    std::pair<std::string, std::string> kmerpair ("AAA", "TAT");

    std::multimap<std::string, size_t> dict {
        {"AAA", 0},
        {"AAA", 1},
        {"TAT", 2},
        {"AGT", 3},
        {"TAT", 4},
        {"AAA", 5}
    };

    

    double out = score_kmer_distances(kmerpair, dict, fragment_length, fragment_stdev, input_k, input_penalty);

    
    fprintf(stderr, "score: %f\n", out);

    std::vector<size_t> tryme = {1,4,6,9,22};
    size_t trymesum = std::accumulate(tryme.begin(), tryme.end(), 0);
    fprintf(stderr, "TEST: %zu\n", trymesum);
    


    //
    // Finished
    //

    return 0;
}