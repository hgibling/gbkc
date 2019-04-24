// count - Score reads based on comparison of k-mer count profiles to those for known alleles

#include <cmath>
#include <getopt.h>
#include <iostream>
#include <map>
#include <math.h>
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


//
// Define functions
//

// Calculate lambda
double calculate_lambda(double read_length, size_t k, double coverage, double sequencing_error)
{
    double lambda = (read_length - k + 1) * (coverage / read_length) * pow((1 - sequencing_error), k);
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
double log_poisson_pmf(double c, double lambda)
{
    double f_c = log_factorial(c);
    double p = (double)c * log(lambda) - lambda - f_c;
    return p;
}

// Score reads k-mer count against allele k-mer count (individual k-mers)
double score_kmer(size_t read_count, size_t allele_count, double lambda, double lambda_error)
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
double score_profile(const kmer_count_map& read_map, const kmer_count_map& allele_map, const std::unordered_set<std::string>& allele_union_kmers, double lambda, double lambda_error)
{
    double score = 0;
    for (auto iter = allele_union_kmers.begin(); iter != allele_union_kmers.end(); ++iter) {
        std::string kmer = *iter;
        auto rc_iter = read_map.find(kmer);
        size_t kmer_count_in_read = rc_iter != read_map.end() ? rc_iter->second : 0;
        auto ac_iter = allele_map.find(kmer);
        size_t kmer_count_in_allele = ac_iter != allele_map.end() ? ac_iter->second : 0;
        score += score_kmer(kmer_count_in_read, kmer_count_in_allele, lambda, lambda_error); 
    }
    return score;
}


//
// Help message
//

static const char *COUNT_USAGE_MESSAGE =
"Score reads based on comparison of k-mer count profiles to those for known alleles.\n\n"
"Usage: gbkc count [options]\n\n"
"Commands:\n"
"		-a       multi-fasta file of alleles/haplotypes of interest\n"
"		-1       multi-fasta/q file of sequencing reads to score (first in pair)\n"
"		-2       [optional] multi-fasta/q file of sequencing reads to score (second in pair)\n"
"		-k       size of k-mers to use\n"
"		-l       read length\n"
"		-e       sequencing error rate\n"
"		-c       sequencing coverage\n"
"		-m       error rate for lambda (default: 1)\n"
"		-o       output file name (default: results.csv)\n";


//
// Program commands
//

int countMain(int argc, char** argv) {
	
	if(argc <= 1) {
        std::cout << COUNT_USAGE_MESSAGE;
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
    double lambda_error = 1;
    std::string output_name = "results.csv";

    for (char c; (c = getopt_long(argc, argv, "a:1:2:k:l:e:c:m:o:", NULL, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'a': arg >> input_alleles_file; break;
            case '1': arg >> input_reads_file1; break;
            case '2': arg >> input_reads_file2; break;
            case 'k': arg >> input_k; break;
            case 'l': arg >> read_length; break;
            case 'e': arg >> sequencing_error; break;
            case 'c': arg >> coverage; break;
            case 'm': arg >> lambda_error; break;
            case 'o': arg >> output_name; break;
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
    std::vector<sequence_record> reads2;
	if (!input_reads_file2.empty()) {
        reads2 = read_sequences_from_file(input_reads_file2);
    }


    //
    // Calculate lambda
    //

    double lambda = calculate_lambda(read_length, input_k, coverage, sequencing_error);


    //
    // Print handy information
    //

    fprintf(stderr, "input reads: %s", input_reads_file1.c_str());
    fprintf(stderr, " %s", input_reads_file2.c_str());
    fprintf(stderr, "\ninput alleles: %s\n", input_alleles_file.c_str());
    fprintf(stderr, "number of alleles: %zu\n", allele_names.size());
    fprintf(stderr, "input k value: %zu\n", input_k);
    fprintf(stderr, "input coverage: %f X, sequencing error: %f %%\n", coverage, sequencing_error);
    fprintf(stderr, "input lambda error: %f\n", lambda_error);
    fprintf(stderr, "lambda calculated as: %f\n", lambda);


    //
    // Count k-mers in alleles
    //

    // All k-mers in each allele
    std::map<std::string, kmer_count_map> allele_kmer_counts; 

    // Union of k-mers from all alleles
    std::unordered_set<std::string> allele_kmers;

    // Iterate over each allele
    for (size_t a = 0; a < alleles.size(); ++a) {
        std::map<std::string, size_t> single_allele_kmer_counts = count_kmers(alleles[a].sequence, input_k);
        allele_kmer_counts[alleles[a].name.c_str()] = single_allele_kmer_counts;
        for (auto iter = single_allele_kmer_counts.begin(); iter != single_allele_kmer_counts.end(); ++iter) {
            allele_kmers.insert(iter->first);
        }
    }


    //
    // Count k-mers in reads
    //

    // k-mers and counts for each read
    std::map<std::string, kmer_count_map> each_read_kmer_counts;

    // k-mers and counts across all reads
    kmer_count_map all_reads_kmer_counts;

    // Union of k-mers across all reads (do I need this?)
    std::unordered_set<std::string> reads_kmers;

    // Iterate over each read
    for (size_t r = 0; r < reads1.size(); ++r) {
        std::map<std::string, size_t> single_read_kmer_counts = count_kmers(reads1[r].sequence, input_k);
        each_read_kmer_counts[reads1[r].name.c_str()] = single_read_kmer_counts;
        for (auto iter = single_read_kmer_counts.begin(); iter != single_read_kmer_counts.end(); ++iter) {
            reads_kmers.insert(iter->first);
            all_reads_kmer_counts[iter->first] += iter->second;
        }
    }
    if (!input_reads_file2.empty()) {
        for (size_t r = 0; r < reads2.size(); ++r) {
            std::map<std::string, size_t> single_read_kmer_counts = count_kmers(reads2[r].sequence, input_k);
            each_read_kmer_counts[reads2[r].name.c_str()] = single_read_kmer_counts;
            for (auto iter = single_read_kmer_counts.begin(); iter != single_read_kmer_counts.end(); ++iter) {
                reads_kmers.insert(iter->first);
                all_reads_kmer_counts[iter->first] += iter->second;
            }
        }  
    }   


    //
    // Score k-mer count profiles for each allele
    //

    std::map<std::string, double> all_scores;
    for (auto iter = allele_names.begin(); iter != allele_names.end(); ++iter) {
        std::string a = *iter;
        all_scores[a] = score_profile(all_reads_kmer_counts, allele_kmer_counts[a], allele_kmers, lambda, lambda_error);
    }


    //
    // Save output
    //

    FILE * output;
    output = fopen(output_name.c_str(), "w");

    for (auto iter = all_scores.begin(); iter != all_scores.end(); ++iter) {
        fprintf(output, "%s,%f\n", iter->first.c_str(), iter->second);
    }

    fclose(output);


    //
    // Finished
    //

	return 0;
}