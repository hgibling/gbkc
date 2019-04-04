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
#include <vector>
#include <zlib.h>
#include "kseq.h"



//
// Define functions
//

// Tell KSEQ (from Heng Li) what functions to use to open/read files
KSEQ_INIT(gzFile, gzread)

// Define sequence record structre
struct sequence_record
{
    std::string name;
    std::string sequence;
};

// Define k-mer count profile structure
typedef std::map<std::string, size_t> kmer_count_map;

// Get complement of a sequence
char complement(char c)
{
    switch(c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: 
            fprintf(stderr, "Error: unrecognized nucleotide %c\n", c);
            exit(EXIT_FAILURE);
    }
}

// Get reverse complement of a sequence
std::string reverse_complement(const std::string& sequence)
{
    size_t l = sequence.size();
    std::string rc(l, 'N');
    for (size_t i = 0; i < l; ++i) {
        rc[i] = complement(sequence[l - i - 1]);
    }
    return rc;
}

// Get the lexicographically lowest k-mer between it and its reverse complement
std::string canonical_kmer(const std::string& kmer)
{
    std::string rc_kmer = reverse_complement(kmer);
    return kmer < rc_kmer ? kmer : rc_kmer;
}

// Obtain all (canonical) k-mers and their counts from a sequence
kmer_count_map count_kmers(const std::string& sequence, size_t k)
{
    kmer_count_map out_map;
    for (size_t i = 0; i < sequence.size() - k + 1; ++i) {
        std::string canon_kmer = canonical_kmer(sequence.substr(i, k));
        out_map[canon_kmer] += 1;
    }
    return out_map;
}   

// Read input files
std::vector<sequence_record> read_sequences_from_file(const std::string& input_filename)
{
    // Open readers
    FILE* read_fp = fopen(input_filename.c_str(), "r");
    if(read_fp == NULL) {
        fprintf(stderr, "error: could not open %s for read\n", input_filename.c_str());
        exit(EXIT_FAILURE);
    }

    gzFile gz_read_fp = gzdopen(fileno(read_fp), "r");
    if(gz_read_fp == NULL) {
        fprintf(stderr, "error: could not open %s using gzdopen\n", input_filename.c_str());
        exit(EXIT_FAILURE);
    }
    
    std::vector<sequence_record> out_sequences;
    int ret = 0;
    kseq_t* seq = kseq_init(gz_read_fp);
    while((ret = kseq_read(seq)) >= 0) {
        sequence_record sr = { seq->name.s, seq->seq.s };
        out_sequences.push_back(sr);
    }    
    kseq_destroy(seq);   

    return out_sequences;
}

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
// Program commands
//

int main(int argc, char** argv) {

    //
    // Read command line arguments
    //
    std::string input_alleles_file;
    std::string input_reads_file;
    size_t input_k = 21;
    double read_length = 0;
    double sequencing_error = 0;
    double coverage = 0;
    double lambda_error = 1;
    std::string output_name = "results.csv";

    for (char c; (c = getopt_long(argc, argv, "a:r:k:l:e:c:m:o:", NULL, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'a':
                arg >> input_alleles_file;
                break;
            case 'r':
                arg >> input_reads_file;
                break;
            case 'k':
                arg >> input_k;
                break;
            case 'l':
                arg >> read_length;
                break;
            case 'e':
                arg >> sequencing_error;
                break;
            case 'c':
                arg >> coverage;
                break;
            case 'm':
                arg >> lambda_error;
                break;
            case 'o':
                arg >> output_name;
                break;
            default:
                exit(EXIT_FAILURE);
        }
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
    fprintf(stderr, "number of alleles: %zu\n", allele_names.size());

    // Get read sequences
    std::vector<sequence_record> reads = read_sequences_from_file(input_reads_file);


    //
    // Calculate lambda
    //
    double lambda = calculate_lambda(read_length, input_k, coverage, sequencing_error);


    //
    // Print handy information
    //
    fprintf(stderr, "input reads: %s\n", input_reads_file.c_str());
    fprintf(stderr, "input alleles: %s\n", input_alleles_file.c_str());
    fprintf(stderr, "input k value: %zu\n", input_k);
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
    for (size_t r = 0; r < reads.size(); ++r) {
        std::map<std::string, size_t> single_read_kmer_counts = count_kmers(reads[r].sequence, input_k);
        each_read_kmer_counts[reads[r].name.c_str()] = single_read_kmer_counts;
        for (auto iter = single_read_kmer_counts.begin(); iter != single_read_kmer_counts.end(); ++iter) {
            reads_kmers.insert(iter->first);
            all_reads_kmer_counts[iter->first] += iter->second;
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

    return 0;
}
