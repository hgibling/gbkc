// check-profiles - Compare k-mer count profiles for all known alleles to check if they are unique

#include <cmath>
#include <getopt.h>
#include <iostream>
#include <math.h>
#include <set>
#include <sstream>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unordered_set>
#include <utility>
#include <zlib.h>

#include "check-profiles.h"
#include "kseq.h"


//
// Define functions
//

// Tell KSEQ (from Heng Li) what functions to use to open/read files
KSEQ_INIT(gzFile, gzread)


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
        fprintf(stderr, "error: could not open file %s. Check file name.\n", input_filename.c_str());
        exit(EXIT_FAILURE);
    }

    gzFile gz_read_fp = gzdopen(fileno(read_fp), "r");
    if(gz_read_fp == NULL) {
        fprintf(stderr, "error: could not open file %s using gzdopen. Check file name.\n", input_filename.c_str());
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

// Compare k-mer count profiles
// from https://stackoverflow.com/questions/8473009/how-to-efficiently-compare-two-maps-of-strings-in-c (user sebastian-mach)
template <typename Map>
bool compare_profiles (Map const &lhs, Map const &rhs)
{
    return lhs.size() == rhs.size()
        && std::equal(lhs.begin(), lhs.end(),
                      rhs.begin());
}

// Generate all possible comparisons of two alleles
std::vector<std::pair<std::string, std::string>> pairwise_comparisons(const std::set<std::string>& alleles) 
{
    std::vector<std::pair<std::string, std::string>> pairs;
    std::vector<std::string> alleles_vector(alleles.begin(), alleles.end());
    for (size_t i = 0; i < alleles_vector.size() - 1; ++i) {
        for (size_t j = 1; j < alleles_vector.size(); ++j) {
            if (i != j){
                pairs.push_back(std::make_pair(alleles_vector[i], alleles_vector[j]));
            }
        }
    }
    return pairs;
}


//
// Help message
//

static const char *CHECK_PROFILES_USAGE_MESSAGE =
"Compare k-mer count profiles for all known alleles to check if they are unique.\n\n"
"Usage: gbkc check-profiles [options]\n\n"
"Commands:\n"
"       -a       multi-fasta file of alleles/haplotypes of interest\n"
"       -l       lower range of k values to be checked (default: 3)\n"
"       -u       upper range of k values to be checked (should not be greater than the read length)\n";


//
// Program commands
//

int checkprofilesMain(int argc, char** argv) {

    if(argc <= 1) {
        std::cout << CHECK_PROFILES_USAGE_MESSAGE;
        return 0;
    };


    //
    // Read command line arguments
    //

    std::string input_alleles_file;
    int lower_range = 3;
    int upper_range = -1;

    for (char c; (c = getopt_long(argc, argv, "a:l:u:", NULL, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'a': arg >> input_alleles_file; break;
            case 'l': arg >> lower_range; break;
            case 'u': arg >> upper_range; break;
            default: exit(EXIT_FAILURE);
        }
    }

    // Check for parameter arguments
    if (input_alleles_file.empty()) {
        fprintf(stderr, "No file for allele sequences. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (upper_range <= 0 || upper_range < lower_range) {
        fprintf(stderr, "Upper range for testing k-mer count profiles must be greater than 0 and greater than the lower range. Check parameters.\n");
        exit(EXIT_FAILURE);
    }


    //
    // Read file
    //

    // TODO: Handle case where file names are null
    
    // Get allele sequences
    std::vector<sequence_record> alleles = read_sequences_from_file(input_alleles_file);

    // Get list of allele names
    std::set<std::string> allele_names;
    for (size_t i = 0; i < alleles.size(); ++i) {
        allele_names.insert(alleles[i].name.c_str());
    }


    //
    // Compare k-mer count profiles
    //

    // Get list of allele comparisons to make
    std::vector<std::pair<std::string, std::string>> allele_pairs = pairwise_comparisons(allele_names);

    // Iterate over all possible values of k given read length l
    for (int k = lower_range; k < upper_range + 1; ++k) {
        printf("Testing k = %d... ", k);

        // All k-mers in each allele
        std::map<std::string, kmer_count_map> allele_kmer_counts; 

        // Get count profiles for each allele
        for (size_t a = 0; a < alleles.size(); ++a) {
            std::map<std::string, size_t> single_allele_kmer_counts = count_kmers(alleles[a].sequence, k);
            allele_kmer_counts[alleles[a].name.c_str()] = single_allele_kmer_counts;
        }

        // Check if k-mer count profiles are unique amongst all alleles
        std::vector<std::pair<std::string, std::string>> identical_profiles;

        for (auto iter = allele_pairs.begin(); iter != allele_pairs.end(); ++iter) {
            bool allele_vs_allele = compare_profiles(allele_kmer_counts[iter->first], allele_kmer_counts[iter->second]);
            if (allele_vs_allele == 1) {
                identical_profiles.push_back(std::make_pair(iter->first, iter->second));
            }
        }

        if (identical_profiles.size() > 0) {
            printf("Some alleles had identical k-mer count profiles:\n");
            for (auto iter = identical_profiles.begin(); iter != identical_profiles.end(); ++iter) {
                printf("%s and %s\n", iter->first.c_str(), iter->second.c_str());
            }
        } 
        else {
            printf("All k-mer count profiles are unique\n");
        }

    }


    //
    // Finished
    //

    return 0;
}