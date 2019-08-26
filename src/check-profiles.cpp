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
char complement(const char nucleotide)
{
    switch(nucleotide) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: 
            fprintf(stderr, "Error: unrecognized nucleotide %c\n", nucleotide);
            exit(EXIT_FAILURE);
    }
}

// Get reverse complement of a sequence
std::string reverse_complement(const std::string& sequence)
{
    size_t length = sequence.size();
    std::string out_reverse_complement(length, 'N');
    for (size_t i = 0; i < length; ++i) {
        out_reverse_complement[i] = complement(sequence[length - i - 1]);
    }
    return out_reverse_complement;
}

// Get the lexicographically lowest k-mer between it and its reverse complement
std::string canonical_kmer(const std::string& kmer)
{
    std::string reverse_complement_kmer = reverse_complement(kmer);
    return kmer < reverse_complement_kmer ? kmer : reverse_complement_kmer;
}

// Obtain all (canonical) k-mers and their counts from a sequence
kmer_count_map count_kmers(const std::string& sequence, const size_t k)
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
bool compare_profiles(const Map& first_map, const Map& second_map)
{
    return first_map.size() == second_map.size()
        && std::equal(first_map.begin(), first_map.end(),
                      second_map.begin());
}

// Generate count differences between two k-mer count profiles
kmer_comparison_map kmer_differences(const kmer_count_map& first_map, const kmer_count_map& second_map)
{
    kmer_comparison_map kmer_differences_out;
    std::set<std::string> all_kmers;
    for (auto iter1 = first_map.begin(); iter1 != first_map.end(); ++iter1) {
        all_kmers.insert(iter1->first);
    }
    for (auto iter2 = second_map.begin(); iter2 != second_map.end(); ++iter2) {
        all_kmers.insert(iter2->first);
    }
    for (auto iter = all_kmers.begin(); iter != all_kmers.end(); ++iter) {
        auto first_iter = first_map.find(iter->c_str());
        size_t first_count = first_iter != first_map.end() ? first_iter->second : 0;
        auto second_iter = second_map.find(iter->c_str());
        size_t second_count = second_iter != second_map.end() ? second_iter->second : 0;
        kmer_differences_out[iter->c_str()] = second_count - first_count;
    }

    return kmer_differences_out;
}

// Generate all possible comparisons of two alleles
std::vector<std::pair<std::string, std::string>> pairwise_comparisons(const std::set<std::string>& alleles, bool same) 
{
    std::vector<std::pair<std::string, std::string>> pairs;
    std::vector<std::string> alleles_vector(alleles.begin(), alleles.end());
    for (size_t i = 0; i < alleles_vector.size(); ++i) {
        for (size_t j = 0; j < alleles_vector.size(); ++j) {
            if ((same == true) and (i <= j)) {
                pairs.push_back(std::make_pair(alleles_vector[i], alleles_vector[j]));
            }
            else if (i < j) {
                pairs.push_back(std::make_pair(alleles_vector[i], alleles_vector[j]));
            }
        }
    }
    return pairs;
}

// Split a diploid genotype into a vector of its alleles
std::vector<std::string> genotype_split(const std::string genotype, const char delimiter)
{
    size_t position = genotype.find(delimiter);
    std::string allele1 = genotype.substr(0, position);
    std::string allele2 = genotype.substr(position + 1, genotype.size());
    std::vector<std::string> out = {allele1, allele2};
    return out;
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
"       -u       upper range of k values to be checked (should not be greater than the read length)\n"
"       -d       check diploid genotypes (if flag is not used, haploid check is performed)\n"
"       -v       verbose printout of comparisons when profiles are not identical (true or false; default: false)\n";


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
    std::string verbose = "false";
    bool is_diploid = false;
    bool is_verbose = false;


    for (char c; (c = getopt_long(argc, argv, "a:l:u:pv", NULL, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'a': arg >> input_alleles_file; break;
            case 'l': arg >> lower_range; break;
            case 'u': arg >> upper_range; break;
            case 'd': is_diploid = true; break;
            case 'v': is_verbose = true; break;
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

    // // Get list of allele comparisons to make
    std::vector<std::pair<std::string, std::string>> allele_pairs;
    std::vector<std::pair<std::string, std::string>> genotype_pairs;
    std::vector<std::pair<std::string, std::string>> genotypes;

    if (is_diploid == false) {
        allele_pairs = pairwise_comparisons(allele_names, false);
    }
    else if (is_diploid == true) {
        genotypes = pairwise_comparisons(allele_names, true);
        std::set<std::string> genotype_names;
        for (auto iter = genotypes.begin(); iter != genotypes.end(); ++iter) {
            std::string name = iter->first + "/" + iter->second;
            genotype_names.insert(name);
        }
        genotype_pairs = pairwise_comparisons(genotype_names, false);
    }


    // Iterate over all possible values of k given read length l
    for (int k = lower_range; k < upper_range + 1; ++k) {
        printf("Testing k = %d... \n", k);

        // All k-mers in each allele or genotype
        std::map<std::string, kmer_count_map> allele_kmer_counts; 
        std::map<std::string, kmer_count_map> genotype_kmer_counts;

        // Get count profiles for each allele
        for (size_t a = 0; a < alleles.size(); ++a) {
            kmer_count_map single_allele_kmer_counts = count_kmers(alleles[a].sequence, k);
            allele_kmer_counts[alleles[a].name] = single_allele_kmer_counts;
        }

        // Combine profiles for diploid genotypes if applicable
        if (is_diploid == true) {
            for (auto iter1 = genotypes.begin(); iter1 != genotypes.end(); ++iter1) {
                kmer_count_map combined_alelle_map = allele_kmer_counts[iter1->first];
                for (auto iter2 = allele_kmer_counts[iter1->second].begin(); iter2 != allele_kmer_counts[iter1->second].end(); ++iter2) {
                    combined_alelle_map[iter2->first] += iter2->second;
                }
                std::string genotype_name = iter1->first + "/" + iter1->second;
                genotype_kmer_counts[genotype_name] = combined_alelle_map;
            }
        }

        // Check if k-mer count profiles are unique amongst all alleles
        std::vector<std::pair<std::string, std::string>> identical_profiles;

        // k-mer comparison map
        std::map<std::string, kmer_comparison_map> kmer_comparisons;

        if (is_diploid == false) {
            for (auto iter = allele_pairs.begin(); iter != allele_pairs.end(); ++iter) {
                bool allele_vs_allele = compare_profiles(allele_kmer_counts[iter->first], allele_kmer_counts[iter->second]);
                if (allele_vs_allele == 1) {
                    identical_profiles.push_back(std::make_pair(iter->first, iter->second));
                }
                if (is_verbose == true && allele_vs_allele == 0) {
                    std::string compared_alleles = iter->first + "," + iter->second;
                    kmer_comparisons[compared_alleles] = kmer_differences(allele_kmer_counts[iter->first], allele_kmer_counts[iter->second]);
                }
            }
        }
        else if (is_diploid == true) {
            for (auto iter = genotype_pairs.begin(); iter != genotype_pairs.end(); ++iter) {
                bool genotype_vs_genotype = compare_profiles(genotype_kmer_counts[iter->first], genotype_kmer_counts[iter->second]);
                if (genotype_vs_genotype == 1) {
                    identical_profiles.push_back(std::make_pair(iter->first, iter->second));
                }
                if (is_verbose == true && genotype_vs_genotype == 0) {
                    std::string compared_genotypes = iter->first + "," + iter->second;
                    kmer_comparisons[compared_genotypes] = kmer_differences(genotype_kmer_counts[iter->first], genotype_kmer_counts[iter->second]);
                }
            }
        }

        if (is_verbose == false) {
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
        else if (is_verbose == true) {
            for (auto iter1 = kmer_comparisons.begin(); iter1 != kmer_comparisons.end(); ++iter1) {
                for (auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2) {
                    printf("%s,%s,%i\n", iter1->first.c_str(), iter2->first.c_str(), iter2->second);
                }

            }
        }

    }


    //
    // Finished
    //

    return 0;
}