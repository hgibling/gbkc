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

using std::cout;
using std::equal;
using std::istringstream;
using std::pair;
using std::make_pair;
using std::map;
using std::set;
using std::string;
using std::vector;


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
	case 'N': return 'N';
        default:
            fprintf(stderr, "Error: unrecognized nucleotide %c\n", nucleotide);
            exit(EXIT_FAILURE);
    }
}

// Get reverse complement of a sequence
string reverse_complement(const string& sequence)
{
    size_t length = sequence.size();
    string out_reverse_complement(length, 'X');
    for (size_t i = 0; i < length; ++i) {
        out_reverse_complement[i] = complement(sequence[length - i - 1]);
    }
    return out_reverse_complement;
}

// Get the lexicographically lowest k-mer between it and its reverse complement
string canonical_kmer(const string& kmer)
{
    string reverse_complement_kmer = reverse_complement(kmer);
    return kmer < reverse_complement_kmer ? kmer : reverse_complement_kmer;
}

// Obtain all (canonical) k-mers and their counts from a sequence
kmer_count_map count_kmers(const string& sequence, const size_t k)
{
    kmer_count_map out_map;
    char N_char = 'N';
    for (size_t i = 0; i < sequence.size() - k + 1; ++i) {
        string canon_kmer = canonical_kmer(sequence.substr(i, k));
        // Only keep k-mers without Ns
        if (canon_kmer.find(N_char) == string::npos) {
            out_map[canon_kmer] += 1;
        }
    }
    return out_map;
}

// Read input files
vector<sequence_record> read_sequences_from_file(const string& input_filename)
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

    vector<sequence_record> out_sequences;
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
        && equal(first_map.begin(), first_map.end(),
                      second_map.begin());
}

// Generate count differences between two k-mer count profiles
kmer_comparison_map kmer_differences(const kmer_count_map& first_map, const kmer_count_map& second_map)
{
    kmer_comparison_map kmer_differences_out;
    set<string> all_kmers;
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
vector<pair<string, string>> pairwise_comparisons(const set<string>& alleles, bool same)
{
    vector<pair<string, string>> pairs;
    vector<string> alleles_vector(alleles.begin(), alleles.end());
    for (size_t i = 0; i < alleles_vector.size(); ++i) {
        for (size_t j = 0; j < alleles_vector.size(); ++j) {
            if ((same == true) and (i <= j)) {
                pairs.push_back(make_pair(alleles_vector[i], alleles_vector[j]));
            }
            else if (i < j) {
                pairs.push_back(make_pair(alleles_vector[i], alleles_vector[j]));
            }
        }
    }
    return pairs;
}

// Split a diploid genotype into a vector of its alleles
vector<string> genotype_split(const string genotype, const char delimiter)
{
    size_t position = genotype.find(delimiter);
    string allele1 = genotype.substr(0, position);
    string allele2 = genotype.substr(position + 1, genotype.size());
    vector<string> out = {allele1, allele2};
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
"       -v       verbose printout of comparisons when profiles are not identical\n"
"       -p       print individual k-mer count profiles instead of comparisons\n";


//
// Program commands
//

int checkprofilesMain(int argc, char** argv) {

    if(argc <= 1) {
        cout << CHECK_PROFILES_USAGE_MESSAGE;
        return 0;
    };


    //
    // Read command line arguments
    //

    string input_alleles_file;
    size_t lower_range = 3;
    size_t upper_range = 0;
    bool is_diploid = false;
    bool is_verbose = false;
    bool is_profiles = false;


    for (char c; (c = getopt_long(argc, argv, "a:l:u:dvp", NULL, NULL)) != -1;) {
        istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'a': arg >> input_alleles_file; break;
            case 'l': arg >> lower_range; break;
            case 'u': arg >> upper_range; break;
            case 'd': is_diploid = true; break;
            case 'v': is_verbose = true; break;
            case 'p': is_profiles = true; break;
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
    vector<sequence_record> alleles = read_sequences_from_file(input_alleles_file);

    // Get list of allele names
    set<string> allele_names;
    for (size_t i = 0; i < alleles.size(); ++i) {
        allele_names.insert(alleles[i].name.c_str());
    }


    //
    // Compare k-mer count profiles
    //

    // // Get list of allele comparisons to make
    vector<pair<string, string>> allele_pairs;
    vector<pair<string, string>> genotype_pairs;
    vector<pair<string, string>> genotypes;

    if (!is_diploid) {
        allele_pairs = pairwise_comparisons(allele_names, false);
    }

    else if (is_diploid) {
        genotypes = pairwise_comparisons(allele_names, true);
        set<string> genotype_names;
        for (auto iter = genotypes.begin(); iter != genotypes.end(); ++iter) {
            string name = iter->first + "/" + iter->second;
            genotype_names.insert(name);
        }
        genotype_pairs = pairwise_comparisons(genotype_names, false);
    }


    // Iterate over all possible values of k given length l
    for (size_t k = lower_range; k < upper_range + 1; ++k) {
        fprintf(stderr, "Testing k = %zu... \n", k);

        // All k-mers in each allele or genotype
        map<string, kmer_count_map> allele_kmer_counts;
        map<string, kmer_count_map> genotype_kmer_counts;

        // Get count profiles for each allele
        for (size_t a = 0; a < alleles.size(); ++a) {
            kmer_count_map single_allele_kmer_counts = count_kmers(alleles[a].sequence, k);
            allele_kmer_counts[alleles[a].name] = single_allele_kmer_counts;
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

        if (is_profiles) {

            if (is_diploid) {
                for (auto iter1 = genotype_kmer_counts.begin(); iter1 != genotype_kmer_counts.end(); ++iter1) {
                    for (auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2) {
                        printf("%zu,%s,%s,%zu\n", k, iter1->first.c_str(), iter2->first.c_str(), iter2->second);
                    }
                }
            }

            else if (!is_diploid) {
               for (auto iter1 = allele_kmer_counts.begin(); iter1 != allele_kmer_counts.end(); ++iter1) {
                    for (auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2) {
                        printf("%zu,%s,%s,%zu\n", k, iter1->first.c_str(), iter2->first.c_str(), iter2->second);
                    }
                }
            }
        }

        else if (!is_profiles) {

            printf("Testing k = %zu... \n", k);

            // Check if k-mer count profiles are unique amongst all alleles
            vector<pair<string, string>> identical_profiles;

            // k-mer comparison map
            map<string, kmer_comparison_map> kmer_comparisons;

            if (!is_diploid) {
                for (auto iter = allele_pairs.begin(); iter != allele_pairs.end(); ++iter) {
                    bool allele_vs_allele = compare_profiles(allele_kmer_counts[iter->first], allele_kmer_counts[iter->second]);
                    if (allele_vs_allele == 1) {
                        identical_profiles.push_back(make_pair(iter->first, iter->second));
                    }
                    if (is_verbose == true && allele_vs_allele == 0) {
                        string compared_alleles = iter->first + "," + iter->second;
                        kmer_comparisons[compared_alleles] = kmer_differences(allele_kmer_counts[iter->first], allele_kmer_counts[iter->second]);
                    }
                }
            }

            else if (is_diploid) {
                for (auto iter = genotype_pairs.begin(); iter != genotype_pairs.end(); ++iter) {
                    bool genotype_vs_genotype = compare_profiles(genotype_kmer_counts[iter->first], genotype_kmer_counts[iter->second]);
                    if (genotype_vs_genotype == 1) {
                        identical_profiles.push_back(make_pair(iter->first, iter->second));
                    }
                    if (is_verbose == true && genotype_vs_genotype == 0) {
                        string compared_genotypes = iter->first + "," + iter->second;
                        kmer_comparisons[compared_genotypes] = kmer_differences(genotype_kmer_counts[iter->first], genotype_kmer_counts[iter->second]);
                    }
                }
            }

            if (!is_verbose) {
                if (identical_profiles.size() > 0) {
                    printf("Some alleles had identical k-mer count profiles:\n");
                    for (auto iter = identical_profiles.begin(); iter != identical_profiles.end(); ++iter) {
                        printf("%s and %s at k=%zu\n", iter->first.c_str(), iter->second.c_str(), k);
                    }
                }
                else {
                    printf("All k-mer count profiles are unique\n");
                }
            }

            else if (is_verbose) {
                // Output is second profile - first profile
                for (auto iter1 = kmer_comparisons.begin(); iter1 != kmer_comparisons.end(); ++iter1) {
                    for (auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2) {
                        printf("%s,%s,%i\n", iter1->first.c_str(), iter2->first.c_str(), iter2->second);
                    }

                }
            }
        }

    }


    //
    // Finished
    //

    return 0;
}
