#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <getopt.h>
#include <vector>
#include <map>
#include <zlib.h>
#include <stdint.h>
#include "kseq.h"
#include <sstream>
#include <iostream>
#include <unordered_set>

//
// Define functions
//

// Tell KSEQ (from Heng Li) what functions to use to open/read files
KSEQ_INIT(gzFile, gzread)

struct SequenceRecord
{
    std::string name;
    std::string sequence;
};

// Define kmer count profile structure
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

// Obtain all (canonical) kmers and their counts from a sequence
kmer_count_map count_kmers(const std::string& sequence, size_t k)
{
    kmer_count_map out_map;
    for (size_t i = 0; i < sequence.size() - k + 1; ++i)
    {
        std::string canon_kmer = canonical_kmer(sequence.substr(i, k));
        out_map[canon_kmer] += 1;
    }
    return out_map;
}   

// Read input files
std::vector<SequenceRecord> read_sequences_from_file(const std::string& input_filename)
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
    
    std::vector<SequenceRecord> out_sequences;
    int ret = 0;
    kseq_t* seq = kseq_init(gz_read_fp);
    while((ret = kseq_read(seq)) >= 0) {
        SequenceRecord sr = { seq->name.s, seq->seq.s };
        out_sequences.push_back(sr);
    }    
    kseq_destroy(seq);   

    return out_sequences;
}

//
// Program commands
//

int main(int argc, char** argv) {

    // Read command line arguments
    std::string input_alleles_file;
    std::string input_reads_file;
    size_t input_k = 21;

    for (char c; (c = getopt_long(argc, argv, "a:r:k:", NULL, NULL)) != -1;) {
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
            default:
                exit(EXIT_FAILURE);
        }
    }

    // Read files
    // TODO: Handle case where file names are null
    fprintf(stderr, "input reads: %s\n", input_reads_file.c_str());
    fprintf(stderr, "input alleles: %s\n", input_alleles_file.c_str());
    fprintf(stderr, "input k value: %zu\n", input_k);

    std::vector<SequenceRecord> alleles = read_sequences_from_file(input_alleles_file);
    // for (size_t i = 0; i < alleles.size(); ++i) {
    //     printf("name: %s seq: %s\n", alleles[i].name.c_str(), alleles[i].sequence.substr(0, 10).c_str());
    // }

    std::vector<SequenceRecord> reads = read_sequences_from_file(input_reads_file);
    // for (size_t i = 0; i < reads.size(); ++i) {
    //     printf("name: %s seq: %s\n", reads[i].name.c_str(), reads[i].sequence.substr(0, 10).c_str());
    // }

    fprintf(stderr, "number of alleles: %zu\n", alleles.size());
    // for (size_t a = 0; a < alleles.size(); ++a) {
    //     fprintf(stderr, "allele name: %s\n", alleles[a].name.c_str());
    // }


    //
    // Count k-mers in alleles
    //

    // All kmers in each allele
    std::map<std::string, kmer_count_map> allele_kmer_counts; 

    // Union of kmers fromm all alleles
    std::unordered_set<std::string> allele_kmers;

    // Iterate over each allele
    for (size_t a = 0; a < alleles.size(); ++a) {
        std::map<std::string, size_t> single_allele_kmer_counts = count_kmers(alleles[a].sequence, input_k);
        allele_kmer_counts[alleles[a].name.c_str()] = single_allele_kmer_counts;
        for (auto iter = single_allele_kmer_counts.begin(); iter != single_allele_kmer_counts.end(); ++iter) {
            allele_kmers.insert(iter->first);
        }

    }

    // Print map contents for debugging
    // for (auto iter = allele_kmer_counts.begin(); iter != allele_kmer_counts.end(); ++iter) {
    //     printf("allele: %s\n", iter->first.c_str());
    //     std::map<std::string, size_t> &single_allele_kmer_counts = iter->second;
    //     for (auto iter2 = single_allele_kmer_counts.begin(); iter2 != single_allele_kmer_counts.end(); ++iter2) {
    //         printf("kmer: %s, count: %zu\n", iter2->first.c_str(), iter2->second);
    //     }
    // }

    // for (auto iter = allele_kmers.begin(); iter != allele_kmers.end(); ++iter) {
    //     std::cout << *iter << '\n';
    // }

    

    //
    // Count k-mers in reads
    //

    // k-mers and counts for each read
    std::map<std::string, kmer_count_map> each_read_kmer_counts;

    // k-mers and counts across all reads
    kmer_count_map all_reads_kmer_counts;

    // Union of k-mers across all reads
    std::unordered_set<std::string> reads_kmers;

    // Iterate over each read
    for (size_t r = 0; r < reads.size(); ++r) {
        std::map<std::string, size_t> single_read_kmer_counts = count_kmers(reads[r].sequence, input_k);
        each_read_kmer_counts[reads[r].name.c_str()] = single_read_kmer_counts;
        for (auto iter = single_read_kmer_counts.begin(); iter != single_read_kmer_counts.end(); ++iter) {
            reads_kmers.insert(iter->first);
            all_reads_kmer_counts[iter->first]  += iter->second;
        }
    }

    // Print map contents for debugging
    // for (auto iter = each_read_kmer_counts.begin(); iter != each_read_kmer_counts.end(); ++iter) {
    //     printf("read: %s\n", iter->first.c_str());
    //     std::map<std::string, size_t> &single_read_kmer_counts = iter->second;
    //     for (auto iter2 = single_read_kmer_counts.begin(); iter2 != single_read_kmer_counts.end(); ++iter2) {
    //         printf("kmer: %s, count: %zu\n", iter2->first.c_str(), iter2->second);
    //     }
    // }

    // for (auto iter = all_reads_kmer_counts.begin(); iter != all_reads_kmer_counts.end(); ++iter) {
    //     printf("kmer: %s, count: %zu\n", iter->first.c_str(), iter->second);
    // }    





}
