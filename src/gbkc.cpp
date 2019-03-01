#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <getopt.h>
#include <vector>
#include <map>
#include <zlib.h>
#include <stdint.h>
#include "kseq.h"

// Tell KSEQ what functions to use to open/read files
KSEQ_INIT(gzFile, gzread)

struct SequenceRecord
{
    std::string name;
    std::string sequence;
};

char complement(char c)
{
    switch(c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: 
            fprintf(stderr, "Error: unrecognized nucleotide %c\n", c);
            exit(1);
    }

}

std::string reverse_complement(const std::string& sequence)
{
    size_t l = sequence.size();
    std::string rc(l, 'N');
    for (size_t i = 0; i < l; ++i) {
        rc[i] = complement(sequence[l - i - 1]);
    }
    return rc;
}

std::string canonical_kmer(const std::string& kmer)
{
    std::string rc_kmer = reverse_complement(kmer);
    return kmer < rc_kmer ? kmer : rc_kmer;
}

std::map<std::string, size_t> count_kmers(const std::string& sequence, int k)
{
    std::map<std::string, size_t> out_map;
    for (size_t i = 0; i < sequence.size() - k + 1; ++i)
    {
        std::string canon_kmer = canonical_kmer(sequence.substr(i, k));
        out_map[canon_kmer] += 1;
    }
    return out_map;
}   

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

int main(int argc, char** argv) {

    int opt;

    std::string input_alleles_file;
    std::string input_reads_file;

    while ((opt = getopt(argc, argv, "a:r:")) != -1) {
        switch (opt) {
            case 'a':
                input_alleles_file = optarg;
                break;
            case 'r':
                input_reads_file = optarg;
                break;
            default: /* '?' */
                exit(EXIT_FAILURE);
        }
    }

    // TODO: Handle case where file names are null
    fprintf(stderr, "input reads: %s\n", input_reads_file.c_str());
    fprintf(stderr, "input alleles: %s\n", input_alleles_file.c_str());

    std::vector<SequenceRecord> alleles = read_sequences_from_file(input_alleles_file);

    for (size_t i = 0; i < alleles.size(); ++i) {
        printf("name: %s seq: %s\n", alleles[i].name.c_str(), alleles[i].sequence.substr(0, 10).c_str());

    }

    std::vector<SequenceRecord> reads = read_sequences_from_file(input_reads_file);

    for (size_t i = 0; i < reads.size(); ++i) {
        printf("name: %s seq: %s\n", reads[i].name.c_str(), reads[i].sequence.substr(0, 10).c_str());

    }

    // Count kmers //

    int k = 21;

    std::map<std::string, size_t> allele_kmer_map = count_kmers(alleles[0].sequence, k);
    //std::map<std::string, size_t>::iterator iter = allele_kmer_map.begin();
    for(auto iter = allele_kmer_map.begin(); iter != allele_kmer_map.end(); ++iter) {
        printf("kmer: %s count: %zu\n", iter->first.c_str(), iter->second);
    }



}
