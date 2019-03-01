#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <getopt.h>
#include <vector>
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


}
