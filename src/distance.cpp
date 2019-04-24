// distance - Score reads based on comparison of k-mer pair distances to those for known alleles

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
#include "distance.h"
#include "kseq.h"


//
// Define functions
//

kmer_position_map kmer_positions(const std::string& sequence, size_t k)
{
	kmer_position_map out_map;
    for (size_t i = 0; i < sequence.size() - k + 1; ++i) {
        out_map.insert(std::pair<std::string, size_t>(sequence.substr(i, k), i));
    }
    return out_map;
}


//
// Help message
//

static const char *DISTANCE_USAGE_MESSAGE =
"Score reads based on comparison of k-mer pair distances to those for known alleles.\n\n"
"Usage: gbkc distance [options]\n\n"
"Commands:\n"
"		-a       multi-fasta file of alleles/haplotypes of interest\n"
"		-1       multi-fasta/q file of sequencing reads to score (first in pair)\n"
"		-2       multi-fasta/q file of sequencing reads to score (second in pair)\n"
"		-k       size of k-mers to use\n"
"		-l       read length\n"
"		-e       sequencing error rate\n"
"		-c       sequencing coverage\n"
"		-f       mean fragment length\n"
"		-s       standard deviation of fragment length\n"
"		-o       output file name (default: results.csv)\n";


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
    std::string output_name = "distance-results.csv";

    for (char c; (c = getopt_long(argc, argv, "a:1:2:k:l:e:c:f:s:o:", NULL, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'a':
                arg >> input_alleles_file;
                break;
            case '1':
                arg >> input_reads_file1;
                break;
            case '2':
                arg >> input_reads_file2;
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
            case 'f':
                arg >> fragment_length;
                break;
            case 's':
                arg >> fragment_stdev;
                break;
            case 'o':
                arg >> output_name;
                break;
            default:
                exit(EXIT_FAILURE);
        }
    }


    // Check for parameter arguments
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


    //
    // Print handy information
    //

    fprintf(stderr, "input reads: %s and %s\n", input_reads_file1.c_str(), input_reads_file2.c_str());
    fprintf(stderr, "input alleles: %s\n", input_alleles_file.c_str());
    fprintf(stderr, "number of alleles: %zu\n", allele_names.size());
    fprintf(stderr, "input k value: %zu\n", input_k);
    fprintf(stderr, "input coverage: %f X, sequencing error: %f %%\n", coverage, sequencing_error);
    fprintf(stderr, "input mean fragment length: %f, standard deviation: %f\n", fragment_length, fragment_stdev);


    //
    // Get k-mer distances for alleles
    //

    std::map<std::string, kmer_position_map> allele_positions;
    for (size_t a = 0; a < alleles.size(); ++a) {
    	kmer_position_map single_allele_positions = kmer_positions(alleles[a].sequence, input_k);
    	allele_positions[alleles[a].name] = single_allele_positions;
    }

    // for (auto iter = allele_positions.begin(); iter != allele_positions.end(); ++iter) {
    // 	fprintf(stderr, "%s:\n", iter->first.c_str());
    // 	kmer_position_map &single_allele_positions = iter->second;
    // 	for (auto iter2 = single_allele_positions.begin(); iter2 != single_allele_positions.end(); ++iter2) {
    // 		fprintf(stderr, "%s, %zu\n", iter2->first.c_str(), iter2->second);
    // 	}
    // }


    //
    // Finished
    //

	return 0;
}