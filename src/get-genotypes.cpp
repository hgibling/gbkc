// get-genotypes - Get all possible diploid genotypes given all known alleles

#include <getopt.h>
#include <iostream>
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
using std::istringstream;
using std::pair;
using std::set;
using std::string;
using std::vector;


//
// Help message
//

static const char *GET_GENOTYPES_USAGE_MESSAGE =
"Get all possible diploid genotypes given all known alleles/haplotypes.\n\n"
"Usage: gbkc get-genotypes -a <file>\n\n"
"Commands:\n"
"       -a       multi-fasta file of alleles/haplotypes of interest\n"
"       -s       separator between alleles of a genotype (default: /)\n"
"       -o       output file name (default: results-get-genotypes.txt)\n";


//
// Program commands
//

int getgenotypesMain(int argc, char** argv) {

    if(argc <= 1) {
        cout << GET_GENOTYPES_USAGE_MESSAGE;
        return 0;
    };


    //
    // Read command line arguments
    //

    string input_alleles_file;
    string separator = "/";
    string output_name = "results-get-genotypes.txt";

    for (char c; (c = getopt_long(argc, argv, "a:s:o:", NULL, NULL)) != -1;) {
        istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'a': arg >> input_alleles_file; break;
            case 's': arg >> separator; break;
            case 'o': arg >> output_name; break;
            default: exit(EXIT_FAILURE);
        }
    }

    // Check for parameter arguments
    if (input_alleles_file.empty()) {
        fprintf(stderr, "No file for allele sequences. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

    if (separator.empty()) {
        fprintf(stderr, "Please use some sort of separator for the two alleles that comprise the genotype\n");
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

    // Get list of all possible allele combinations
    vector<pair<string, string>> genotypes;
    genotypes = pairwise_comparisons(allele_names, true);
    set<string> genotype_names;
    for (auto iter = genotypes.begin(); iter != genotypes.end(); ++iter) {
        string name = iter->first + separator + iter->second;
        genotype_names.insert(name);
    }

    // Print list of all possible allele combinations/genotypes
    FILE * output;
    output = fopen(output_name.c_str(), "w");

    for (auto iter = genotype_names.begin(); iter != genotype_names.end(); ++iter) {
        fprintf(output, "%s\n", iter->c_str());
    }

    fclose(output);


    //
    // Finished
    //

    return 0;
}
