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


//
// Help message
//

static const char *GET_GENOTYPES_USAGE_MESSAGE =
"Get all possible diploid genotypes given all known alleles/haplotypes.\n\n"
"Usage: gbkc get-genotypes -a <file>\n\n"
"Commands:\n"
"       -a       multi-fasta file of alleles/haplotypes of interest\n"
"       -s       single-character separator (default: ,)\n";


//
// Program commands
//

int getgenotypesMain(int argc, char** argv) {

    if(argc <= 1) {
        std::cout << GET_GENOTYPES_USAGE_MESSAGE;
        return 0;
    };


    //
    // Read command line arguments
    //

    std::string input_alleles_file;
    char separator = ',';

    for (char c; (c = getopt_long(argc, argv, "a:", NULL, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'a': arg >> input_alleles_file; break;
            case 's': arg >> separator; break;
            default: exit(EXIT_FAILURE);
        }
    }

    // Check for parameter arguments
    if (input_alleles_file.empty()) {
        fprintf(stderr, "No file for allele sequences. Check parameters.\n");
        exit(EXIT_FAILURE);
    }

//    if (separator == '') {
//        fprintf(stderr, "Please use some sort of separator for the two alleles that comprise the genotype\n");
//        exit(EXIT_FAILURE);
//    }


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

    // Get list of all possible allele combinations
    std::vector<std::pair<std::string, std::string>> genotypes;
    genotypes = pairwise_comparisons(allele_names, true);
    std::set<std::string> genotype_names;
    for (auto iter = genotypes.begin(); iter != genotypes.end(); ++iter) {
        std::string name = iter->first + separator + iter->second;
        genotype_names.insert(name);
    }

    // Print list of all possible allele combinations/genotypes
    for (auto iter = genotype_names.begin(); iter != genotype_names.end(); ++iter) {
        printf("%s\n", iter->c_str());
    }


    //
    // Finished
    //

    return 0;
}
