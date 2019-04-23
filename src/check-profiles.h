// check-profiles - Compare k-mer count profiles for all known alleles to check if they are unique

#ifndef CHECK_PROFILES_H
#define CHECK_PROFILES_H

#include <iostream>
#include <map>
#include <string>
#include <vector>


//
// Type defs
//

// Define k-mer count profile structure
typedef std::map<std::string, size_t> kmer_count_map;

// Define sequence record structre
struct sequence_record
{
    std::string name;
    std::string sequence;
};


//
// Functions
//

int checkprofilesMain(int argc, char** argv);

char complement(char c);
std::string reverse_complement(const std::string& sequence);
std::string canonical_kmer(const std::string& kmer);
kmer_count_map count_kmers(const std::string& sequence, size_t k);
std::vector<sequence_record> read_sequences_from_file(const std::string& input_filename);


#endif