// distance - Score reads based on comparison of k-mer pair distances to those for known alleles

#ifndef DISTANCE_H
#define DISTANCE_H

#include <iostream>
#include <string>

#include "check-profiles.h"
#include "count.h"


//
// Type defs
//

// Define k-mer position structure
typedef std::multimap<std::string, size_t> kmer_position_map;
typedef std::multimap<std::string, std::pair<std::string, std::string>> kmer_pairs_map;

int distanceMain(int argc, char** argv);


#endif