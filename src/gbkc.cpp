#include <iostream>
#include <string>

#include "check-profiles.h"
#include "count.h"
#include "distance.h"


static const char *GBKC_USAGE_MESSAGE =
"Program: gbkc\n"
"Usage: gbkc <command> [options]\n\n"
"Commands:\n"
"          check-profiles           Compare k-mer count profiles for all known alleles to check if they are unique\n"
"          count                    Score reads based on comparison of k-mer count profiles to those for known alleles\n"
"          distance                 Score reads based on comparison of k-mer pair distances to those for known alleles\n";

int main(int argc, char** argv){


    //
    // Read command line arguments
    //
    if(argc <= 1) {
        std::cout << GBKC_USAGE_MESSAGE;
        return 0;
    }
    else {
        std::string command(argv[1]);
        if(command == "help" || command == "--help")
        {
            std::cout << GBKC_USAGE_MESSAGE;
            return 0;
        }

        if(command == "check-profiles") {
            checkprofilesMain(argc - 1, argv + 1);
        }
        else if(command == "count") {
            countMain(argc - 1, argv + 1);
        }
        else if(command == "distance") {
            distanceMain(argc - 1, argv + 1);
        }

        else {
            std::cerr << "Unrecognized command: " << command << "\n";
            return 1;
        }
    }

    //
    // Finished
    //
    return 0;
}
