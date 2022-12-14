## gbkc: genotyping by *k*-mer counts

Allele calling/genotyping approaches using *k*-mers for repetitive polymorphic loci.


### Commands for calling alleles/genotypes:
Each approach outputs a list of alleles/genotypes and the probabilistic scores for each, where the largest score indiciates the called allele/genotype.
#### count
Generate *k*-mer count proifles for each known allele and compare to *k*-mer count profile from a set of sequencing reads with unknown genotype.
#### distance
Get outermost *k*-mer pairs for each read pair in a set of sequencing reads with unknown genotype and locate the *k*-mer pairs in known alleles. Using the expected average fragment length for the reads, determine similarity to distance between *k*-mer pairs in known alleles.


### Additional useful commands:
#### check-profiles
Get information about *k*-mer count profiles for a list of known alleles. Print out the count profiles for each allele, determine which alleles share count profiles, or determine per-*k*-mer count differences between pairs of alleles.
#### get-genotypes
Get a list of all possible diploid genotypes given a list of known alleles.
