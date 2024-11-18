// haplotype_counter.h
#ifndef HAPLOTYPE_COUNTER_H
#define HAPLOTYPE_COUNTER_H

#include <vector>
#include <unordered_set>
#include <string>

/**
 * Compte le nombre d'haplotypes distincts dans une liste d'haplotypes.
 *
 * @param haplotypes Un vecteur de vecteurs représentant les haplotypes.
 * @return Le nombre d'haplotypes distincts.
 */
int countDistinctHaplotypes(const std::vector<std::vector<int>> &haplotypes) {
    std::unordered_set<std::string> distinct_haplotypes;
    
    for (const auto& haplotype : haplotypes) {
        std::string haplotype_str;
        for (int allele : haplotype) {
            haplotype_str += std::to_string(allele);
        }
        distinct_haplotypes.insert(haplotype_str);
    }
    
    return distinct_haplotypes.size();
}

#endif // HAPLOTYPE_COUNTER_H
