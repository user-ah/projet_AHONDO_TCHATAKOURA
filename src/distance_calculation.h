// distance_calculation.h
#ifndef DISTANCE_CALCULATION_H
#define DISTANCE_CALCULATION_H

#include <vector>

/**
 * Calcule la distance entre les haplotypes inférés et les haplotypes de référence.
 *
 * @param ground_truth Un vecteur de vecteurs représentant les haplotypes de référence.
 * @param inferred Un vecteur de vecteurs représentant les haplotypes inférés.
 * @return La distance totale entre les haplotypes de référence et les haplotypes inférés.
 */
int calculateHaplotypeDistance(const std::vector<std::vector<int>> &ground_truth,
                               const std::vector<std::vector<int>> &inferred) {
    int distance = 0;
    for (size_t i = 0; i < ground_truth.size(); ++i) {
        for (size_t j = 0; j < ground_truth[i].size(); ++j) {
            if (ground_truth[i][j] != inferred[i][j]) {
                ++distance;
            }
        }
    }
    return distance;
}

#endif // DISTANCE_CALCULATION_H
