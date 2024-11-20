#ifndef DISTANCE_CALCULATION_H
#define DISTANCE_CALCULATION_H

#include <vector>
using namespace std;

/**
 * Calcule la distance de Hamming entre deux haplotypes.
 *
 * @param haplotype1 Le premier haplotype (vecteur de 0 et 1).
 * @param haplotype2 Le second haplotype (vecteur de 0 et 1).
 * @return La distance de Hamming entre les deux haplotypes.
 */
inline int calculateHammingDistance(const vector<int>& haplotype1, const vector<int>& haplotype2) {
    int distance = 0;
    for (size_t i = 0; i < haplotype1.size(); ++i) {
        if (haplotype1[i] != haplotype2[i]) {
            ++distance;
        }
    }
    return distance;
}

/**
 * Calcule la distance de Hamming entre les haplotypes de référence et les haplotypes inférés.
 *
 * @param ground_truth Un vecteur de vecteurs représentant les haplotypes de référence.
 * @param inferred Un vecteur de vecteurs représentant les haplotypes inférés.
 * @return La distance totale de Hamming entre les haplotypes de référence et les haplotypes inférés.
 */
inline int calculateHaplotypeDistance(const vector<vector<int>>& ground_truth,
                               const vector<vector<int>>& inferred) {
    int total_distance = 0;
    // Assurez-vous que les tailles des deux vecteurs sont les mêmes
    for (size_t i = 0; i < ground_truth.size(); ++i) {
        total_distance += calculateHammingDistance(ground_truth[i], inferred[i]);
    }
    return total_distance;
}

#endif // DISTANCE_CALCULATION_H
