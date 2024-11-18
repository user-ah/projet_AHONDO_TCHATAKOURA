// estimation_esperance.h
#ifndef ESTIMATION_ESPERANCE_H
#define ESTIMATION_ESPERANCE_H

#include <vector>
#include <cmath>

/**
 * Calcule le log-vraisemblance pour l'étape actuelle de l'algorithme EM.
 *
 * Préconditions :
 * - Les fréquences des haplotypes sont connues pour l'étape actuelle (étape i).
 * - Le nombre d'individus possédant chaque génotype (Ng) est connu.
 *
 * Postconditions :
 * - Retourne la valeur du log-vraisemblance (loglikelihood).
 *
 * @param genotypes Un vecteur de génotypes.
 * @param haplotype_frequencies Les fréquences actuelles des haplotypes.
 * @return double La valeur du log-vraisemblance calculée.
 */
double estimation_esperance(const std::vector<std::vector<int>> &genotypes, 
                            const std::vector<double> &haplotype_frequencies) {
    double loglikelihood = 0.0;

    for (const auto &genotype : genotypes) {
        double proba = 0.0;

        for (size_t h1 = 0; h1 < haplotype_frequencies.size(); ++h1) {
            for (size_t h2 = h1; h2 < haplotype_frequencies.size(); ++h2) {
                double p1 = haplotype_frequencies[h1];
                double p2 = haplotype_frequencies[h2];
                double contrib_proba = (h1 == h2) ? (p1 * p1) : (2 * p1 * p2);

                proba += contrib_proba;
            }
        }

        loglikelihood += std::log(proba);
    }

    return loglikelihood;
}

#endif // ESTIMATION_ESPERANCE_H
