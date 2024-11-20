// calcul_proba_genotypes.h
#ifndef CALCUL_PROBA_GENOTYPES_H
#define CALCUL_PROBA_GENOTYPES_H

#include <vector>
using namespace std;

/**
 * Calcule la probabilité de chaque génotype en fonction des fréquences des paires d'haplotypes possibles.
 *
 * @param genotypes Un vecteur de génotypes.
 * @param haplotype_frequencies Un vecteur de fréquences des haplotypes.
 * @param genotype_probabilities Un vecteur pour stocker les probabilités calculées pour chaque génotype.
 */
inline void calcul_proba_genotypes(const std::vector<std::vector<int>> &genotypes, 
                            const std::vector<double> &haplotype_frequencies, 
                            std::vector<double> &genotype_probabilities) {
    genotype_probabilities.clear();
    
    for (const auto &genotype : genotypes) {
        double proba = 0.0;

        // Itération sur toutes les paires d'haplotypes explicatives pour le génotype
        for (size_t h1 = 0; h1 < haplotype_frequencies.size(); ++h1) {
            for (size_t h2 = h1; h2 < haplotype_frequencies.size(); ++h2) {
                double p1 = haplotype_frequencies[h1];
                double p2 = haplotype_frequencies[h2];
                double contrib_proba;

                if (h1 == h2) {
                    contrib_proba = p1 * p1;  // Cas où h1 et h2 sont identiques
                } else {
                    contrib_proba = 2 * p1 * p2;  // Cas où h1 et h2 sont différents
                }

                proba += contrib_proba;
            }
        }
        genotype_probabilities.push_back(proba);
    }
}

#endif // CALCUL_PROBA_GENOTYPES_H
