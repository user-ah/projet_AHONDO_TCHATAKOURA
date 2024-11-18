// maximisation.h
#ifndef MAXIMISATION_H
#define MAXIMISATION_H

#include <vector>

/**
 * Met à jour les fréquences des haplotypes en fonction des génotypes observés et de leurs probabilités.
 *
 * Préconditions :
 * - Chaque haplotype est associé à la liste des génotypes pour lesquels il peut être explicatif.
 * - Les fréquences des haplotypes (freq_prec) et les probabilités des génotypes (proba_prec) sont connues à l'étape i-1.
 * - Le nombre total d'individus (nb_ind) et le nombre d'individus possédant chaque génotype (Ng) sont connus.
 *
 * Postconditions :
 * - Les fréquences des haplotypes sont mises à jour pour l'étape actuelle (étape i).
 *
 * @param genotypes Un vecteur de génotypes.
 * @param haplotype_frequencies Référence pour les fréquences actuelles des haplotypes.
 * @param genotype_probabilities Les probabilités des génotypes pour l'étape précédente.
 * @param n_ind Le nombre total d'individus.
 */
void maximisation(const std::vector<std::vector<int>> &genotypes, 
                  std::vector<double> &haplotype_frequencies, 
                  const std::vector<double> &genotype_probabilities, 
                  int n_ind) {
    std::vector<double> updated_frequencies(haplotype_frequencies.size(), 0.0);

    for (size_t h1 = 0; h1 < haplotype_frequencies.size(); ++h1) {
        double freq = 0.0;
        
        for (size_t g = 0; g < genotypes.size(); ++g) {
            double Ng = 1.0;  // Supposons un individu par génotype pour simplification

            // Calculer la contribution de l'haplotype à la probabilité du génotype
            if (h1 == h1) {
                freq += ((2 * haplotype_frequencies[h1] * haplotype_frequencies[h1]) / genotype_probabilities[g]) * (Ng / n_ind);
            } else {
                for (size_t h2 = h1; h2 < haplotype_frequencies.size(); ++h2) {
                    freq += (2 * haplotype_frequencies[h1] * haplotype_frequencies[h2] / genotype_probabilities[g]) * (Ng / n_ind);
                }
            }
        }

        updated_frequencies[h1] = freq / 2.0;
    }

    haplotype_frequencies = updated_frequencies;
}

#endif // MAXIMISATION_H
