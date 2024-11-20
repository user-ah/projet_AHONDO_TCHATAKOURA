// initialisation_frequences.h
#ifndef INITIALISATION_FREQUENCES_H
#define INITIALISATION_FREQUENCES_H

#include <vector>
#include <iostream>

/**
 * Initialise les fréquences des haplotypes avec une probabilité uniforme.
 *
 * @param haplotype_frequencies Un vecteur pour stocker les fréquences des haplotypes.
 * @param nb_haplotypes Le nombre total d'haplotypes possibles.
 */
inline void initialisation_frequences_haplotypes(std::vector<double> &haplotype_frequencies, int nb_haplotypes) {
    double initial_freq = 1.0 / nb_haplotypes;  // Fréquence uniforme pour chaque haplotype
    haplotype_frequencies.assign(nb_haplotypes, initial_freq);
    std::cout << "Fréquences des haplotypes initialisées à " << initial_freq << " pour chaque haplotype." << std::endl;
}

#endif // INITIALISATION_FREQUENCES_H
