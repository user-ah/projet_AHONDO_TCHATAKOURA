// haplotype_generator.h
#ifndef HAPLOTYPE_GENERATOR_H
#define HAPLOTYPE_GENERATOR_H

#include <vector>
#include <cstdlib>  // pour rand()

/**
 * Génère deux haplotypes à partir d'un génotype donné.
 * 
 * @param genotype Vecteur d'entiers représentant le génotype.
 * @param haplotype1 Référence pour stocker le premier haplotype.
 * @param haplotype2 Référence pour stocker le deuxième haplotype.
 */
void generateHaplotypes(const std::vector<int> &genotype, std::vector<int> &haplotype1, std::vector<int> &haplotype2) {
    for (size_t j = 0; j < genotype.size(); ++j) {
        if (genotype[j] == 0) {
            haplotype1.push_back(0);
            haplotype2.push_back(0);
        } else if (genotype[j] == 1) {
            haplotype1.push_back(1);
            haplotype2.push_back(1);
        } else if (genotype[j] == 2) {
            // Pour le génotype "2", on attribue "1 0" ou "0 1" de manière aléatoire
            if (rand() % 2 == 0) {
                haplotype1.push_back(1);
                haplotype2.push_back(0);
            } else {
                haplotype1.push_back(0);
                haplotype2.push_back(1);
            }
        }
    }
}

#endif // HAPLOTYPE_GENERATOR_H