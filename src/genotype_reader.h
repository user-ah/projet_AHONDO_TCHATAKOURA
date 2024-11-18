// genotype_reader.h
#ifndef GENOTYPE_READER_H
#define GENOTYPE_READER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

/**
 * Lit le fichier de génotypes et détermine le nombre d'individus (n_ind) et le nombre de loci (n_loci).
 *
 * @param filename Le nom du fichier de génotypes.
 * @param genotypes Un vecteur de vecteurs pour stocker les génotypes.
 * @param n_ind Référence pour stocker le nombre d'individus.
 * @param n_loci Référence pour stocker le nombre de loci.
 * @return bool Vrai si la lecture a réussi, faux sinon.
 */
bool readGenotypeFile(const std::string &filename, std::vector<std::vector<int>> &genotypes, int &n_ind, int &n_loci) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier de génotypes " << filename << std::endl;
        return false;
    }
    
    std::string line;
    n_ind = 0;
    while (getline(file, line)) {
        std::vector<int> genotype;
        for (char ch : line) {
            if (ch != ' ' && ch != ',') {  // Ignorer les espaces et les virgules
                genotype.push_back(ch - '0');  // Convertir le caractère en entier
            }
        }
        if (n_ind == 0) {
            n_loci = genotype.size();  // Définir n_loci en fonction de la première ligne
        } else if (genotype.size() != n_loci) {
            std::cerr << "Erreur : Nombre de loci incohérent dans le fichier de génotypes" << std::endl;
            return false;
        }
        genotypes.push_back(genotype);
        ++n_ind;
    }
    file.close();
    return true;
}

#endif // GENOTYPE_READER_H
