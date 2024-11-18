// parametre.h
#ifndef PARAMETRE_H
#define PARAMETRE_H

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

// Structure pour stocker les paramètres
struct Parameters {
    int n_ind;             // Nombre d'individus
    int n_loci;            // Nombre de loci
    int n_distinct_haplo;  // Nombre d'haplotypes distincts dans la population
    int n_twos;	   // Nombre de loci hétérozygotes
};

/**
 * Lit les paramètres nécessaires à partir d'un fichier texte.
 * Le fichier doit contenir les valeurs pour n_ind, n_loci, et n_distinct_haplo
 * sous la forme de lignes commentées avec des valeurs par défaut.
 *
 * @param filename Le chemin du fichier à lire.
 * @return Une structure Parameters contenant les valeurs lues.
 */

Parameters readParameters(const std::string &filename) {
    Parameters params;
    std::ifstream file(filename);
    
    if (!file) {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier " << filename << std::endl;
        exit(1);
    }

    std::string line;
    while (getline(file, line)) {
        if (line.find("n_ind") != std::string::npos) {
            params.n_ind = std::stoi(line.substr(line.find_last_of(" ") + 1));
        } else if (line.find("n_loci") != std::string::npos) {
            params.n_loci = std::stoi(line.substr(line.find_last_of(" ") + 1));
        } else if (line.find("n_distinct_haplo") != std::string::npos) {
            params.n_distinct_haplo = std::stoi(line.substr(line.find_last_of(" ") + 1));
        } else if (line.find("n_twos") != std::string::npos) {
            params.n_twos = std::stoi(line.substr(line.find_last_of(" ") + 1));
        }
        
    }


    file.close();
    return params;
}

#endif // PARAMETRE_H
