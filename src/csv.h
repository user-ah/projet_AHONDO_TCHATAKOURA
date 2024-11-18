// csv.h
#ifndef CSV_H
#define CSV_H

#include <vector>
#include <string>
#include <fstream>

/**
 * Écrit les données dans un fichier CSV.
 *
 * @param filename Le nom du fichier de sortie.
 * @param data Un vecteur de vecteurs représentant les données à écrire.
 */
void writeCSV(const std::string &filename, const std::vector<std::vector<int>> &data) {
    std::ofstream file(filename);
    
    if (!file) {
        std::cerr << "Erreur : Impossible de créer le fichier " << filename << std::endl;
        return;
    }

    for (const auto &row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) file << ",";
        }
        file << "\n";
    }
    
    file.close();
}

#endif // CSV_H
