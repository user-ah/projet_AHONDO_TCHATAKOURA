#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>

void readParameters(const std::string& parametersFile, float& percentMiss) {
    std::ifstream file(parametersFile);
    if (!file) {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier " << parametersFile << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.find("percent_miss") != std::string::npos) {
            percentMiss = std::stof(line.substr(line.find('<') + 1));
            break;
        }
    }
    file.close();
}

void readGenotypes(const std::string& genotypesFile, std::vector<std::vector<int>>& genotypes) {
    std::ifstream file(genotypesFile);
    if (!file) {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier " << genotypesFile << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<int> row;
        for (char c : line) {
            if (std::isdigit(c)) {
                row.push_back(c - '0');
            }
        }
        genotypes.push_back(row);
    }
    file.close();
}

void introduceMissingData(std::vector<std::vector<int>>& genotypes, float percentMiss) {
    int totalElements = genotypes.size() * genotypes[0].size();
    int missingCount = static_cast<int>(percentMiss * totalElements / 100);

    std::srand(std::time(nullptr));
    while (missingCount > 0) {
        int i = std::rand() % genotypes.size();
        int j = std::rand() % genotypes[0].size();

        if (genotypes[i][j] != -1) { // -1 indique une donnée manquante
            genotypes[i][j] = -1;
            --missingCount;
        }
    }
}

void writeGenotypes(const std::string& outputFile, const std::vector<std::vector<int>>& genotypes) {
    std::ofstream file(outputFile);
    if (!file) {
        std::cerr << "Erreur : Impossible d'écrire dans le fichier " << outputFile << std::endl;
        exit(EXIT_FAILURE);
    }

    for (const auto& row : genotypes) {
        for (size_t i = 0; i < row.size(); ++i) {
            if (row[i] == -1) {
                file << "NA";
            } else {
                file << row[i];
            }
            if (i < row.size() - 1) {
                file << " ";
            }
        }
        file << "\n";
    }
    file.close();
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: ./alter_data <parameters.txt> <genotypes.csv>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string parametersFile = argv[1];
    std::string genotypesFile = argv[2];
    std::string outputFile = genotypesFile.substr(0, genotypesFile.find_last_of('.')) + "_altered.csv";

    float percentMiss = 0.0f;
    readParameters(parametersFile, percentMiss);

    std::vector<std::vector<int>> genotypes;
    readGenotypes(genotypesFile, genotypes);

    introduceMissingData(genotypes, percentMiss);

    writeGenotypes(outputFile, genotypes);

    std::cout << "Fichier de génotypes altéré écrit dans " << outputFile << std::endl;

    return 0;
}
