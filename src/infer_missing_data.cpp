#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <unordered_map>

using namespace std;

// Lecture du fichier parameters.txt pour obtenir la valeur de k
bool readParameters(const string &parameters_file, int &k) {
    ifstream file(parameters_file);
    if (!file.is_open()) {
        cerr << "Erreur : Impossible d'ouvrir " << parameters_file << endl;
        return false;
    }
    string line;
    while (getline(file, line)) {
        if (line.find("k") != string::npos) {
            istringstream iss(line);
            string key;
            iss >> key >> k;
        }
    }
    file.close();
    return true;
}

// Lecture du fichier CSV contenant les génotypes
// Lecture du fichier CSV contenant les génotypes avec des espaces comme séparateur
bool readGenotypes(const string &genotypes_file, vector<vector<int>> &genotypes) {
    ifstream file(genotypes_file);
    if (!file.is_open()) {
        cerr << "Erreur : Impossible d'ouvrir " << genotypes_file << endl;
        return false;
    }
    string line;
    while (getline(file, line)) {
        vector<int> row;
        stringstream ss(line);
        string cell;
        
        // Lire chaque cellule séparée par un espace
        while (ss >> cell) { // Utilise l'opérateur '>>' pour découper la ligne par espaces
            if (cell == "NA") {  // Valeur manquante
                row.push_back(-1);
            } else {
                row.push_back(stoi(cell)); // Convertir la cellule en entier
            }
        }
        genotypes.push_back(row);
    }
    file.close();
    return true;
}

// Écriture du fichier CSV des génotypes imputés
bool writeGenotypes(const string &output_file, const vector<vector<int>> &genotypes) {
    ofstream file(output_file);
    if (!file.is_open()) {
        cerr << "Erreur : Impossible d'ouvrir " << output_file << endl;
        return false;
    }
    for (const auto &row : genotypes) {
        for (size_t i = 0; i < row.size(); ++i) {
            if (row[i] == -1) {
                file << "NA";
            } else {
                file << row[i];
            }
            if (i < row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }
    file.close();
    return true;
}

// Calcul de la distance Hamming entre deux vecteurs (ignorer les valeurs manquantes)
int hammingDistance(const vector<int> &a, const vector<int> &b) {
    int distance = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != -1 && b[i] != -1 && a[i] != b[i]) {
            ++distance;
        }
    }
    return distance;
}

// KNN pour imputer une ligne manquante
void imputeMissingData(vector<int> &target, const vector<vector<int>> &genotypes, int k) {
    vector<pair<int, int>> distances; // {distance, index}

    for (size_t i = 0; i < genotypes.size(); ++i) {
        int distance = hammingDistance(target, genotypes[i]);
        distances.emplace_back(distance, i);
    }
    // Trier par distance croissante
    sort(distances.begin(), distances.end());

    // Parcourir chaque locus pour imputer
    for (size_t i = 0; i < target.size(); ++i) {
        if (target[i] == -1) { // Si la donnée est manquante
            vector<int> values; // Valeurs pour imputer
            for (int j = 0; j < k && j < distances.size(); ++j) {
                int neighbor_index = distances[j].second;
                if (genotypes[neighbor_index][i] != -1) { // Ne prendre que les données valides
                    values.push_back(genotypes[neighbor_index][i]);
                }
            }
            if (!values.empty()) {
                // Imputer avec la valeur majoritaire
                int ones = count(values.begin(), values.end(), 1);
                int zeros = count(values.begin(), values.end(), 0);
                target[i] = (ones >= zeros) ? 1 : 0;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <parameters.txt> <genotypes.csv> <output_genotypes.csv>" << endl;
        return 1;
    }

    string parameters_file = argv[1];
    string genotypes_file = argv[2];
    string output_file = argv[3];

    // Lire les paramètres
    int k = 0;
    if (!readParameters(parameters_file, k)) {
        return 1;
    }
    cout << "Valeur de k : " << k << endl;

    // Lire les génotypes
    vector<vector<int>> genotypes;
    if (!readGenotypes(genotypes_file, genotypes)) {
        return 1;
    }

    // Imputer les données manquantes pour chaque ligne
    for (auto &row : genotypes) {
        imputeMissingData(row, genotypes, k);
    }

    // Enregistrer le fichier avec les données imputées
    if (!writeGenotypes(output_file, genotypes)) {
        return 1;
    }

    cout << "Imputation terminée. Fichier enregistré : " << output_file << endl;
    return 0;
}
