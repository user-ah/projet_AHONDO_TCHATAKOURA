#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <string>
using namespace std;

// Fonction pour lire les paramètres du fichier
void read_parameters(const string& filepath, int& n_ind, int& n_loci, int& n_distinct_haplo) {
    ifstream infile(filepath);
    if (!infile) {
        cerr << "Erreur d'ouverture du fichier de paramètres.\n";
        exit(1);
    }
    string line;
    while (infile >> line) {
        if (line == "n_ind") {
            infile >> n_ind;
        } else if (line == "n_loci") {
            infile >> n_loci;
        } else if (line == "n_distinct_haplo") {
            infile >> n_distinct_haplo;
        }
    }
    infile.close();
}

// Fonction pour générer des haplotypes aléatoires
vector<vector<int>> generate_haplotypes(int n_ind, int n_loci, int n_distinct_haplo) {
    vector<vector<int>> distinct_haplotypes;
    
    // Générer des haplotypes distincts
    for (int i = 0; i < n_distinct_haplo; ++i) {
        vector<int> haplotype;
        for (int j = 0; j < n_loci; ++j) {
            haplotype.push_back(rand() % 2);  // 0 ou 1 aléatoire
        }
        distinct_haplotypes.push_back(haplotype);
    }

    // Assignation aléatoire des haplotypes
    vector<vector<int>> haplotypes;
    for (int i = 0; i < n_ind; ++i) {
        int haplo1_idx = rand() % n_distinct_haplo;
        int haplo2_idx = rand() % n_distinct_haplo;
        haplotypes.push_back(distinct_haplotypes[haplo1_idx]);
        haplotypes.push_back(distinct_haplotypes[haplo2_idx]);
    }

    return haplotypes;
}

// Fonction pour générer des génotypes à partir des haplotypes
vector<vector<int>> generate_genotypes(const vector<vector<int>>& haplotypes, int n_loci) {
    vector<vector<int>> genotypes;

    for (size_t i = 0; i < haplotypes.size(); i += 2) {
        vector<int> genotype;
        const vector<int>& haplo1 = haplotypes[i];
        const vector<int>& haplo2 = haplotypes[i + 1];

        for (int j = 0; j < n_loci; ++j) {
            if (haplo1[j] == 0 && haplo2[j] == 0) {
                genotype.push_back(0);
            } else if (haplo1[j] == 1 && haplo2[j] == 1) {
                genotype.push_back(1);
            } else {
                genotype.push_back(2);
            }
        }
        genotypes.push_back(genotype);
    }

    return genotypes;
}

// Fonction pour sauvegarder un fichier CSV
void save_csv(const string& filename, const  vector< vector<int>>& data) {
     ofstream file(filename);
    
    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
         cerr << "Usage: ./create_geno_haplo_data <chemin parameters.txt> <chemin genotypes.csv> <chemin haplotypes.csv>\n";
        return 1;
    }

    // Fichier de paramètres et chemins pour les fichiers de sortie
     string param_file = argv[1];
     string geno_file = argv[2];
     string haplo_file = argv[3];

    srand(time(0));  // Initialisation du générateur de nombres aléatoires

    // Variables pour les paramètres
    int n_ind, n_loci, n_distinct_haplo;

    // Lecture des paramètres
    read_parameters(param_file, n_ind, n_loci, n_distinct_haplo);

    // Génération des haplotypes et génotypes
     vector< vector<int>> haplotypes = generate_haplotypes(n_ind, n_loci, n_distinct_haplo);
     vector< vector<int>> genotypes = generate_genotypes(haplotypes, n_loci);

    // Sauvegarde des fichiers CSV
    save_csv(haplo_file, haplotypes);
    save_csv(geno_file, genotypes);

     cout << "Fichiers '" << geno_file << "' et '" << haplo_file << "' générés avec succès.\n";

    return 0;
}
