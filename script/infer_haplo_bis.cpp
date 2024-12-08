#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <cmath>
#include <functional>
#include <unordered_map>
#include <algorithm>

// Custom hash function for vector<int>
struct VectorHash {
    size_t operator()(const std::vector<int>& v) const {
        std::hash<int> hasher;
        size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

// Custom equality comparator for vector<int>
struct VectorEqual {
    bool operator()(const std::vector<int>& lhs, const std::vector<int>& rhs) const {
        return lhs == rhs;
    }
};

// Déclaration des structures
using Haplotype = std::vector<int>;
using HaploPair = std::pair<Haplotype, Haplotype>;
using HaploList = std::vector<HaploPair>;
using HaploFreqMap = std::unordered_map<Haplotype, double, VectorHash, VectorEqual>;
using GenotypeProbMap = std::unordered_map<std::vector<int>, double, VectorHash, VectorEqual>;
using GenotypeCountMap = std::unordered_map<std::vector<int>, int, VectorHash, VectorEqual>;

// Lire les génotypes depuis un fichier CSV
std::vector<std::vector<int>> readGenotypesFromCSV(const std::string& filename) {
    std::vector<std::vector<int>> genotypes;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::vector<int> genotype;
        std::stringstream ss(line);
        std::string val;
        while (std::getline(ss, val, ',')) {
            genotype.push_back(std::stoi(val));  // Convertir chaque valeur en int
        }
        genotypes.push_back(genotype);
    }
    return genotypes;
}

// Générer toutes les paires d'haplotypes possibles pour un génotype donné
HaploList generateHaplotypes(const std::vector<int>& genotype) {
    HaploList result;

    std::function<void(size_t, Haplotype, Haplotype)> backtrack =
        [&](size_t index, Haplotype h1, Haplotype h2) {
            if (index == genotype.size()) {
                result.emplace_back(h1, h2);
                return;
            }

            if (genotype[index] == 0) {
                h1.push_back(0);
                h2.push_back(0);
                backtrack(index + 1, h1, h2);
            } else if (genotype[index] == 1) {
                h1.push_back(1);
                h2.push_back(1);
                backtrack(index + 1, h1, h2);
            } else if (genotype[index] == 2) {
                Haplotype h1_copy = h1, h2_copy = h2;

                h1.push_back(0);
                h2.push_back(1);
                backtrack(index + 1, h1, h2);

                h1_copy.push_back(1);
                h2_copy.push_back(0);
                backtrack(index + 1, h1_copy, h2_copy);
            }
        };

    backtrack(0, {}, {});
    return result;
}

// Supprimer les doublons d'une liste de paires d'haplotypes
HaploList removeDuplicates(const HaploList& haplotypes) {
    std::set<HaploPair> unique_haplotypes;

    for (const auto& pair : haplotypes) {
        auto sorted_pair = (pair.first < pair.second) ? pair : std::make_pair(pair.second, pair.first);
        unique_haplotypes.insert(sorted_pair);
    }

    return HaploList(unique_haplotypes.begin(), unique_haplotypes.end());
}

// Initialiser les fréquences uniformes
HaploFreqMap initializeFrequencies(const HaploList& haplotypes) {
    HaploFreqMap frequencies;
    double uniform_frequency = 1.0 / (2 * haplotypes.size());

    for (const auto& pair : haplotypes) {
        const Haplotype& h1 = pair.first;
        const Haplotype& h2 = pair.second;

        if (frequencies.count(h1) == 0) {
            frequencies[h1] = 0.0;
        }
        if (frequencies.count(h2) == 0) {
            frequencies[h2] = 0.0;
        }
        frequencies[h1] += uniform_frequency;
        frequencies[h2] += uniform_frequency;
    }

    return frequencies;
}

// Trier les haplotypes par fréquence décroissante
std::vector<std::pair<Haplotype, double>> sortHaplotypeFrequencies(const HaploFreqMap& haplo_freq) {
    std::vector<std::pair<Haplotype, double>> sorted_haplotypes(haplo_freq.begin(), haplo_freq.end());

    // Trier par fréquence décroissante
    std::sort(sorted_haplotypes.begin(), sorted_haplotypes.end(),
              [](const std::pair<Haplotype, double>& a, const std::pair<Haplotype, double>& b) {
                  return a.second > b.second;
              });

    return sorted_haplotypes;
}

// Écrire les haplotypes triés dans un fichier
void writeHaplotypeFrequenciesToFile(const std::vector<std::pair<Haplotype, double>>& sorted_haplotypes,
                                     const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier " << filename << " pour écrire." << std::endl;
        return;
    }

    for (const auto& [haplotype, freq] : sorted_haplotypes) {
        for (int val : haplotype) {
            outFile << val << " ";
        }
        outFile << ": " << freq << std::endl;
    }

    outFile.close();
    std::cout << "Les haplotypes triés ont été sauvegardés dans le fichier : " << filename << std::endl;
}

// Maximisation de l'algorithme EM
void maximisation(
    const HaploFreqMap& freq_prec,
    const GenotypeProbMap& proba_prec,
    const GenotypeCountMap& genotype_count,
    int N,
    HaploFreqMap& new_freq,
    double seuil,
    int max_iter
) {
    double diff = seuil + 1.0; // Initialiser diff à une valeur plus grande que le seuil
    int iteration = 0;

    // Initialiser new_freq avec des zéros
    for (const auto& [h, freq] : freq_prec) {
        new_freq[h] = 0.0;
    }

    while (diff > seuil && iteration < max_iter) {
        iteration++;
        diff = 0.0;

        // Réinitialiser new_freq à zéro
        for (auto& [h, freq] : new_freq) {
            freq = 0.0;
        }

        for (const auto& [h1, freq_h1] : freq_prec) {
            for (const auto& [genotype, count] : genotype_count) {
                auto genotype_haplotypes = generateHaplotypes(genotype);

                for (const auto& [h_a, h_b] : genotype_haplotypes) {
                    if ((h_a == h1 || h_b == h1) && proba_prec.count(genotype) > 0) {
                        double p1 = freq_prec.at(h1);
                        double p2 = freq_prec.at(h_a == h1 ? h_b : h_a);
                        double prob_genotype = proba_prec.at(genotype);

                        double contrib = (h_a == h_b ? std::pow(p1, 2) : 2 * p1 * p2) 
                                         / prob_genotype * (count / static_cast<double>(N));

                        new_freq[h1] += contrib;
                    }
                }
            }
        }

        // Calculer la différence maximale
        for (const auto& [h, freq] : new_freq) {
            if (std::fabs(freq - freq_prec.at(h)) > diff) {
                diff = std::fabs(freq - freq_prec.at(h));
            }
        }
    }
}

int main() {
    try {
        // Lire les génotypes depuis un fichier
        std::vector<std::vector<int>> genotypes = readGenotypesFromCSV("genotypes.csv");

        // Générer les paires d'haplotypes possibles pour chaque génotype
        HaploList haplotypes;
        for (const auto& genotype : genotypes) {
            HaploList genotype_haplotypes = generateHaplotypes(genotype);
            haplotypes.insert(haplotypes.end(), genotype_haplotypes.begin(), genotype_haplotypes.end());
        }

        // Retirer les doublons des paires d'haplotypes
        HaploList unique_haplotypes = removeDuplicates(haplotypes);

        // Initialiser les fréquences des haplotypes
        HaploFreqMap frequencies = initializeFrequencies(unique_haplotypes);

        // Maximisation (EM)
        HaploFreqMap new_freq;
        maximisation(frequencies, {}, {}, genotypes.size(), new_freq, 1e-3, 100);

        // Trier les haplotypes par fréquence décroissante
        auto sorted_haplotypes = sortHaplotypeFrequencies(new_freq);

        // Sauvegarder les haplotypes triés dans un fichier
        writeHaplotypeFrequenciesToFile(sorted_haplotypes, "haplotype_frequencies.txt");

        // Afficher les haplotypes triés
        std::cout << "Haplotypes triés par fréquence décroissante : " << std::endl;
        for (const auto& [haplotype, freq] : sorted_haplotypes) {
            for (int val : haplotype) {
                std::cout << val << " ";
            }
            std::cout << ": " << freq << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Erreur: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
