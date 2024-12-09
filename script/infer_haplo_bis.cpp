#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <cmath>
#include <functional>
#include <unordered_map>
#include <algorithm>
#include <climits>
#include <limits>

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
            genotype.push_back(std::stoi(val));
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

// Calculer la distance de Hamming
int calculateHammingDistance(const Haplotype& h1, const Haplotype& h2) {
    if (h1.size() != h2.size()) {
        throw std::invalid_argument("Les haplotypes doivent avoir la même taille pour calculer la distance de Hamming.");
    }

    int distance = 0;
    for (size_t i = 0; i < h1.size(); ++i) {
        if (h1[i] != h2[i]) {
            ++distance;
        }
    }
    return distance;
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

// Initialisation des fréquences des haplotypes
void initialisation_frequences_haplotypes(HaploFreqMap& frequencies, const HaploList& haplotypes) {
    double initial_freq = 1.0 / haplotypes.size();
    for (const auto& pair : haplotypes) {
        frequencies[pair.first] += initial_freq;
        frequencies[pair.second] += initial_freq;
    }
}

// Calculer la probabilité des génotypes
void calcul_proba_genotypes(
    const std::vector<std::vector<int>>& genotypes,
    const HaploFreqMap& frequencies,
    GenotypeProbMap& proba_genotypes
) {
    for (const auto& genotype : genotypes) {
        auto pairs = generateHaplotypes(genotype);
        double proba = 0.0;
        for (const auto& pair : pairs) {
            auto it1 = frequencies.find(pair.first);
            auto it2 = frequencies.find(pair.second);
            if (it1 != frequencies.end() && it2 != frequencies.end()) {
                double p1 = it1->second;
                double p2 = it2->second;
                proba += (pair.first == pair.second) ? p1 * p1 : 2 * p1 * p2;
            }
        }
        proba_genotypes[genotype] = proba;
    }
}

// Maximisation
void maximisation(
    const HaploFreqMap& freq_prec,
    const GenotypeProbMap& proba_prec,
    const GenotypeCountMap& genotype_count,
    int total_genotypes,
    HaploFreqMap& new_freq
) {
    for (auto& [haplo, freq] : new_freq) {
        freq = 0.0;
    }

    for (const auto& [genotype, count] : genotype_count) {
        auto pairs = generateHaplotypes(genotype);
        for (const auto& pair : pairs) {
            auto it1 = freq_prec.find(pair.first);
            auto it2 = freq_prec.find(pair.second);
            if (it1 != freq_prec.end() && it2 != freq_prec.end()) {
                double contrib = (pair.first == pair.second)
                                 ? it1->second * it1->second
                                 : 2 * it1->second * it2->second;
                auto it_proba = proba_prec.find(genotype);
                if (it_proba != proba_prec.end() && it_proba->second > 0) {
                    contrib /= it_proba->second;
                }
                contrib *= (count / static_cast<double>(total_genotypes));
                new_freq[pair.first] += contrib;
                new_freq[pair.second] += contrib;
            }
        }
    }

    for (auto& [haplo, freq] : new_freq) {
        freq /= 2.0;
    }
}

// Comparer les haplotypes inférés avec les haplotypes de base
void compareHaplotypes(
    const HaploFreqMap& inferred_haplotypes,
    const HaploList& base_haplotypes,
    const std::string& filename
) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier " << filename << " pour écrire." << std::endl;
        return;
    }

    for (const auto& [inferred_haplo, freq] : inferred_haplotypes) {
        int min_distance = INT_MAX;
        HaploPair closest_pair;

        for (const auto& base_pair : base_haplotypes) {
            int dist1 = calculateHammingDistance(inferred_haplo, base_pair.first);
            int dist2 = calculateHammingDistance(inferred_haplo, base_pair.second);

            int min_pair_distance = std::min(dist1, dist2);
            if (min_pair_distance < min_distance) {
                min_distance = min_pair_distance;
                closest_pair = base_pair;
            }
        }

        outFile << "Haplotype inféré : ";
        for (int val : inferred_haplo) {
            outFile << val << " ";
        }
        outFile << "(Fréquence : " << freq << ") -> Distance minimale : " << min_distance << std::endl;
    }

    outFile.close();
}

// Fonction principale
int main() {
    try {
        std::vector<std::vector<int>> genotypes = readGenotypesFromCSV("genotypes.csv");

        HaploList haplotypes;
        for (const auto& genotype : genotypes) {
            auto pairs = generateHaplotypes(genotype);
            haplotypes.insert(haplotypes.end(), pairs.begin(), pairs.end());
        }

        HaploList unique_haplotypes = removeDuplicates(haplotypes);

        HaploFreqMap frequencies;
        initialisation_frequences_haplotypes(frequencies, unique_haplotypes);

        GenotypeProbMap proba_genotypes;
        calcul_proba_genotypes(genotypes, frequencies, proba_genotypes);

        compareHaplotypes(frequencies, unique_haplotypes, "hamming_distances.txt");

        std::cout << "Calcul terminé." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Erreur: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
