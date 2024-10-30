// number of individuals
n_ind 10
// number of loci
n_loci 8
// number of haplotypes in the population of individuals
n_distinct_haplo 5



#include <iostream>
#include <vector>
#include <fstream>
#include <random>

using namespace std;

// Fonction pour générer des haplotypes aléatoires
vector<vector<int>> generate_haplotypes(int n_ind, int n_loci, int n_distinct_haplo) {
    vector<vector<int>> haplotypes;

    // Générer n_distinct_haplo haplotypes distincts aléatoirement
    for (int i = 0; i < n_distinct_haplo; i++) {
        vector<int> haplotype;
        for (int j = 0; j < n_loci; j++) {
            int allele = rand() % 2; // Valeur aléatoire 0 ou 1
            haplotype.push_back(allele);
        }
        haplotypes.push_back(haplotype); // Ajouter l'haplotype au vecteur
    }

    // Associer deux haplotypes pour chaque individu
    vector<vector<int>> paired_haplotypes;
    for (int i = 0; i < n_ind; i++) {
        int h1 = rand() % n_distinct_haplo;
        int h2 = rand() % n_distinct_haplo;
        paired_haplotypes.push_back(haplotypes[h1]); // Haplotype 1
        paired_haplotypes.push_back(haplotypes[h2]); // Haplotype 2
    }

    return paired_haplotypes;
}

// Fonction pour générer des génotypes à partir de deux haplotypes
vector<vector<int>> generate_genotypes(const vector<vector<int>>& haplotypes, int n_ind, int n_loci) {
    vector<vector<int>> genotypes;

    // Pour chaque individu (chaque paire d'haplotypes)
    for (int i = 0; i < n_ind; i++) {
        vector<int> genotype;
        const vector<int>& haplotype1 = haplotypes[2 * i];     // Haplotype 1
        const vector<int>& haplotype2 = haplotypes[2 * i + 1]; // Haplotype 2

        // Générer le génotype à partir des deux haplotypes
        for (int j = 0; j < n_loci; j++) {
            if (haplotype1[j] == 0 && haplotype2[j] == 0) {
                genotype.push_back(0); // 0 et 0 -> génotype 0
            } else if (haplotype1[j] == 1 && haplotype2[j] == 1) {
                genotype.push_back(1); // 1 et 1 -> génotype 1
            } else {
                genotype.push_back(2); // 0 et 1 ou 1 et 0 -> génotype 2
            }
        }

        genotypes.push_back(genotype); // Ajouter le génotype à la liste
    }

    return genotypes;
}
