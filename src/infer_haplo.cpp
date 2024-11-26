#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <chrono>
#include "genotype_reader.h"         
#include "initialisation_frequences.h"
#include "calcul_proba_genotypes.h"
#include "maximisation.h"
#include "estimation_esperance.h"
#include "haplotype_generator.h"
#include "csv.h"
#include "haplotype_counter.h"
#include "distance_calculation.h"
#include "logging.h"

using namespace std;
using namespace std::chrono;

// Fonction pour générer le fichier parameters.txt
void generateParameters(const string& genotypes_file, const string& parameters_file, int& n_ind, int& n_loci) {
    vector<vector<int>> genotypes;

    // Lecture des données de génotypes
    if (!readGenotypeFile(genotypes_file, genotypes, n_ind, n_loci)) {
        throw runtime_error("Impossible de lire le fichier de génotypes.");
    }

    // Calcul du nombre d'individus et de loci
    n_ind = genotypes.size();
    n_loci = genotypes.empty() ? 0 : genotypes[0].size();

    // Calcul du nombre d'haplotypes distincts
    int n_distinct_haplo = countDistinctHaplotypes(genotypes);

    // Écriture des paramètres dans parameters.txt
    ofstream params_out(parameters_file);
    if (!params_out) {
        throw runtime_error("Impossible d'écrire le fichier parameters.txt.");
    }

    params_out << "// number of individuals\nn_ind " << n_ind << "\n";
    params_out << "// number of loci\nn_loci " << n_loci << "\n";
    
    params_out.close();
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        cerr << "Usage: " << argv[0] << " <parameters.txt> <genotypes.csv> <haplotypes.csv> <log.txt>" << endl;
        return 1;
    }

    // Récupération des arguments de fichier
    string parameters_file = argv[1];
    string genotypes_file = argv[2];
    string haplotypes_file = argv[3];
    string log_file = argv[4];

    // Génération de parameters.txt à partir de genotypes.csv
    int n_ind, n_loci;
    try {
        generateParameters(genotypes_file, parameters_file, n_ind, n_loci);
    } catch (const exception& e) {
        cerr << e.what() << endl;
        return 1;
    }

    cout << "Fichier " << parameters_file << " généré avec succès." << endl;
    cout << "Nombre d'individus (n_ind) : " << n_ind << endl;
    cout << "Nombre de loci (n_loci) : " << n_loci << endl;

    // Lecture des génotypes
    vector<vector<int>> genotypes;
    if (!readGenotypeFile(genotypes_file, genotypes, n_ind, n_loci)) {
        cerr << "Échec de la lecture du fichier de génotypes." << endl;
        return 1;
    }

    // Initialisation des fréquences d'haplotypes
    int nb_haplotypes = 1 << n_loci;
    vector<double> haplotype_frequencies;
    initialisation_frequences_haplotypes(haplotype_frequencies, nb_haplotypes);

    // Initialisation des probabilités des génotypes
    vector<double> genotype_probabilities;
    calcul_proba_genotypes(genotypes, haplotype_frequencies, genotype_probabilities);

    // Démarrage du chronomètre pour la durée d'exécution
    auto start_time = high_resolution_clock::now();

    // Boucle de l'algorithme EM
    bool convergence = false;
    int nb_etapes = 1;
    double loglikelihood_prec = -std::numeric_limits<double>::infinity();
    double loglikelihood;

    while (!convergence && nb_etapes <= 100) { // Limite maximale de 100 étapes pour l'EM
        cout << "Étape " << nb_etapes << " :" << endl;

        // Étape de maximisation
        maximisation(genotypes, haplotype_frequencies, genotype_probabilities, n_ind);

        // Étape d'estimation de l'espérance et calcul de la log-vraisemblance
        loglikelihood = estimation_esperance(genotypes, haplotype_frequencies);
        cout << "Log-vraisemblance : " << loglikelihood << endl;

        // Vérification de la convergence
        double difference = std::abs(loglikelihood_prec - loglikelihood);
        convergence = (difference / std::abs(loglikelihood_prec)) < 0.001; // Seuil de convergence

        if (!convergence) {
            loglikelihood_prec = loglikelihood;
            calcul_proba_genotypes(genotypes, haplotype_frequencies, genotype_probabilities);
            ++nb_etapes;
        }
    }

    // Fin du chronomètre
    auto end_time = high_resolution_clock::now();
    double duration = duration_cast<milliseconds>(end_time - start_time).count() / 1000.0;

    // Génération des haplotypes inférés en utilisant les fréquences finales
    vector<vector<int>> inferred_haplotypes(2 * n_ind, vector<int>(n_loci));
    for (int i = 0; i < n_ind; ++i) {
        vector<int> haplotype1, haplotype2;
        generateHaplotypes(genotypes[i], haplotype1, haplotype2);
        inferred_haplotypes[2 * i] = haplotype1;
        inferred_haplotypes[2 * i + 1] = haplotype2;
    }

    // Enregistrement des haplotypes inférés dans inferred_haplotypes.csv
    string inferred_haplotypes_file = "../data/inferred_haplotypes.csv";
    writeCSV(inferred_haplotypes_file, inferred_haplotypes);
    cout << "Haplotypes inférés enregistrés dans " << inferred_haplotypes_file << endl;

    // Calcul du nombre d'haplotypes distincts dans la référence et les haplotypes inférés
    vector<vector<int>> ground_truth_haplotypes;
    readGenotypeFile(haplotypes_file, ground_truth_haplotypes, n_ind, n_loci);

    int n_distinct_ground_truth = countDistinctHaplotypes(ground_truth_haplotypes);
    int n_distinct_inferred = countDistinctHaplotypes(inferred_haplotypes);

    // Calcul de la distance entre les haplotypes de référence et les haplotypes inférés
    int distance = calculateHaplotypeDistance(ground_truth_haplotypes, inferred_haplotypes);

    // Écriture des informations dans le fichier log.txt
    writeLog(log_file, n_distinct_ground_truth, n_distinct_inferred, distance, duration);

    return 0;
}
