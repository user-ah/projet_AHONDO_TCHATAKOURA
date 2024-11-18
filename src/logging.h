// logging.h
#ifndef LOGGING_H
#define LOGGING_H

#include <string>
#include <fstream>

/**
 * Écrit les informations de log dans un fichier.
 *
 * @param log_file Le chemin vers le fichier de log.
 * @param n_distinct_ground_truth Nombre d'haplotypes distincts dans la référence.
 * @param n_distinct_inferred Nombre d'haplotypes distincts dans la solution inférée.
 * @param distance La distance entre les haplotypes de référence et inférés.
 * @param duration Durée d'exécution de l'algorithme d'inférence.
 */
void writeLog(const std::string &log_file,
              int n_distinct_ground_truth,
              int n_distinct_inferred,
              int distance,
              double duration) {
    std::ofstream log(log_file);
    if (log.is_open()) {
        log << "n_distinct_haplo (référence): " << n_distinct_ground_truth << std::endl;
        log << "n_distinct_haplo (inféré): " << n_distinct_inferred << std::endl;
        log << "Distance entre haplotypes de référence et inférés: " << distance << std::endl;
        log << "Durée d'exécution: " << duration << " secondes" << std::endl;
        log.close();
        std::cout << "Log écrit dans " << log_file << std::endl;
    } else {
        std::cerr << "Échec de l'ouverture du fichier de log." << std::endl;
    }
}

#endif // LOGGING_H
