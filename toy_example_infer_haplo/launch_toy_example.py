import os
import sys

# Variables pour les fichiers et paramètres
parameters_file = "../src/parameters.txt"
genotypes_file = "../data/genotypes.csv"
haplotypes_file = "../data/haplotypes.csv"
log_file = "../src/log.txt"

# Nom de l'exécutable
executable = "../src/infer_haplo"

# Création de la commande
cmd = f"{executable} {parameters_file} {genotypes_file} {haplotypes_file} {log_file}"

# Affichage de la commande pour validation
print(f"Exécution de la commande : {cmd}")

# Exécution de la commande
exit_code = os.system(cmd)

# Vérification de la réussite de l'exécution
if exit_code != 0:
    sys.stderr.write(f"Erreur lors de l'exécution de {cmd}\n")
    sys.exit(1)
else:
    print(f"Exécution réussie. Résultats et log enregistrés dans {log_file}")
