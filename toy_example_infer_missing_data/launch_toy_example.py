import os
import sys

# Variables pour les fichiers et paramètres
parameters_file = "../src/parameters.txt"
genotypes_file = "../data/genotypes_altered.csv"
output_genotypes_file = "../data/imputed_genotypes.csv"

# Nom de l'exécutable et code source
executable = "../src/infer_missing_data"
source_file = "../src/infer_missing_data.cpp"

# Compilation de l'exécutable
compile_cmd = f"g++ -o {executable} {source_file} -std=c++11"
print(f"Compilation de l'exécutable : {compile_cmd}")
compile_exit_code = os.system(compile_cmd)

# Vérification de la réussite de la compilation
if compile_exit_code != 0:
    sys.stderr.write("Erreur lors de la compilation du programme C++.\n")
    sys.exit(1)
else:
    print("Compilation réussie.")

# Création de la commande pour exécuter le programme
cmd = f"{executable} {parameters_file} {genotypes_file} {output_genotypes_file}"
print(f"Exécution de la commande : {cmd}")

# Exécution de la commande
exit_code = os.system(cmd)

# Vérification de la réussite de l'exécution
if exit_code != 0:
    sys.stderr.write(f"Erreur lors de l'exécution de {cmd}\n")
    sys.exit(1)
else:
    print(f"Exécution réussie. Fichier des génotypes imputés enregistré dans {output_genotypes_file}")
