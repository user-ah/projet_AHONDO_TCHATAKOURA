import os
import sys
import subprocess

# Variables pour les fichiers et paramètres
parameters_file = "../src/parameters.txt"
genotypes_file = "../data/genotypes.csv"
haplotypes_file = "../data/haplotypes.csv"
log_file = "../src/log.txt"
cpp_file = "../src/infer_haplo.cpp"  # Fichier source C++

# Nom de l'exécutable
executable = "../src/infer_haplo"

# Commande de compilation pour le fichier C++
compile_command = f"g++ -o {executable} {cpp_file}"

# Création de la commande
cmd = f"{executable} {parameters_file} {genotypes_file} {haplotypes_file} {log_file}"

# Vérification de la compilation
try:
    print(f"Compilation de {cpp_file}...")
    subprocess.run(compile_command, shell=True, check=True)
    print("Compilation réussie.")
except subprocess.CalledProcessError as e:
    sys.stderr.write(f"Erreur lors de la compilation du fichier C++ : {cpp_file}\n{e}\n")
    sys.exit(1)

# Affichage de la commande pour validation
print(f"Exécution de la commande : {cmd}")

# Exécution de la commande
try:
    subprocess.run(cmd, shell=True, check=True)
    print("Exécution réussie.")
    print(f"Fichiers générés :\n - {log_file}\n - {parameters_file}")
except subprocess.CalledProcessError as e:
    sys.stderr.write(f"Erreur lors de l'exécution de la commande : {cmd}\n{e}\n")
    sys.exit(1)