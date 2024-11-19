import os
import sys
import subprocess

# Variables pour les fichiers et paramètres
parameters_file = "../src/parameters.txt"
genotypes_file = "../data/genotypes.csv"
cpp_file = "../src/alter_data.cpp"

# Nom de l'exécutable
executable = "../src/alter_data"  # L'exécutable après compilation

# Commande de compilation pour le fichier C++
compile_command = f"g++ -std=c++11 -o {executable} {cpp_file}"

# Vérification de la compilation
try:
    print(f"Compilation de {cpp_file}...")
    subprocess.run(compile_command, shell=True, check=True)
    print("Compilation réussie.")
except subprocess.CalledProcessError:
    sys.stderr.write(f"Erreur lors de la compilation du fichier C++ : {cpp_file}\n")
    sys.exit(1)

# Création de la commande pour exécuter l'exécutable compilé
cmd = f"{executable} {parameters_file} {genotypes_file}"

# Affichage de la commande pour validation
print(f"Exécution de la commande : {cmd}")

# Exécution de la commande
try:
    subprocess.run(cmd, shell=True, check=True)
    print("Exécution réussie.")
except subprocess.CalledProcessError:
    sys.stderr.write(f"Erreur lors de l'exécution de la commande : {cmd}\n")
    sys.exit(1)
