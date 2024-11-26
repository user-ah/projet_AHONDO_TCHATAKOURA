import os
import sys
import subprocess

# Variables pour les fichiers et paramètres
parameters_file = "../script/parameters.txt"  # Fichier de paramètres
genotypes_file = "../data/genotypes.csv"  # Fichier de sortie pour les génotypes
haplotypes_file = "../data/haplotypes.csv"  # Fichier de sortie pour les haplotypes
cpp_file = "../script/create_geno_haplo_data.cpp"  # Fichier source C++

# Nom de l'exécutable
executable = "../script/create_geno_haplo_data"  # Nom de l'exécutable après compilation

# Commande de compilation pour le fichier C++
compile_command = f"g++ -o {executable} {cpp_file}"

# Vérification de la compilation
try:
    print(f"Compilation de {cpp_file}...")
    subprocess.run(compile_command, shell=True, check=True)
    print("Compilation réussie.")
except subprocess.CalledProcessError as e:
    sys.stderr.write(f"Erreur lors de la compilation du fichier C++ : {cpp_file}\n{e}\n")
    sys.exit(1)

# Création de la commande pour exécuter l'exécutable compilé
cmd = f"{executable} {parameters_file} {genotypes_file} {haplotypes_file}"

# Affichage de la commande pour validation
print(f"Exécution de la commande : {cmd}")

# Exécution de la commande
try:
    subprocess.run(cmd, shell=True, check=True)
    print("Exécution réussie.")
    print(f"Fichiers générés :\n - {genotypes_file}\n - {haplotypes_file}")
except subprocess.CalledProcessError as e:
    sys.stderr.write(f"Erreur lors de l'exécution de la commande : {cmd}\n{e}\n")
    sys.exit(1)