import subprocess
import os

# Spécifiez les chemins des fichiers
parameters_file = os.path.expanduser("~/Documents/projet_AHONDO_TCHATAKOURA/script/parameters.txt")
genotypes_file = os.path.expanduser("~/Documents/projet_AHONDO_TCHATAKOURA/script/genotypes.csv")
haplotypes_file = os.path.expanduser("~/Documents/projet_AHONDO_TCHATAKOURA/script/haplotypes.csv")

# Chemin vers le fichier C++
cpp_file = os.path.expanduser("~/Documents/projet_AHONDO_TCHATAKOURA/script/create_geno_haplo_data.cpp")

# Commande pour compiler le programme C++
compile_command = f'g++ -o create_geno_haplo_data {cpp_file}'

# Exécution de la commande de compilation
subprocess.run(compile_command, shell=True, check=True)

# Exécuter le programme compilé
run_command = f'./create_geno_haplo_data {parameters_file} {genotypes_file} {haplotypes_file}'
subprocess.run(run_command, shell=True, check=True)
