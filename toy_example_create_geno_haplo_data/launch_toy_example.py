import subprocess
import os
import sys

# Utilisez des arguments de ligne de commande pour les chemins des fichiers
if len(sys.argv) != 5:
    print("Usage: python script.py <parameters_file> <genotypes_file> <haplotypes_file> <cpp_file>")
    sys.exit(1)

parameters_file = os.path.expanduser(sys.argv[1])
genotypes_file = os.path.expanduser(sys.argv[2])
haplotypes_file = os.path.expanduser(sys.argv[3])
cpp_file = os.path.expanduser(sys.argv[4])

# Commande pour compiler le programme C++
compile_command = f'g++ -o create_geno_haplo_data {cpp_file}'

# Exécution de la commande de compilation
subprocess.run(compile_command, shell=True, check=True)

# Exécuter le programme compilé
run_command = f'./create_geno_haplo_data {parameters_file} {genotypes_file} {haplotypes_file}'
subprocess.run(run_command, shell=True, check=True)
