import os

# Chemin vers le fichier parameters.txt
parameters_file = "~/Documents/projet_AHONDO_TCHATAKOURA/script/parameters.txt"

# Chemin vers les fichiers de sortie pour les génotypes et haplotypes
genotypes_file = "~/Documents/projet_AHONDO_TCHATAKOURA/script/genotypes.csv"
haplotypes_file = "~/Documents/projet_AHONDO_TCHATAKOURA/script/haplotypes.csv"

# Ligne de commande pour exécuter create_geno_haplo_data avec les arguments
command = f"./create_geno_haplo_data {parameters_file} {genotypes_file} {haplotypes_file}"

# Exécuter la commande dans le shell
os.system(command)
