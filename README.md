# master2_stage

Reprend les codes faits pendant le stage de M2.

- assignation_de_spectre : tout les scripts pour l'assignation de spectres de masse
- comparaison_taxonomie_phylogenie : le script permettant la comparaison d'une phylogénie et d'une taxonomie
- outils_peptide_marker_validation : les scripts côté html et cgi-bin de l'outil peptide_marker_validation à mettre sur le server bioinfodev 

### Prérequis

- Python 3 ou plus
- [Biopython](https://github.com/biopython/biopython)
- [RPG](https://rapid-peptide-generator.readthedocs.io/en/latest/)
- [Pyteomics](https://pyteomics.readthedocs.io/en/latest/index.html) (pour l'assignation de spectres de masse)
- [PTM_mass_simulator](https://gitlab.univ-lille.fr/bilille/ptm_mass_simulator) (pour le calcule des masses avec PTM : optionnel)


## Utilisation

### assignation_de_spectre

Dans main_pep_liste_search_in_spectra.py:

    python main_pep_liste_search_in_spectra.py -i PEPTIDESCSVFILE -o RESULTSPATH -t TAXONOMICTREE -s SEQUENCEFASTA -p PEAKMS [-r PTM] [-m THRESHOLD] 

- PEPTIDESCSVFILE : le chemin du fichier csv peptide généré par RPG
- RESULTSPATH : le dossier résultat
- TAXONOMICTREE : le chemin du fichier de l'arbre taxonomique au format Newick
- SEQUENCEFASTA : le chemin du fichier multifasta contenant les séquences protéiques
- PEAKMS : le chemin du fichier du spectre de masse
- PTM : le chemin du fichier généré par PTM_mass_simulator (optionnel)
- THRESHOLD : la précision demandée lors du calcule et de la comparaison des masses sous la forme ".{}f" (défaut : 3 chiffres après la virgule, sous la forme : ".3f")

Fichiers résultat:
- peptides.csv / peptides.fasta : les peptides avec leur séquence, leurs masses, les LCA des espèces dans lesquelles le peptide apparait 
- taxonomy.xml : le fichier de l'arbre taxonomique converti au format phyloXML avec les séquences protéiques associées aux feuilles de l'arbre taxonomie
- taxonomy_with_peak.xml : l'arbre taxonomique avec les pics du spectre de masse associés aux peptides dans les séquences 
- terminal_nodes.csv : les noeuds terminaux dans l'arbre taxonomique qui n'ont aucun pic associé à un peptide en dessous d'eux, leur niveau dans la taxonomie, et le nombre de peptides associés à un pic entre la racine et le clade

### outils_peptide_marker_validation

- dossier html : le formulaire + page html + css
- dossier cgi-bin : wrapper + scripts de l'outil

Possibilité de lancer les scripts sans passer par le formulaire :

    python main_peptide_marker_validation.py -type TYPEOFSEQ -fpep NAMEPEP -fprot NAMESEQ -output NAMEOUTPUT -model NAMEMODEL [-e ENZYME] [-m MASSTHRESHOLD][-PTM MASSPTM]

- TYPEOFSEQ : type de séquences utilisé ("proteins", "RNA" ou "DNA")
- NAMEPEP : chemin vers le fichier csv des peptides qui doit être formaté comme "peptide_input.csv"
- NAMESEQ : le chemin vers les séquences au format fasta
- NAMEOUTPUT : le nom du dossier résultat
- NAMEMODEL : le chemin vers le modèle utilisé pour avoir les protéines matures (dans le dossier files/) pour l'instant 2 modèles sont proposés : col1a1 et col1a2 de l'humain
- ENZYME : le numéro de l'enzyme utilisé par RPG 
- MASSTHRESHOLD : la précision demandé pour le calcule et la comparaison des masses ("1" ou "3" pour 1 chiffre ou 3 chiffres après la virgule)
- MASSPTM : pour la prise en compte de la masse des PTM (non pris en compte pour le moment dans le script)

Fichiers résultat:

- peptides.csv : la séquence du peptide, la masse, le nom des séquences dans lesquelles il apparait, la position de début et de fin de l'occurence du peptide dans la protéine mature + le nom de la séquence, le label (le nom du peptide donné par l'utilisateur), la séquence target + le nom de la séquence, la masse de la séquence target + l'espèce
- mature_prot.fasta : (seulement si utilisation de RPG) les séquences protéiques matures au format fasta
- rpg_converted.csv : (seulement si utilisation de RPG) le fichier résultat donné par RPG convertit dans le format CSV utilisé en input 
- rpg.csv : (seulement si utilisation de RPG) le fichier de sortie de RPG

### comparaison_taxonomie_phylogenie

    python comp_tree.py -p PHYLOGENETICTREE -t TAXONOMICTREE -s SEQUENCEFASTA -e PEPTIDESEQ
  
- PHYLOGENETICTREE : le chemin du fichier de l'arbre phylogénétique au format Newick
- TAXONOMICTREE : le chemin de l'arbre taxonomique au format Newick
- SEQUENCEFASTA : le chemin des séquences protéiques au format fasta
- PEPTIDESEQ : le chemin du fichier des peptides à vérifié dans la taxonomie (csv ou une ligne = combinaison de peptides uniques)

Exemple d'output : output_taxonomie_complete_specie.txt

