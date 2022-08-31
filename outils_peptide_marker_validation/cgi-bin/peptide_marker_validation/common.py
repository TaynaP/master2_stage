#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

#Remplacer exemple_web_server par le nom de votre logiciel
HTML_PATH = "/bio2/www/html/peptide_marker_validation"
CGI_PATH = "/bio2/www/cgi-bin/peptide_marker_validation"

TEMP_DIR = os.path.join(CGI_PATH, 'tmp/')
RESULT_DIR = os.path.join(HTML_PATH, "result/")



#Renvoie le pwd pour le dossier tmp/
def tmp_dir(run_id):    
    tmpdir = os.path.join(TEMP_DIR, run_id)
    return tmpdir


#Renvoie le pwd du dossier result/
def result_dir(run_id):    
    resdir = os.path.join(RESULT_DIR, run_id)
    return resdir


#Retourne le pwd du dossier où sont stockés les résultats de la requête
def result_dir_html(run_id):    
    resdir = os.path.join(RESULT_DIR, run_id)
    return resdir

#Retourne le pwd du fichier fasta contenant les séquences soumises
def fasta_path_result(run_id, type):
    result_dir = result_dir_html(run_id)
    fpath = os.path.join(result_dir, run_id + type + ".fasta")
    return fpath

#Retourne le pwd du fichier csv contenant les séquences des peptides soumis
def csv_pep_path_result(run_id):
    result_dir = result_dir_html(run_id)
    csv_pep_path = os.path.join(result_dir, run_id + ".csv")
    return csv_pep_path
