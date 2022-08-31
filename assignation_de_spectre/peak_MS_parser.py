from pep_parser import *
from peak import Peak
from pyteomics import mzml, auxiliary
import csv
import os


def peak_parser_csv(peak_file, threshold):
    """Parser for the spectra in csv 

    Args:
        peak_file (str): path to the peak file
        threshold (str): the precision chosen ".3f" or ".1f"

    Raises:
        NameError: if the file is not a csv

    Returns:
        list: list of peaks obj corresponding to the peaks in the spectra
    """
    _, ext = os.path.splitext(peak_file)
    if ext != ".csv":
        raise NameError("the input file must be a csv file")
    list_peak = []
    f= open (peak_file)
    myReader = csv.reader(f)
    next(myReader)
    for row in myReader:
        mass = format(float(row[0]), threshold)
        intensity = row[1]
        list_peak.append(Peak(mass, intensity))

    return list_peak

def peak_parser_mgf(peak_file, threshold):
    """Parser for the spectra in mgf

    Args:
        peak_file (str): path to the peak file
        threshold (str): the precision chosen ".3f" or ".1f"

    Returns:
        list: list of peaks obj corresponding to the peaks in the spectra
    """
    list_peak = []
    with open(peak_file, "r") as filin: #we read the file
        l = filin.readline() #we begin to read the 1st line
        while l != "":
            ligne = l.rstrip()
            if ligne != "": 
                if ligne[0].isdigit():
                    mass = format(float(ligne.split()[0]), threshold)
                    intensity = ligne.split()[1]
                    list_peak.append(Peak(mass, intensity))
            l = filin.readline() #we begin to read the 2nd line
            
    return list_peak


def peak_parser_mzml(peak_file, threshold):
    """Parser for the spectra in mzML

    Args:
        peak_file (str): path to the peak file
        threshold (str): the precision chosen ".3f" or ".1f"

    Returns:
        list: list of peaks obj corresponding to the peaks in the spectra
    """
    liste_peak = []
    f = mzml.MzML(peak_file, use_index=True) #utilisation du module mzml de pyteomics
    mass_array = f[0]["m/z array"] #pour avoir les m/z des pics
    intensity_array = f[0]["intensity array"] #pour avoir l'intensité des pics
 
    for mass, intensity in zip(mass_array.tolist(), intensity_array.tolist()):
        peak = Peak(format(float(mass), threshold), intensity) #on fait en sorte d'avoir la masse dans la même précision que demandé
        if peak not in liste_peak:
            liste_peak.append(peak)
          
    return liste_peak


def main_peak_parser(peak_file, threshold):
    """AI is creating summary for main_peak_parser

    Args:
        peak_file (str): the path to the spectra file
        threshold (str): the precision chosen ".3f" or ".1f"

    Returns:
        list: list of peak obj corresponding to the peak in the spectra
    """

    _, ext = os.path.splitext(peak_file)
    if ext == ".csv": # si le fichier de spectre est un csv
        list_peak = peak_parser_csv(peak_file, threshold)
    elif ext == ".mgf": # si le fichier de spectre est un mgf
        list_peak = peak_parser_mgf(peak_file, threshold)
    elif ext == ".mzML": # si le fichier de spectre est un mzML
        list_peak = peak_parser_mzml(peak_file, threshold)

    return list_peak
    
