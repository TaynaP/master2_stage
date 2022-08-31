#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import time
import random
import json
import common


#Permet de créer le dossier de résultats dans /html/nom_logiciel/result/
def new_run():    
    runid = new_run_id()
    resdir = common.result_dir(runid)
    tmpdir = common.tmp_dir(runid)
    os.mkdir(resdir)
    os.mkdir(tmpdir)
    return runid


#Création d'un nombre au hasard avec la date et l'heure de la requête
def new_run_id():
    """ generate a job id """
    return time.strftime('%b_%d_%Y_%H_%M_%S') + "_%i" % random.randint(1,32768)


def main():
    run_id=""
    run_id = new_run()
    result={}
    if run_id!="":
        result['success'] = True
        result['run_id'] = run_id
    else:
        result['sucess'] = False
        
    sys.stdout.write("Content-Type: application/json")
    sys.stdout.write("\n\n")
    sys.stdout.write(json.dumps(result, indent=1))
    sys.stdout.write("\n")
    
main()
    
