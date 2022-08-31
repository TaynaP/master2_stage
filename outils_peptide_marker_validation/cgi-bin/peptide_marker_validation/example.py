#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys

def main():        
    chemin_result = sys.argv[1]
    run_id = sys.argv[2]
    case = sys.argv[3]
    file_fasta = open(chemin_result+"/"+run_id+".fasta", "r")
    sequence = file_fasta.read()
    file_fasta.close()
    my_file = open(chemin_result+"results.php", "w")
    my_file.write("<div id='center'>")
    if case=="U":
        my_file.write("<p id='seq'>"+sequence.upper()+"</p>")
    else:
        my_file.write("<p id='seq'>"+sequence+"</p>")
    if len(sys.argv)==5:
        cg = sequence.count("g") + sequence.count("c")
        my_file.write("<br/><br/>")
        my_file.write("<p>The number of GC is %s</p>"%(cg))
    my_file.write("</div>")
    my_file.close()

main()
    
