import numpy as np
import matplotlib.pyplot as plt
import re
from Bio import SeqIO
from Bio.Seq import Seq
from openpyxl import load_workbook
from Bio import motifs
MC = []

print("G4PromFinder searches for promoters in intergenic regions of a bacterial genome. It is recommended for GC rich genomes\n\n")
print("These G4PromFinder version consider all the elements identified in a determined IR")
print("Input: Genome file and annotation file")
print("Output: a text file containing promoter coordinates and a file containing informations about them")
print("Genome file must be in fasta format")
print("Annotation file must be in txt format. It must be organized in 3 columns. The first column must contain CDS start, the second CDS end and the third the strand (+1 for strand+ and -1 for strand-)")
File_genoma = input("Enter the input genome file name: ")

try:
    for seq_record in SeqIO.parse(File_genoma, "fasta"):
        MC.append(str(seq_record.seq))
    genoma = MC[0] 
    gc = 100*(genoma.count("G")+genoma.count("C"))/len(genoma)
    File_geni = input("Enter the input annotation file name  : ")
    try:
        n = np.loadtxt(File_geni)
        m = np.shape(n)
        m1 = m[0]
        pattern1 = "G{2,4}\D{1,10}G{2,4}\D{1,10}G{2,4}\D{1,10}G{2,4}"
        pattern2 = "C{2,4}\D{1,10}C{2,4}\D{1,10}C{2,4}\D{1,10}C{2,4}"
        patternA = "TA\D{3}T"
        patternB = "A\D{3}TA"
        patternC = "TTGAC"
        patternD = "GTCAA"
        lista1 = []
        lista2 = []
        distanza1 = []
        distanza2 = []
        distanze = []
        intergeniche = 0
        
        for i in range (m1):
            u = int(n[i,0])
            v = int(n[i,1])
            if n[i,2] == -1:
                n[i,0] = v
                n[i,1] = u
                lista2.append(i)
            else:
                lista1.append(i)
                
        
        
        promotori = 0 
        regione = []
        number_predictions = []
        SI1 = 0
        SI2 = 0
        NO1 = 0
        NO2 = 0
        cod1 = 1
        cod2 = 1
        ATregions = []
        print ("G4PromFinder is working, please wait...")
        
        q = open("Promoter coordinates (Method 1).txt","a")
        q.write("Region         Gene           Strand              Start               End\n")
        
        for i in lista1:
            u = int(n[i,0])
            v = int(n[i,1])
            if i == 0:
                z1 = genoma[:u]  # z1 = intergenic region
                inizio = 0
            else:
                u1 = int(n[i-1,0])
                v1 = int(n[i-1,1])
                if u1 < v1:
                    z1 = genoma[v1:u]
                    inizio = v1
                if u1 > v1:
                    z1 = genoma[u1:u]
                    inizio = u1
            if len(z1) >= 50:
                intergeniche += 1
            scale = 0
            count = 0
            for j in range(100000):
                x1 = scale + j
                x2 = scale + j + 25
                if x2 > len(z1) or len(z1) < 50:
                    break
                else:
                    w = z1[x1:x2]
                    at =100*(w.count("A")+w.count("T"))/len(w)
                    if x1 >= 50:
                        prom = z1[x1-50:x2]
                    else:
                        prom = z1[:x2]
                    if at >= 40:
                        if re.search(pattern1,prom) or re.search(pattern2,prom):
                            count += 1
                            if regione.count(i+1) == 0:
                                regione.append(i+1)
                            scale = x2 - j - 1
                            dati = []
                            for g in range(25):
                                w = z1[x1+g:x2+g]
                                at =100*(w.count("A")+w.count("T"))/len(w)
                                if x1+g >= 50:
                                    prom = z1[x1+g-50:x2+g]
                                else:
                                    prom = z1[:x2+g]
                                if at >= 40:
                                    if re.search(pattern2,prom) or re.search(pattern1,prom):
                                        dati.append(at)
                                    else:
                                        dati.append(0)
                            maxP = np.argmax(dati)
                            at = np.max(dati)
                            x1 = x1 + maxP
                            x2 = x2 + maxP
                            if x1 >= 50:
                                prom = z1[x1-50:x2]
                                x = x1-50+inizio
                                y = x2+inizio
                            else:
                                prom = z1[:x2]
                                x = inizio
                                y = x2+inizio
                            q.write("\nP")
                            q.write(str(cod1))
                            q.write("plus             ")
                            q.write(str(i+1))
                            q.write("             positive           ")
                            q.write(str(x))
                            q.write("           ")
                            q.write(str(y))
                            cod1 += 1
                            ATregions.append(at)
                            promotori += 1
                            if re.search(patternA,prom):
                                SI1 += 1
                            else:
                                NO1 += 1
                            if re.search(patternC,prom):
                                SI2 += 1
                            else:
                                NO2 += 1
            if len(z1) >= 50:
                number_predictions.append(count)        
            
        for i in lista2:
            u = int(n[i,0])
            v = int(n[i,1])
            if i < m1-1:
                u2 = int(n[i+1,0])
                v2 = int(n[i+1,1])
                if u2 < v2:
                    z1 = genoma[u:u2]
                else:
                    z1 = genoma[u:v2]
            else:
                z1 = genoma[u:]
            if len(z1) >= 50:
                intergeniche += 1
            scale = 0
            count = 0
            for j in range(10000):
                x1 = scale + j
                x2 = scale + j + 25
                if x2 > len(z1) or len(z1) < 50:
                    break
                else:
                    w = z1[x1:x2]
                    at =100*(w.count("A")+w.count("T"))/len(w)
                    if x2 + 50 <= len(z1):
                        prom = z1[x1:x2+50]
                    else:
                        prom = z1[x1:]
                    if at >= 40:
                        if re.search(pattern1,prom) or re.search(pattern2,prom):
                            count += 1
                            if regione.count(i+1) == 0:
                                regione.append(i+1)
                            scale = x2 - j - 1
                            dati = []
                            for g in range(25):
                                w = z1[x1+g:x2+g]
                                at =100*(w.count("A")+w.count("T"))/len(w)
                                if x2 + g + 50 <= len(z1):
                                    prom = z1[x1+g:x2+g+50]
                                else:
                                    prom = z1[x1+g:]
                                if at >= 40:
                                    if re.search(pattern1,prom) or re.search(pattern2,prom):
                                        dati.append(at)
                                    else:
                                        dati.append(0)
                            maxP = np.argmax(dati)
                            at = np.max(dati)
                            x1 = x1 + maxP
                            x2 = x2 + maxP 
                            if x2 + 50 <= len(z1):
                                prom = z1[x1:x2+50]
                                x = x1+u
                                y = x2+50+u
                            else:
                                prom = z1[x1:]
                                x = x1+u
                                y = len(z1)+u
                            q.write("\nP")
                            q.write(str(cod2))
                            q.write("minus             ")
                            q.write(str(i+1))
                            q.write("             negative           ")
                            q.write(str(x))
                            q.write("           ")
                            q.write(str(y))
                            cod2 += 1
                            ATregions.append(at)
                            promotori += 1
                            if re.search(patternB,prom):
                                SI1 += 1
                            else:
                                NO1 += 1
                            if re.search(patternD,prom):
                                SI2 += 1
                            else:
                                NO2 += 1
            if len(z1) >= 50:
                number_predictions.append(count)
        
        
        q.close()
        
        regione.sort()
        with open("ABOUT PROMOTERS (METHOD 1).txt", "a") as file:
            file.write("Genes for which it was predicted a promoter\n\n\n")
            for s in regione:
                file.write("\nGene number ")
                file.write(str(s))
                file.write("\n")
                
            file.write("\n\nMean GC% content of genome sequence: ")
            file.write(str(gc))
            file.write("%")
            r = 100*(len(regione)/intergeniche)
            file.write("\n\n\nNumber of intergenic regions >=50bp in length >= 50 nt: ")
            file.write(str(intergeniche))
            file.write("\n\n\nPercentage of genes in which in the upstream intergenic region>=50bp at least a putative promoter was predicted ")
            file.write(str(r))
            file.write("%")
            file.write("\nThe number of these intergenic regions is ")
            file.write(str(len(regione)))
            file.write("\nTotal number of predicted regions: ")
            file.write(str(promotori))
            file.write("\n\n\n\nPercentage of putative promoters cointaining the pattern TANNNT: ")
            file.write(str(100*SI1/(SI1+NO1)))
            file.write("%\nThe total number of these regions is ")
            file.write(str(SI1))
            file.write("\n\nPercentage of putative promoters cointaining the pattern TTGAC ")
            file.write(str(100*SI2/(SI2+NO2)))
            file.write("%\nThe total number of these regions is ")
            file.write(str(SI2))
            file.write("\n\nMean AT% content of AT-rich element in predicted promoters: ")
            file.write(str(np.mean(ATregions)))
            file.write("%\nStandard deviation: ")
            file.write(str(np.std(ATregions)))
            file.write("\n\n\n")
            for i in range(np.max(number_predictions)+1):
                file.write("Intergenic regions with ")
                file.write(str(i))
                file.write(" predictions: ")
                file.write(str(number_predictions.count(i)))
                file.write("\n")
        
    
        print("Work finished, see output files in the current directory")
    except IOError:   
        print ("File %s inexistent in the current directory!" %(File_geni))   

except IOError:   
    print ("File %s inexistent in the current directory!" %(File_genoma))