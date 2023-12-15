#######This script will use OLINK_proteomics file from UK-biobank data and transpose in matrix in way #that PID (patient ID) will be in first colum and the protein id will be in second colums
#Olink_proteomics_data_transposed.txt. #This script will also use the UKB data code 143 and decode #the protein to the UNIPORT name
#Firoj Mahmud, Post-doc karolinska Institute, firoj.mahmud@ki.se


import re

fl = open("olink_testFile.txt", "r") ## input file for reading downloaded from UKbiobank

Dpid = {}
pid = {}

while True:
    line = fl.readline()
        if not line:
        break
        if re.search("eid", line):
        continue
        line = line.rstrip()
    arr = line.split("\t")
    try:
        Dpid[arr[2]]
    except KeyError:
        Dpid[arr[2]] = {}
    
    Dpid[arr[2]][arr[0]] = arr[3]
    pid[arr[0]] = 1

fl.close()
pid = list(pid.keys())
with open("Olink_proteomics_data_transposed.txt", "w") as f:
    f.write("PID\t")
    for i in Dpid.keys():
        f.write(i + "\t")
    f.write("\n")
        for i in range(len(pid)):
        f.write(pid[i] + "\t")
        for j in Dpid.keys():
            f.write(Dpid[j].get(pid[i], "") + "\t")
        f.write("\n")

print("Proteomics data has been processed and written to 'Olink_proteomics_data_transposed.txt'.")

##########################################################################
#####
##This part of the code require the data code 143 from UKBB
#
################################################################################

dictionary = {}
with open("coding143_Olinks.tsv", "r") as f: ###code143 from UKBB
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            dictionary[parts[0]] = parts[1].split(";")[0]
with open("Olink_proteomics_data_transposed2UNIPORT.txt", "w") as f2:
    with open("Olink_proteomics_data_transposed.txt", "r") as f1:
        for line in f1:
            parts = line.strip().split("\t")
                   if parts[1] in dictionary:
                   parts[1] = dictionary[parts[1]]
                        f2.write("\t".join(parts) + "\n")
print("Protein names have been decoded and written to 'Olink_proteomics_data_transposed2UNIPORT.txt'.")
