#!/usr/bin/python3

"""
#====================================================================================
#
#          FILE: cisVar (python 3)
#        AUTHOR: Ashley Tehranchi, tehranchi576@gmail.com 
#  ORGANIZATION: Stanford University
#       CREATED: 2015-12-12 11:33
# Last modified: 2016-11-11 10:56
#
#
#====================================================================================

######################################################################################
genotypes file format as follows, sorted by individuals(columns) and by snps (rows)

chr   position    ref alt NA18486 NA18489 ...
chr1  10583   g   a   0   1 ...
chr1  10611   c   g   0   0 ...
chr1  13302   c   t   1   0 ... 

######################################################################################

cisVar mpileup -F <SampleName> -f <fastaFile> -p <mpileupBEDfile> -B <sortedBam>
cisVar post -F <SampleName> -r <readDepth> -a <allelesFile>
cisVar geno -F <SampleName> -r <readDepth> -i <individualsFile> -g <genotypesFile>
cisVar qtls -F <SampleName> -r <readDepth> -n <numberIndividuals>

"""

######################################################################################
# ARGUEMENTS
######################################################################################


import argparse
import sys
import os
import operator
import subprocess
import scipy
import numpy as np
from scipy import stats


parser = argparse.ArgumentParser(description='Find cis QTLs based on an experimental selection method')
parser.add_argument('inputCommand', help='<mpileup>  <post>  <geno>  <qtls> ')
parser.add_argument('-F', '--SampleName', required=False, dest='prefix_name', help='sample/population name')
parser.add_argument('-f', '--fasta', required=False, dest='allCHRfasta', help='fasta file with all chromosomes')
parser.add_argument('-B', '--BAMfile', required=False, dest='sortedBam', help='sorted BAM file')
parser.add_argument('-r', '--readDepth', required=False, type=int, dest='trial_depths', help='minimum read depth per variant')
parser.add_argument('-i', '--individualsFile', required=False, dest='individualslist', help='list of individuals matching genotype matrix; one indv per line')
parser.add_argument('-a', '--allelesFile', required=False, dest='allelesFileName', help='format: chr1    10583   G   A')
parser.add_argument('-n', '--numberIndividuals', required=False, dest='numIndv', help='The number of individuals in the pool')
parser.add_argument('-p', '--mpileupBEDfile', required=False, dest='mpileupBEDfile', help='The  mpileup BED file')
parser.add_argument('-g', '--genoFile', required=False, dest='genosFile', help='The  genotypes file')
args = parser.parse_args()


WorkingDirectory = os.getcwd() 

######################################################################################
# MPILEUP FROM BAM FILES USING HG19 MASKED GENOME WITH BLACKLISTED REGIONS REMOVED
######################################################################################

def mpileup(allCHRfasta, mpileupBEDfile, sortedBam, prefix_name):
    os.system("samtools mpileup -Q 25 -f " + allCHRfasta +" -l "+ mpileupBEDfile +" " + sortedBam + " > " + prefix_name + ".mpileup.txt") 
    ## outputs prefix.mpileup.txt
    return()


######################################################################################
#POST-CHIP CALCULATION AND FILE TRIMMIMG PLUS HEADER ADDITION
######################################################################################

def postcalc(prefix_name, trial_depths, allelesFileName):
    os.system("module load samtools")
    print (prefix_name)
    print ("Post-ChIP calculation for min reads ",trial_depths)
    infilename = prefix_name + ".mpileup.txt"
    inalleles = allelesFileName
    outfilename = prefix_name + "." + str(trial_depths) + ".POST.txt"
    depth = int(trial_depths)
    alleles = open(inalleles, 'r')
    pileup = open(infilename, 'r')

    pileup_dictionary = {}
    count = 1
    outfile = open(outfilename, 'w')
    with open(infilename, 'r') as pileup:
        for line in pileup:
            #print count
            count = count + 1
            rows = line.rstrip('\n').split("\t")  
            if (int(rows[3]) < depth):
                Acount = 'NA'
                Ccount = 'NA'
                Gcount = 'NA'
                Tcount = 'NA'
                ALTallele = 'NA'
                ALTdepth = 'NA'
                REFdepth = 'NA'
                RefFreq = 'NA'
                AltFreq = 'NA'
            else:
                if (int(rows[3]) >= depth):
                    read = str(rows[4]).lower()                              # isolate overlapping positions from reads, make all letters lowercase to count them properly
                    Acount = read.count('a')            # count all A's, A is plus strand, a in negative strand, count both
                    Ccount = read.count('c')
                    Gcount = read.count('g')
                    Tcount = read.count('t')
                    countslist = [Acount, Ccount, Gcount, Tcount]       #make list of all counts
                    index, ALTdepth = max(enumerate(countslist), key=operator.itemgetter(1))    #find index of largest number and the value
                    #print("index ", index)
                    indels = rows[4].count('-') + rows[4].count('+')
                    REFdepth = int(rows[3]) - int(Acount) - int(Ccount) - int(Gcount) - int(Tcount) - int(indels)
                    #print("indels, REFdepth ", indels, REFdepth)
                    ##added 9/2/14
                    #print(rows)
                    rows[3] = str(int(REFdepth) + int(ALTdepth)) ## Adjusts totaldepth to remove third alleles/sequencing errors
                    #print(rows)
                    ##
                    ALTallele = "NA"
                    if (index == 0):
                        ALTallele = 'A'
                    if (index == 1):
                        ALTallele = 'C'
                    if (index == 2):
                        ALTallele = 'G'
                    if (index == 3):
                        ALTallele = 'T'
                    RefFreq = "NA"
                    AltFreq = "NA"
                    if (rows[3] == 0):
                        RefFreq = str(0)
                        print ("mistake")
                    else:
                        if (float(rows[3]) > depth) and ((float(REFdepth) / float(rows[3])) >= 0) :    ## depth greater than 0 and postfreq > 0 
                            RefFreq = str(REFdepth / float(rows[3]))            # to get ref freq, totaldepth - altdepth = refdepth / totaldepth
                            AltFreq = str(1 - float(RefFreq))
            pileup_dictionary[rows[0] + "." + rows[1]] = rows[2:] + [Acount, Ccount, Gcount, Tcount, ALTdepth, REFdepth, ALTallele, RefFreq, AltFreq]
            #print("ALTallele, row ",ALTallele, rows)
            assert len(pileup_dictionary[rows[0] + "." + rows[1]]) == 13
    linecounter =1  
    with open(inalleles) as alleles:
        for line in alleles:
            row = line.rstrip().split('\t')
            genoRefAllele = row[2]
            genoAltAllele = row[3]
            try: matchingVector = pileup_dictionary[row[0] + '.' + row[1]]
            except: continue
            matchingVector = [str(x) for x in matchingVector]
            if (matchingVector[8] == str(0)):                 #for alt alleles with 0 reads, force alt allele to be imputed geno alt allele
                matchingVector[10] = genoAltAllele
            if (matchingVector[9] == str(0)):                     #for ref alleles with 0 reads, force ref allele to be imputed geno ref allele  
                matchingVector[0] = genoRefAllele
            postChipRefAllele = matchingVector[0].upper()
            postChipRefFreq = matchingVector[11]
            postChipAltAllele = matchingVector[10].upper()
            postChipAltFreq = matchingVector[12]
            if (postChipRefAllele == genoRefAllele) and (postChipAltAllele == genoAltAllele):
                outAllele = postChipRefAllele
                outfreq = postChipRefFreq
            elif postChipRefAllele == genoAltAllele and (postChipAltAllele == genoRefAllele):
                outAllele = postChipAltAllele
                outfreq = postChipAltFreq
                linecounter = linecounter +1
            else:
                outAllele = 'NA'
                outfreq = 'NA'
            finalout = row[0] + "\t" + row[1] + "\t" + '\t'.join(matchingVector) + "\t" + str(outAllele) + "\t" + str(outfreq) + "\n"
            outfile.write(finalout)
    outfile.close()
    alleles.close()
    print("post-frequency calculation DONE")
    return ()    


######################################################################################
# POST.TXT FILE TRIMMING, ADDITION OF HEADER, BED FILE CREATION
######################################################################################

def postTrim(prefix_name, trial_depths):
    headerOUT = open("mpileup.header.txt", 'w')
    headertemp = 'Chr'+'\t'+'position'+'\t'+'REFallele'+'\t'+'Depth'+'\t'+'Acount'+'\t'+'Ccount'+'\t'+'Gcount'+'\t'+'Tcount'+'\t'+'ALTdepth'+'\t'+'REFDepth'+'\t'+'ALTallele'+'\t'+'REFfreq'+'\t'+'ALTfreq'+'\t'+'POSTallele'+'\t'+'POSTfreq'+'\n'
    headerOUT.write(headertemp)
    headerOUT.close()

    print ("File trimming ")
    os.system("sed '/NA/d' " + prefix_name + "." + str(trial_depths) +  ".POST.txt | cut -f-4,7- > " + prefix_name + "." + str(trial_depths) + ".POSTt.txt") 
    os.system("cat mpileup.header.txt " + prefix_name + "."  + str(trial_depths) + ".POSTt.txt > " + prefix_name + "." + str(trial_depths) + ".POSTth.txt")
    print (os.system("head " + prefix_name + "." + str(trial_depths) + ".POSTth.txt"))
    bedOUT = open(prefix_name + "." + str(trial_depths) + ".bed", 'w')
    openfile = open(prefix_name + "." + str(trial_depths) + ".POSTt.txt", 'r')
    with open(prefix_name + "." + str(trial_depths) + ".POSTt.txt", 'r') as openfile:
        for line in openfile:
            rowtmp = line.rstrip().split('\t')
            chrom = rowtmp[0]
            snp = rowtmp[1]
            startpos = int(snp) - 1
            OUTrow = str(chrom) + '\t' + str(startpos) + '\t' + str(snp) + '\n'
            bedOUT.write(OUTrow)
    bedOUT.close()
    openfile.close()
    print ("File trimming DONE")
    return()

######################################################################################
# GENOTYPE EXTRACTION OF POSTCALC SNPS
######################################################################################
              
def genoExtract(prefix_name, trial_depths, individualslist, genosFile):
    print ("Extracting genotypes ")
    new_prefix = prefix_name + "." + str(trial_depths)
    postdataIN = new_prefix + ".POSTth.txt"
    outName = new_prefix + ".genotypes.txt"
    print("file: ", new_prefix)
    postdata = open(postdataIN, 'r')
    postlocs_dict = {}
    count = 1
    with open(postdataIN, 'r') as postdata:
        for line in postdata:
            #print 'line # ', count
            count = count + 1
            rows = line.rstrip('\n').split('\t')
            postlocs_dict[rows[0] + '.' + rows[1]] = rows[2:]
    filename = genosFile
    genosIN = open(filename, 'r')
    outfile = open(prefix_name+"geno.txt", 'w')
    with open(filename) as genosIN:
        for line in genosIN:
            matchingVector = 'NA'
            row = line.rstrip().split('\t')
            namestring = row[0] + '.' + row[1]
            if (namestring in postlocs_dict):
                matchingVector = row
                matchingVector = [str(x) for x in matchingVector]
                fileOUT = '\t'.join(matchingVector) + '\n'
                outfile.write(fileOUT)
    genosIN.close()
    outfile.close()
    print ("Extraction done")
    print ("transpose geno file with snps matching the postchip data")
    
    filename = prefix_name+"geno.txt"
    tempgeno = open(filename, "r")
    total_list = list()  #Big list for the entire file
    flag=0
    counter =0
    for line in tempgeno:
        counter +=1
        line_list = line.rstrip().split("\t") #line_list holds the elements of one line
        if flag==0:
            lengthLineCheck = len(line_list)
            flag +=1
        if len(line_list)!=lengthLineCheck:
            print (line)
            print ("** WARNING ** Line has weird number of elements at line ", counter) 
            print ("Number of elements: ", len(line_list)) 
            print ("expected number : ", lengthLineCheck) 
        total_list.append(line_list) 
    tempgeno.close()
    
    line_of_interest = total_list[0]
    length_of_line = len(line_of_interest)
    length_of_file = len(total_list)
    print ("File Read In complete")
    #print ("Length of line is ", length_of_line)
    #print ("Length of file is ", length_of_file)
    fileOUT = open(prefix_name+"transpose.txt", 'w')
    for x in range(length_of_line):
        #print "Working on line: " , x , " of ", length_of_line
        val_list = [i[x] for i in total_list]   #for i in total_list, keep the xth element in i
        out_line = '\t'.join(val_list)  + "\n"
        fileOUT.write(out_line)
    fileOUT.close()
    
    # read in individuals
    indv = open(individualslist , 'r')
    indv_list= list()
    print ("Reading file Locations")
    for line in indv:
        row = line.rstrip().split("\t")
        indv_list.append(row)
    print ("Individual Locations Read Complete")
    genosIN1 = open(prefix_name+"transpose.txt", 'r')
    genos_dict = {}
    count = 1
    with open(prefix_name+"transpose.txt", 'r') as genosIN1:
        for line in genosIN1:
            #print count
            count = count + 1
            rows = line.rstrip('\n').split('\t')    #strips \n spaces and splits column by tabs
            genos_dict[rows[0]] = rows[1:]

    # subset pooled individuals from total list of all individuals
    print ("Extracting individuals")
    indv = (individualslist , 'r')
    outfinal = open(outName, 'w')
    with open(individualslist) as indv:
        for line in indv:
            matchingVector = 'NA'
            indrow = line.rstrip().split('\t')
            try: matchingVector = genos_dict[indrow[0]]
            except: continue
            matchingVector = [str(x) for x in matchingVector]
            #print (line)
            fileOUT2 = '\t'.join(matchingVector) + '\n'
            outfinal.write(fileOUT2)
    outfinal.close()
    genos_dict.clear()
    os.system("rm " +prefix_name+ "geno.txt " +prefix_name+ "transpose.txt")
    print ("Extracting genotypes DONE")
    return()


######################################################################################
# REGRESSION WITH Z-SCORE, P-VALUE
######################################################################################

def regPVAL(prefix_name, trial_depths, numIndv):
    new_prefix = prefix_name + "." + str(trial_depths) 
    postdataIN = new_prefix + ".POSTth.txt"
    genosoutName = new_prefix + ".genotypes.txt"
### make sure location of R script is correct    
    commandLine = "sed -e 's!postdataIN!" + postdataIN + "!g'" + " regression_qtls.R | sed -e 's!genosoutName!" + genosoutName + "!g' | sed -e 's!prefix!" + new_prefix + "!g' | sed -e 's!numIndiv!" + numIndv + "!g'> " +prefix_name+ "1.R"
    os.system(commandLine)
    commandout = "Rscript " +prefix_name+ "1.R"
    os.system("chmod u+x " +prefix_name+ "1.R")
    os.system(commandout)
    os.system("rm " +prefix_name+ "1.R")
    print("regression DONE")
    return()


##################
# COMMANDS #
##################

if args.inputCommand == 'mpileup':
    mpileup(args.allCHRfasta, args.mpileupBEDfile, args.sortedBam, args.prefix_name)
elif args.inputCommand == 'post':
    postcalc(args.prefix_name, args.trial_depths, args.allelesFileName)
    postTrim(args.prefix_name, args.trial_depths)
elif args.inputCommand == 'geno':
    genoExtract(args.prefix_name, args.trial_depths, args.individualslist, args.pop_name)
elif args.inputCommand == 'qtls':
    regPVAL(args.prefix_name, args.trial_depths, args.numIndv)
