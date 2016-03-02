# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 10:34:26 2016

@author: Chester

Note: You should download the ADNI dataset & the HapmapPhased database at first!!!
"""

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
import argparse
import os
import sys
from MySQLdb.constants import FIELD_TYPE
import _mysql

""""""""""""""""""""""""""""""
# set parameters
""""""""""""""""""""""""""""""
str_inputGeneName = ""
str_inputGeneChromosome = ""
str_inputGeneStartLoci = ""
str_inputGeneEndLoci = ""
str_inputFileName = ""
str_outputFilePath = ""
db = None

""""""""""""""""""""""""""""""
# define functions
""""""""""""""""""""""""""""""
### get args from command line
def queryUCSC_argParser():
    global str_inputGeneName
    global str_inputGeneChromosome
    global str_inputGeneStartLoci
    global str_inputGeneEndLoci
    global str_inputFileName
    global str_outputFilePath
    str_description = "This Script is for Query UCSC Database. There are three kinds of input for query. 1:GeneSymbol; 2:User-defined segment; 3:String filename of gene list."
    argparser_thisParser = argparse.ArgumentParser(description=str_description)
    argparser_thisParser.add_argument("-g", "--str_inputGeneName", type=str, help="a string of gene name")
    argparser_thisParser.add_argument("-c", "--int_inputGeneChromosome", type=int, help="chromosome num of user-defined segment")
    argparser_thisParser.add_argument("-s", "--int_inputGeneStartLoci", type=int, help="a integer loci of user-defined segment's start site")
    argparser_thisParser.add_argument("-e", "--int_inputGeneEndLoci", type=int, help="a integer loci of user-defined segment's end site")
    argparser_thisParser.add_argument("-i", "--str_inputFileName", type=str, help="a string filename of input file of gene list")
    argparser_thisParser.add_argument("-o", "--str_outputFilePath", default=os.getcwd()+"/", help="a string filepath of output file")
    args = argparser_thisParser.parse_args()
    if not(args.str_inputGeneName or (args.int_inputGeneChromosome and args.int_inputGeneStartLoci and args.int_inputGeneEndLoci) or args.str_inputFileName):
        print "Lack of input parameters."
    else:
        if args.str_inputGeneName:
            print "A string of gene name has gotten."
            str_inputGeneName = args.str_inputGeneName
        elif (args.int_inputGeneChromosome and args.int_inputGeneStartLoci and args.int_inputGeneEndLoci):
            print "A user-defined segment has gotten."
            str_inputGeneChromosome = str(args.int_inputGeneChromosome)
            str_inputGeneStartLoci = str(args.int_inputGeneStartLoci)
            str_inputGeneEndLoci = str(args.int_inputGeneEndLoci)
        elif args.str_inputFileName:
            print "A string filename of input file of gene list has gotten."
            str_inputFileName = args.str_inputFileName
    if str(args.str_outputFilePath)!=str(os.getcwd())+"/":
        print "A string filepath of output file has gotten."
        str_outputFilePath = args.str_outputFilePath

### write results to .gen file
def queryUCSC_outputResults(str_outputFileName, list_outputSNPs):
    with open(str_outputFileName, 'w') as file_outputFile:
        for idxl in range(0, len(list_outputSNPs)):
            file_outputFile.write(list_outputSNPs[idxl])

### query UCSC database by sql command
def queryUCSC_mySQLConnector(str_hgVersion, str_geneSymbol="", str_rsIDOfSNP="", str_chromOfSNP=""):
    global db
    ### db is global to prevent reconnecting.
    if db is None:
        print "SQL connection between UCSC database is re-opened. (" + str(str_hgVersion) + ")" 
        conv = { FIELD_TYPE.LONG: int }
        db = _mysql.connect(host='genome-mysql.cse.ucsc.edu',user='genome',passwd='',db=str_hgVersion,conv=conv)
    if str_geneSymbol!="":
        db.query("""SELECT * FROM kgXref INNER JOIN knownGene ON kgXref.kgID=knownGene.name WHERE kgXref.geneSymbol = '%s'""" % str_geneSymbol)
        r = db.use_result().fetch_row(how=1,maxrows=0)
        try:
            if len(r)>=1:
                return r[0]['geneSymbol'], r[0]['chrom'], r[0]['txStart'], r[0]['txEnd'], r[0]['strand']
            else:
                print str(str_geneSymbol) + " miss in UCSC database."
                return ""
        except:
            pass
    elif (str_rsIDOfSNP!="" and str_chromOfSNP!=""):
        db.query("SELECT name, func FROM snp130 WHERE name='" + str_rsIDOfSNP + "' AND chrom='" + str_chromOfSNP + "' AND func REGEXP 'nonsense|missense|frameshift' LIMIT 1")
        r = db.use_result().fetch_row(how=1,maxrows=0)
        try:
            if len(r)>=1:
                return True
            else:
                return False
        except:
            pass
    else:
        print "SQL connection between UCSC database is abandoned"
        return

### search snp in gene from local HapmapPhased database
def queryUCSC_ADNIParser(tuple_thisUCSCDBGene,list_HapmapPhasedSNPInGene=""):
    list_outputADNISNPs = []
    list_outputDBSNPSNPs = []
    list_outputHapMapSNPs = []
    list_outputPromoterSNPs = []
    dict_HapmapPhasedSNPInGene = {}
    if list_HapmapPhasedSNPInGene!="":
        list_HapmapPhasedSNPInGene = str(list_HapmapPhasedSNPInGene[0]).split()
        list_HapmapPhasedSNPInGene = list_HapmapPhasedSNPInGene[1].split(":")[1].split(",")
        for item in list_HapmapPhasedSNPInGene:
            dict_HapmapPhasedSNPInGene[item] = 1
    try:
        if tuple_thisUCSCDBGene != "":
            ### get number of SNP from .gen file
            int_snpCount = 1
            with open(os.getcwd() + "/ADNI/ADNI_merged_" + str(tuple_thisUCSCDBGene[1]).replace("chr","") + ".gen", 'r') as file_genFile:
                int_snpCount = sum(1 for _ in file_genFile)
            with open(os.getcwd() + "/ADNI/ADNI_merged_" + str(tuple_thisUCSCDBGene[1]).replace("chr","") + ".gen", 'r') as file_ADNISNP:
                int_lineCount = 0            
                for line in file_ADNISNP:
                    list_thisSNP = line.strip().split()
                    if int(list_thisSNP[3]) >= int(tuple_thisUCSCDBGene[2]) and int(list_thisSNP[3]) <= int(tuple_thisUCSCDBGene[3]):
                        ### ADNI dataset
                        list_outputADNISNPs.append(line)
                        ### ADNI dataset with dbSNP filter
                        try:
                            if(queryUCSC_mySQLConnector("hg18", str_rsIDOfSNP=str(list_thisSNP[2]), str_chromOfSNP=str(tuple_thisUCSCDBGene[1]))):
                                list_outputDBSNPSNPs.append(line)
                        except:
                            global db
                            db = None
                            if(queryUCSC_mySQLConnector("hg18", str_rsIDOfSNP=str(list_thisSNP[2]), str_chromOfSNP=str(tuple_thisUCSCDBGene[1]))):
                                list_outputDBSNPSNPs.append(line)
                        ### ADNI dataset with dbSNP filter
                        if list_HapmapPhasedSNPInGene!="":
                            try: 
                                if (dict_HapmapPhasedSNPInGene[str(list_thisSNP[2])]==1):
                                    list_outputHapMapSNPs.append(line)
                            except:
                                pass
                    if str(tuple_thisUCSCDBGene[4])=='+' and int(list_thisSNP[3]) >= int(tuple_thisUCSCDBGene[2])-2000 and int(list_thisSNP[3]) <= int(tuple_thisUCSCDBGene[2]):
                        list_outputPromoterSNPs.append(line)
                    elif str(tuple_thisUCSCDBGene[4])=='-' and int(list_thisSNP[3]) >= int(tuple_thisUCSCDBGene[3]) and int(list_thisSNP[3]) <= int(tuple_thisUCSCDBGene[3])+2000:
                        list_outputPromoterSNPs.append(line)
                    int_lineCount = int_lineCount + 1
                    str_print = "Parse .gen file: " + str(tuple_thisUCSCDBGene[0]) + " - " + "{0:.2f}".format(float(int_lineCount)/int_snpCount*100) + "%"
                    sys.stdout.write('%s\r' % str_print)
                    sys.stdout.flush()
                print ("")
                queryUCSC_outputResults(str_outputFilePath + str(tuple_thisUCSCDBGene[0]) + "_ADNI_" + str(len(list_outputADNISNPs)) + "_" + str(len(list_outputADNISNPs)) + ".gen",list_outputADNISNPs)
                queryUCSC_outputResults(str_outputFilePath + str(tuple_thisUCSCDBGene[0]) + "_dbSNP_" + str(len(list_outputDBSNPSNPs)) + "_" + str(len(list_outputADNISNPs)) + ".gen",list_outputDBSNPSNPs)
                if list_HapmapPhasedSNPInGene!="":
                    queryUCSC_outputResults(str_outputFilePath + str(tuple_thisUCSCDBGene[0]) + "_HapMapPhased_"+ str(len(list_outputHapMapSNPs)) + "_" + str(len(list_outputADNISNPs)) +".gen",list_outputHapMapSNPs)
                queryUCSC_outputResults(str_outputFilePath + str(tuple_thisUCSCDBGene[0]) + "_promoter_" + str(len(list_outputPromoterSNPs)) + "_" + str(len(list_outputADNISNPs)) + ".gen",list_outputPromoterSNPs)
    except:
        pass

### search snp in gene from local HapmapPhased database
def queryUCSC_hapmapPhasedSearcher(tuple_thisUCSCDBGene):
    list_HapmapPhasedSNP = []    
    str_HapmapPhasedSNPInGene = ""
    list_HapmapPhasedSNPInGene = []
    try:
        if tuple_thisUCSCDBGene != "":
            with open(os.getcwd() + "/HapmapPhased/Phased3/hapmap3_r2_b36_fwd.consensus.qc.poly." + tuple_thisUCSCDBGene[1] + "_ceu.unr.phased.gz", 'r') as file_HapmapPhasedSNP:
            #with open(os.getcwd() + "/HapmapPhased/Phased2/genotypes_" + tuple_thisUCSCDBGene[1] + "_CEU_r21_nr_fwd_legend.txt.gz", 'r') as file_HapmapPhasedSNP:
                for line in file_HapmapPhasedSNP:
                    list_thisSNP = line.strip().split()
                    list_HapmapPhasedSNP.append([list_thisSNP[0],list_thisSNP[1]])
            list_HapmapPhasedSNP.pop(0)
            for item in list_HapmapPhasedSNP:
                if int(str(item[1]).strip()) > int(tuple_thisUCSCDBGene[2]) and int(str(item[1]).strip()) < int(tuple_thisUCSCDBGene[3]):
                    str_HapmapPhasedSNPInGene = str_HapmapPhasedSNPInGene + str(item[0]) + ","
            str_HapmapPhasedSNPInGene = str_HapmapPhasedSNPInGene[:-1]
            list_HapmapPhasedSNPInGene.append(str(tuple_thisUCSCDBGene[0])+" "+ str(tuple_thisUCSCDBGene[1]) + ":" + str_HapmapPhasedSNPInGene)
        return list_HapmapPhasedSNPInGene
    except:
        return ""

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def main():
    global db
    queryUCSC_argParser()
    list_HapmapPhasedSNPInGene = []
    if str_inputGeneName!="":
        db = None        
        tuple_UCSCDBGene_hg19 = queryUCSC_mySQLConnector("hg19",str_geneSymbol=str_inputGeneName)
        db = None
        tuple_UCSCDBGene_hg18 = queryUCSC_mySQLConnector("hg18",str_geneSymbol=str_inputGeneName)
        if tuple_UCSCDBGene_hg19!="" and tuple_UCSCDBGene_hg18!="":
            list_HapmapPhasedSNPInGene = queryUCSC_hapmapPhasedSearcher(tuple_UCSCDBGene_hg18)
            queryUCSC_ADNIParser(tuple_UCSCDBGene_hg19, list_HapmapPhasedSNPInGene)
    elif (str_inputGeneChromosome!="" and str_inputGeneStartLoci!="" and str_inputGeneEndLoci!=""):
        tuple_UCSCDBGene_hg19 = ("User-Defined-Segment","chr"+str_inputGeneChromosome,str_inputGeneStartLoci,str_inputGeneEndLoci)
        queryUCSC_ADNIParser(tuple_UCSCDBGene_hg19)
    elif str_inputFileName!="":
        list_inputGeneSymbol = []
        with open(str_inputFileName, 'r') as file_inputFileName:
            for line in file_inputFileName:
                list_inputGeneSymbol.append(line.strip())
        for gene in list_inputGeneSymbol:
            db = None
            tuple_UCSCDBGene_hg19 = queryUCSC_mySQLConnector("hg19",str_geneSymbol=gene)
            db = None
            tuple_UCSCDBGene_hg18 = queryUCSC_mySQLConnector("hg18",str_geneSymbol=gene)
            if tuple_UCSCDBGene_hg19!="" and tuple_UCSCDBGene_hg18!="":
                list_HapmapPhasedSNPInGene = queryUCSC_hapmapPhasedSearcher(tuple_UCSCDBGene_hg18)
                queryUCSC_ADNIParser(tuple_UCSCDBGene_hg19, list_HapmapPhasedSNPInGene)
    
if __name__ == "__main__":
    main()