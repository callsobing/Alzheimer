# -*- coding: utf-8 -*-
"""
Created on Tue Feb 02 01:23:57 2016

@author: Chester
"""

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
import argparse
import os
import sys
import glob
import numpy as np
import scipy as sp
from scipy.sparse import coo_matrix
from sklearn.utils import shuffle
from sklearn.cross_validation import KFold
from liblinearutil import *

""""""""""""""""""""""""""""""
# set parameters
""""""""""""""""""""""""""""""
str_inputFileName = ""
str_inputFilePath = ""
str_problemDefinition = ""
str_kFoldCV = ""
str_outputFilePath = ""

""""""""""""""""""""""""""""""
# define functions
""""""""""""""""""""""""""""""
### get args from command line
def parseFeatureAndModeling_argParser():
    global str_inputFileName
    global str_inputFilePath
    global str_problemDefinition
    global str_kFoldCV
    global str_outputFilePath
    str_description = "This Script is for Parsing Feature and Modeling. There are two kinds of input for implementing. 1:Input single file; 2:Input a file path of files"
    argparser_thisParser = argparse.ArgumentParser(description=str_description)
    argparser_thisParser.add_argument("-f", "--str_inputFileName", type=str, help="a string of input file name")
    argparser_thisParser.add_argument("-i", "--str_inputFilePath", type=str, help="a file path of input files")
    argparser_thisParser.add_argument("-s", "--int_problemDefinition", type=int, help="a liblinear code of modeling algo. (Logistic:0; SVC:1(default); SVR:11)",default=1)
    argparser_thisParser.add_argument("-k", "--int_kFoldCV", type=int, help="a integer number k of k-fold cross-validation",default=2)
    argparser_thisParser.add_argument("-o", "--str_outputFilePath", default=os.getcwd()+"/", help="a string filepath of output file")
    args = argparser_thisParser.parse_args()
    if not(args.str_inputFileName or args.str_inputFilePath):
        print "Lack of input parameters."
    else:
        if args.str_inputFileName:
            print "A string of input file name has gotten."
            str_inputFileName = args.str_inputFileName
        elif args.str_inputFilePath:
            print "A file path of input files has gotten."
            str_inputFilePath = args.str_inputFilePath
    if args.int_problemDefinition==1:
        print "Use default modeling algorithm (Support Vector Classification)"
        str_problemDefinition = "1"
    elif args.int_problemDefinition==0:
        print "Change modeling algorithm to logistic regression (Classification)"
        str_problemDefinition = "0"
    elif args.int_problemDefinition==11:
        print "Change modeling algorithm to support vector regression (Regression)"
        str_problemDefinition = "11"
    if args.int_kFoldCV>=2:
        print "Use " + str(args.int_kFoldCV) + "-folds cross validation"
        str_kFoldCV = str(args.int_kFoldCV)
    if str(args.str_outputFilePath)!=str(os.getcwd())+"/":
        print "A string filepath of output file has gotten."
        str_outputFilePath = args.str_outputFilePath

### write results to .gen file
def parseFeatureAndModeling_outputResults(str_outputFileName, list_results):
    with open(str_outputFileName, 'w') as file_outputFile:
        for idxl in range(0, len(list_results)):
            file_outputFile.write(list_results[idxl]+"\n")

### detect allel type of subject
def parseFeatureAndModeling_allelTypeParser(list_SNPOfSubject):
    list_allelType = [0,0,0]
    list_allelType[np.argmax(list_SNPOfSubject)] = 1
    
    return list_allelType

### parse each SNP in .gen file
def parseFeatureAndModeling_snpPairParser(np_ADNIClinical_ID, np_ADNIClinical_Phenotype, np_ADNIClinical_Progress, list_inputFileContent):
    np_allelType = np.empty([len(np_ADNIClinical_ID),len(list_inputFileContent)*3],dtype='int')    
    for idxSNP in range(0,len(list_inputFileContent)):
        list_splitInputFileContent = list_inputFileContent[idxSNP].split(" ")
        #list_InputFileContent_ID = list_splitInputFileContent[2]
        list_InputFileContent_SNP = list_splitInputFileContent[6:]
        for idxSubject in range(0,len(np_ADNIClinical_ID)):
            np_allelType[idxSubject,idxSNP*3:idxSNP*3+3] = parseFeatureAndModeling_allelTypeParser([list_InputFileContent_SNP[idxSubject*3],list_InputFileContent_SNP[idxSubject*3+1],list_InputFileContent_SNP[idxSubject*3+2]])
    
    ### select available subject
    list_availableSujectIdx = np.where((np_ADNIClinical_ID=="1"))[0]
    np_allelType = np_allelType[list_availableSujectIdx,:]    
    np_ADNIClinical_Phenotype = np_ADNIClinical_Phenotype[list_availableSujectIdx,:]
    np_ADNIClinical_Progress = np_ADNIClinical_Progress[list_availableSujectIdx]
    
    ### parse label
    np_feature_snp = np_allelType
    np_label_phenotype = np.argmax(np_ADNIClinical_Phenotype, axis=1)
    np_label_progress = np_ADNIClinical_Progress.astype(float)
    #np_label_progress = np.zeros(len(np_ADNIClinical_Progress))
    #list_ProgressInFastSujectIdx = np.where(np_ADNIClinical_Progress.astype(float)<=0)[0]
    #np_label_progress[list_ProgressInFastSujectIdx] = 1
    
    return np_feature_snp,np_label_phenotype,np_label_progress

def parseFeatureAndModeling_classificationCrossValidation(np_X, np_y, str_model, kf):
    try:    
        X_func = np_X
        y_func = np_y
        X_sparse = coo_matrix(X_func)
        X_func, X_sparse, y_func = shuffle(X_func, X_sparse, y_func, random_state=0)
        
        list_target = []
        list_predict = []
           
        for idxTr, idxTe in kf:
            prob  = problem(list(y_func[idxTr]), [list(item) for item in list(X_func[idxTr])])
            param = parameter("-s " + str_model + " -c 1 -q")
            m = train(prob, param)
            label, acc, val = predict(list(y_func[idxTe]), [list(item) for item in list(X_func[idxTe])], m)
            for idxy,idxl in zip(list(y_func[idxTe]),label):
                list_target.append(idxy)
                list_predict.append(idxl) 
        
        float_accuracy = float(len([i for i, j in zip(list_target, list_predict) if i == j]))/len(list_target)
        
        return float_accuracy
    except:
        return ""

def parseFeatureAndModeling_regressionCrossValidation(np_X, np_y, str_model, kf):
    try:
        X_func = np_X
        y_func = np_y
        X_sparse = coo_matrix(X_func)
        X_func, X_sparse, y_func = shuffle(X_func, X_sparse, y_func, random_state=0)
        
        list_target = []
        list_predict = []
           
        list_target = []
        list_predict = []    
        for idxTr, idxTe in kf:
            prob  = problem(list(y_func[idxTr]), [list(item) for item in list(X_func[idxTr])])
            param = parameter("-s " + str_model + " -c 1 -q")
            m = train(prob, param)
            label, acc, val = predict(list(y_func[idxTe]), [list(item) for item in list(X_func[idxTe])], m)
            for idxy,idxl in zip(list(y_func[idxTe]),label):
                list_target.append(float(idxy))
                list_predict.append(idxl)
        
        float_Pearson = sp.stats.stats.pearsonr(list_target,list_predict)[0]
        float_Spearman = sp.stats.stats.spearmanr(list_target,list_predict)[0]
        
        return float_Pearson,float_Spearman
    except:
        return "",""
        

        
""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def main():
    parseFeatureAndModeling_argParser()
    ### get filename clinical data
    str_filename_ADNIClinical = os.getcwd() + "/ADNI/Alzheimer_Data_Train_MergeGenotypeID.csv"

    ### get ADNI file
    global list_ADNIClinical
    list_ADNIClinical = []
    file_ADNIClinical = open(str_filename_ADNIClinical, 'r')
    for line in file_ADNIClinical.readlines():
        list_ADNIClinical.append(line.strip().split(','))
    file_ADNIClinical.close()
    np_ADNIClinical = np.array(list_ADNIClinical)

    global int_lineCount
    int_inputFileSubjectCount = 0
    if str_inputFileName!="":
        ### get number of subject in .gen file
        with open(str_inputFileName, 'r') as file_inputFile:
            for line in file_inputFile:
                int_inputFileSubjectCount = (len(line.strip().split())-6)/3
        ### get .gen file
        list_inputFileContent = []
        with open(str_inputFileName, 'r') as file_inputFile:
            for line in file_inputFile:
                list_inputFileContent.append(line)
        np_feature_snp,np_label_phenotype,np_label_progress = parseFeatureAndModeling_snpPairParser(np_ADNIClinical[:,1], np_ADNIClinical[:,2:6], np_ADNIClinical[:,13],list_inputFileContent)
        kf = KFold(np_feature_snp.shape[0], n_folds=int(str_kFoldCV))
        if (str_problemDefinition=="0" or str_problemDefinition=="1"):
            list_results = ["Accuracy"]
            float_accuracy = parseFeatureAndModeling_classificationCrossValidation(np_feature_snp, np_label_phenotype, str_problemDefinition, kf)
            str_outputFileName = str_inputFileName.split("/")[-1].split(".")[0] + "_s" + str_problemDefinition + "_k" + str_kFoldCV + ".txt"
            list_results.append(str(float_accuracy))
            parseFeatureAndModeling_outputResults(str_outputFilePath+str_outputFileName,list_results)
        elif str_problemDefinition=="11":
            list_results = ["Pearson Cor.,Spearman Cor.,AVG(P+S)"]
            float_pearson,float_spearman = parseFeatureAndModeling_regressionCrossValidation(np_feature_snp, np_label_progress, str_problemDefinition, kf)
            float_avg = (float_pearson+float_spearman)/2
            str_outputFileName = str_inputFileName.split("/")[-1].split(".")[0] + "_s" + str_problemDefinition + "_k" + str_kFoldCV + ".txt"
            list_results.append(str(float_pearson) + "," + str(float_spearman) + "," + str(float_avg))
            parseFeatureAndModeling_outputResults(str_outputFilePath+str_outputFileName,list_results)
        print "Completed"
    
    elif str_inputFilePath!="":
        int_lineCount = 0
        str_outputFileName = str_inputFilePath.split("/")[-2] + "_s" + str_problemDefinition + "_k" + str_kFoldCV + ".csv"
        list_results = []        
        if (str_problemDefinition=="0" or str_problemDefinition=="1"):
            list_results.append("GeneSymbol,Filter,Accuracy")
        elif str_problemDefinition=="11":
            list_results.append("GeneSymbol,Filter,Pearson Cor.,Spearman Cor.,AVG(P+S)")
        
        ### get all gene file names from input file path
        list_fileNameOfGenFile = []
        for str_fileName in os.listdir(str_inputFilePath):
            if str_fileName.endswith(".gen"):
                list_fileNameOfGenFile.append(str_fileName)

        kf = KFold(767, n_folds=int(str_kFoldCV))
        for idxf in range(0,len(list_fileNameOfGenFile)):
            ### get number of subject in .gen file
            with open(str_inputFilePath + str(list_fileNameOfGenFile[idxf]), 'r') as file_inputFile:

                for line in file_inputFile:
                    int_inputFileSubjectCount = (len(line.strip().split())-6)/3
            ### get .gen file
            list_inputFileContent = []
            with open(str_inputFilePath + str(list_fileNameOfGenFile[idxf]), 'r') as file_inputFile:
                for line in file_inputFile:
                    list_inputFileContent.append(line)
                    
            np_feature_snp,np_label_phenotype,np_label_progress = parseFeatureAndModeling_snpPairParser(np_ADNIClinical[:,1], np_ADNIClinical[:,2:6], np_ADNIClinical[:,13],list_inputFileContent)
            
            #kf = KFold(np_feature_snp.shape[0], n_folds=int(str_kFoldCV))
            if np_feature_snp!="":
                if (str_problemDefinition=="0" or str_problemDefinition=="1"):
                    try:                
                        float_accuracy = parseFeatureAndModeling_classificationCrossValidation(np_feature_snp, np_label_phenotype, str_problemDefinition, kf)
                        str_results = str(list_fileNameOfGenFile[idxf]).split("_")[0] + "," + str(list_fileNameOfGenFile[idxf]).split("_")[1] + "," + str(float_accuracy)
                        list_results.append(str_results)
                    except:
                        pass
                elif str_problemDefinition=="11":
                    try:
                        float_pearson,float_spearman = parseFeatureAndModeling_regressionCrossValidation(np_feature_snp, np_label_progress, str_problemDefinition, kf)
                        float_avg = (float_pearson+float_spearman)/2
                        str_results = str(list_fileNameOfGenFile[idxf]).split("_")[0] + "," + str(list_fileNameOfGenFile[idxf]).split("_")[1] + "," + str(float_pearson) + "," + str(float_spearman) + "," + str(float_avg)
                        list_results.append(str_results)
                    except:
                        pass
            
            int_lineCount = int_lineCount + 1
            str_print = "Modeling .gen file: " + "{0:.2f}".format(float(int_lineCount)/len(list_fileNameOfGenFile)*100) + "%"
            sys.stdout.write('%s\r' % str_print)
            sys.stdout.flush()
        parseFeatureAndModeling_outputResults(str_outputFilePath+str_outputFileName,list_results)
                
        

if __name__ == "__main__":
    main()