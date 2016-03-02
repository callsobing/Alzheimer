# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 15:23:08 2016

@author: Chester
"""

import os
import sys

#str_inputFilePath = "D:\\Phd\\Grade_04\\Alzheimer\\Log\\0223\\SelectedGenes\\SelectedGenes\\"
str_inputFilePath = "/home/chester/AD/Alzheimer/Log/0204/"
list_fileNameOfGenFile = []
for idx in range(1,23):
    for str_fileName in os.listdir(str_inputFilePath+"chr"+str(idx)+"/"):
        if "dbSNP" in str_fileName:
            list_fileNameOfGenFile.append("chr"+str(idx)+"/"+str_fileName)

print len(list_fileNameOfGenFile)
int_lineCount = 0
with open(str_inputFilePath+'WholeGenome_dbSNP.gen', 'wb') as file_outputFile:
    for idx in range(0,len(list_fileNameOfGenFile)):
        with open(str_inputFilePath+list_fileNameOfGenFile[idx], 'rb') as file_inputFile:
            file_outputFile.write(file_inputFile.read())
            int_lineCount = int_lineCount + 1
            str_print = "Modeling .gen file: " + "{0:.2f}".format(float(int_lineCount)/len(list_fileNameOfGenFile)*100) + "%"
            sys.stdout.write('%s\r' % str_print)
            sys.stdout.flush()