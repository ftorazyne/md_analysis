#-*- coding: utf-8 -*-
import sys
from matrixops import *

def normalize_data(data,filename):
    def empty(point):
        return point != ' ' and point != '\t' and point != ''
    #Загрузка диапазонов нормировок
    comments = ['#', '%', '$', '@', ';']
    ranges = []
    inf1 = open(filename, 'r')
    lines = inf1.readlines()
    inf1.close()
    for s in lines:
        if ((s[0] not in comments) and (s != '')):
            s = s.strip()
            line = s.split(' ')
            line = filter(empty, line)
            if len(line)<2:
                print("Во файле нормировок нет нужного количества столбцов данных.\n"+\
                      "Нехорошо издеваться над программой, она и так написана на Python")
                raise SystemExit(1)
            ranges.append([float(field.strip()) for field in line])
    #Вычисление параметров нормировочных функций
    functions=[]
    for range in ranges:
        koef=1/(range[1]-range[0])
        const=(range[0]+range[1])*koef
        functions.append((koef,const))
    #Нормировка данных
    data_norm=[]
    for exemple in data:
        data_norm.append([2*koef*value-const for value in exemple])
    return data_norm



def calcaverages(data,data_dimension,public=0):
    value_quantity=len(data)
    if (public==1): sys.stdout.write("Вычисление средних       {0:3d}%".format(0))
    medians = [0.0 for i in range(data_dimension)]
    for i in xrange(0, data_dimension):
        for j in xrange(0, value_quantity):
            medians[i] = medians[i] + data[j][i]
        medians[i] = medians[i] / value_quantity
        if (public==1): sys.stdout.write(" \b\b\b\b\b{0:3d}%".format(100 * (i + 1) / data_dimension))
    if (public==1): sys.stdout.write("  Готово\n")
    return medians

def covariance_matrix(data,medians,data_dimension,filename,public=0):
    #Вычисление дисперсий
    value_quantity=len(data)
    covmatrix = [[0.0 for j in range(data_dimension)] for i in range(data_dimension)]
    if (public==1): sys.stdout.write("Вычисление дисперсий     {0:3d}%".format(0))
    for i in xrange(0, data_dimension):
        for j in xrange(0, value_quantity):
            covmatrix[i][i] = covmatrix[i][i] + (data[j][i] - medians[i]) ** 2
        covmatrix[i][i] = covmatrix[i][i] / (value_quantity - 1)
        if (public==1): sys.stdout.write(" \b\b\b\b\b{0:3d}%".format(100 * (i + 1) / data_dimension))
    if (public==1): sys.stdout.write("  Готово\n")
    #Вычисление матрицы ковариаций
    if (public==1): sys.stdout.write("Вычисление ковариаций    {0:3d}%".format(0))
    counter = 0
    for i in xrange(1, data_dimension):
        for j in xrange(0, i):
            current_covariance = 0.0
            for k in xrange(0, value_quantity):
                current_covariance = current_covariance + (data[k][i] - medians[i]) * (data[k][j] - medians[j])
            current_covariance = current_covariance / value_quantity
            covmatrix[i][j] = current_covariance
            covmatrix[j][i] = current_covariance
            counter += 1
        sys.stdout.write(" \b\b\b\b\b{0:3d}%".format(counter * 200 / (data_dimension * (data_dimension - 1))))
    if (public==1): sys.stdout.write("  Готово\n")
    #Вывод матрицы ковариаций
    covf=open(filename+".covariance","w")
    covf.write(matrixtostring(covmatrix,prec=16,field=24))
    covf.close()
    return covmatrix

def load_covmatrix(filename):
    covfile=open(filename,"r")
    covmatrix_lines=covfile.readlines()
    covfile.close()
    return readmatrix(covmatrix_lines)

def load_qrdecomp(filename):
    qrfile=open(filename,'r')
    qr_lines=qrfile.readlines()
    qrfile.close()
    qr_raw=readmatrix(qr_lines)
    qr_dim=len(qr_raw)
    return [[qr_raw[i] for i in range(0,qr_dim/3)],[qr_raw[i] for i in range(qr_dim/3,2*qr_dim/3)],[qr_raw[i] for i in range(2*qr_dim/3,qr_dim)]]

def load_eigenvalues(filename):
    eigvaluefile=open(filename,"r")
    eigenvalues_lines=eigvaluefile.readlines()
    eigvaluefile.close()
    return readmatrix(eigenvalues_lines)

def load_eigenvectors(filename):
    eigvaluefile=open(filename,"r")
    eigenvalues_lines=eigvaluefile.readlines()
    eigvaluefile.close()
    return readmatrix(eigenvalues_lines)