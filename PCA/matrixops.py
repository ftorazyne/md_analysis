#-*- coding: utf-8 -*-

from math import *
from copy import *
import decimal
import sys


def readmatrix(inp):
    dim = len(inp)
    matrix = []
    for i in range(dim):
        matrix.append([])
        string = inp[i].split()
        for s in string:
            matrix[i].append(float(s))
    return matrix


def multiplymatrixx(matr1, matr2):
    result = []
    dim = len(matr1)
    try:
        for i in range(dim):
            result.append([])
            for j in range(dim):
                cij = 0.0
                for k in range(dim):
                    cij = cij + matr1[i][k] * matr2[k][j]
                result[i].append(cij)
    except IndexError:
        print "Dimensions of matrix are not equal!"
        raise SystemExit(1)
    return result


def multiplymatrix(inmatr1, inmatr2, precise=0):
    result = []
    dim1h = len(inmatr1[0])
    dim1v = len(inmatr1)
    dim2h = len(inmatr2[0])
    if (precise == 1):
        matr1 = [[decimal.Decimal(c) for c in l] for l in inmatr1]
        matr2 = [[decimal.Decimal(c) for c in l] for l in inmatr2]
        try:
            for i in range(dim1v):
                result.append([])
                for j in range(dim2h):
                    cij = decimal.Decimal(0.0)
                    for k in range(dim1h):
                        cij = cij + matr1[i][k] * matr2[k][j]
                    result[i].append(cij)
        except IndexError:
            print "Incorrect dimensions of matrices!"
            raise SystemExit(1)
    else:
        matr1 = inmatr1
        matr2 = inmatr2
        try:
            for i in range(dim1v):
                result.append([])
                for j in range(dim2h):
                    cij = 0.0
                    for k in range(dim1h):
                        cij = cij + matr1[i][k] * matr2[k][j]
                    result[i].append(cij)
        except IndexError:
            print "Incorrect dimensions of matrices!"
            raise SystemExit(1)
    return result


def lincombmatrix(inmatr1, inmatr2, k1=1.0, k2=1.0, precise=0):
    result = []
    dimv = len(inmatr1)
    dimh = len(inmatr1[0])
    if (precise == 1):
        matr1 = [[decimal.Decimal(c) for c in l] for l in inmatr1]
        matr2 = [[decimal.Decimal(c) for c in l] for l in inmatr2]
        for i in range(dimv):
            result.append([])
            for j in range(dimh):
                result[i].append(decimal.Decimal(k1) * matr1[i][j] + decimal.Decimal(k2) * matr2[i][j])
    else:
        matr1 = inmatr1
        matr2 = inmatr2
        for i in range(dimv):
            result.append([])
            for j in range(dimh):
                result[i].append(k1 * matr1[i][j] + k2 * matr2[i][j])
    return result


def ludecomp(matrix):
    dim = len(matrix)
    lmatrix = [[0.0 for i in range(dim)] for j in range(dim)]
    umatrix = [[0.0 for i in range(dim)] for j in range(dim)]
    for i in range(dim):
        umatrix[0][i] = matrix[0][i]
        lmatrix[i][0] = matrix[i][0] / umatrix[0][0]
        lmatrix[i][i] = 1.0

    for i in range(1, dim):
        for j in range(i, dim):
            sum = 0.0
            for k in range(i):
                sum = sum + lmatrix[i][k] * umatrix[k][j]
            umatrix[i][j] = matrix[i][j] - sum
        for j in range(i + 1, dim):
            sum = 0.0
            for k in range(i):
                sum = sum + lmatrix[j][k] * umatrix[k][i]
            lmatrix[j][i] = (matrix[j][i] - sum) / umatrix[i][i]
    return [lmatrix, umatrix]


def qrdecomp(inmatrix, copy=0):
    if (copy == 1):
        matrix = deepcopy(inmatrix)
    else:
        matrix = inmatrix
    dim = len(matrix)
    R = [[0.0 for i in range(dim)] for j in range(dim)]
    Q = [[0.0 for i in range(dim)] for j in range(dim)]

    for i in range(dim):
        R[i][i] = euclidnorm([[matrix[k][i]] for k in range(dim)])
        for k in range(dim):
            Q[k][i] = matrix[k][i] / R[i][i]
        for j in range(i + 1, dim):
            R[i][j] = 0.0
            for k in range(dim):
                R[i][j] = R[i][j] + Q[k][i] * matrix[k][j]
            for k in range(dim):
                matrix[k][j] = matrix[k][j] - Q[k][i] * R[i][j]
    return [Q, R]


def eigenvalue(inmatrix, prec=5,initqr=None,backupfile=None,backup_freq=10): #Расчет собственных значений методом QR-разложения
    matrix = deepcopy(inmatrix)
    dim = len(matrix)
    e = 0.5 * (0.1 ** prec)
    oldl = [1.0 for i in range(dim)]
    diff = oldl[:]
    flag = 1
    if (initqr == None):
        qr = qrdecomp(matrix)
    else:
        qr=deepcopy(initqr)
    if (backupfile != None):
        eigfile=open(backupfile,'w')
    counter=0
    while (flag != 0):
        counter+=1
        matrix = multiplymatrix(qr[1], qr[0])
        diff = [matrix[i][i] - oldl[i] for i in range(dim)]
        oldl = [matrix[i][i] for i in range(dim)]
        qr = qrdecomp(matrix)
        flag = 0
        for i in range(dim):
            if (abs(diff[i]) > e):
                flag += 1
        if (counter == backup_freq) and (backupfile != None):
            counter=0
            eigfile.seek(0)
            eigfile.write(matrixtostring(qr[0],prec=16,field=20)+matrixtostring(qr[1],prec=16,field=20))
            eigfile.flush()
    if (backupfile != None):
        eigfile.close()
    result=[round(matrix[i][i],prec) for i in range(dim)]
    result=list(set(result))
    result.sort(reverse=1)
    return [[value] for value in result]


def eigenvalue_sim(inmatrix, prec=5,initqr=None,backupfile=None,backup_freq=10): #Расчет собственных значений методом QR-разложения
    matrix = deepcopy(inmatrix)
    dim = len(matrix)
    e = 0.5 * (0.1 ** prec)
    oldl = [1.0 for i in range(dim)]
    diff = oldl[:]
    flag = 1
    if (initqr == None):
        qr = qrdecomp(matrix)
        eigvec=deepcopy(qr[0])
    else:
        qr=deepcopy(initqr[0:2])
        eigvec=deepcopy(initqr[2])
    if (backupfile != None):
        eigfile=open(backupfile,'w')
    counter=0
    while (flag != 0):
        counter+=1
        matrix = multiplymatrix(qr[1], qr[0])
        diff = [matrix[i][i] - oldl[i] for i in range(dim)]
        oldl = [matrix[i][i] for i in range(dim)]
        qr = qrdecomp(matrix)
        eigvec=multiplymatrix(eigvec,qr[0])
        flag = 0
        for i in range(dim):
            if (abs(diff[i]) > e):
                flag += 1
        if (counter == backup_freq) and (backupfile != None):
            counter=0
            eigfile.seek(0)
            eigfile.write(matrixtostring(qr[0],prec=16,field=20)+matrixtostring(qr[1],prec=16,field=20)+matrixtostring(eigvec,prec=16,field=20))
            eigfile.flush()
    if (backupfile != None):
        eigfile.close()
    result=[round(matrix[i][i],prec) for i in range(dim)]
    return (eigvec,[[value] for value in result])


def gaussmethod(inmatrix, infreeterm, copy=0, precise=0, precision=13,verbose=0):
    dim = len(inmatrix)
    rank=dim
    if (precise == 1):
        decimal.getcontext().prec=precision
        epsilon=decimal.Decimal(0.0)
        matrix = [[decimal.Decimal(inmatrix[i][j])/decimal.Decimal(1.0) for j in range(dim)] for i in range(dim)]
        freeterm = [[decimal.Decimal(infreeterm[i][0])] for i in range(dim)]
    else:
        epsilon=10.0**(-precision) #Граница точности типа float. Если число меньше epsilon, то можно считать, что это нуль
        if (copy == 1):
            matrix = deepcopy(inmatrix)
            freeterm = deepcopy(infreeterm)
        else:
            matrix = inmatrix
            freeterm = infreeterm

    #Прямой ход
    for k in range(dim):
        #Поиск лидирующего элемента
        s = [abs(matrix[i][k]) for i in range(k, dim)]
        indexmax = s.index(max(s)) + k
        #Проверка на наличие нулевых строк - т.е. линейно зависимых
        if (abs(matrix[indexmax][k]) <= epsilon):#Если нашлись линейно зависимые строки
            rank=k
            if (verbose==1): print("Обнаружены линейно зависимые строки! Ранг матрицы {0:d}".format(rank))
            break#Прерываем прямой ход
        #Перестановка строк местами
        for i in range(dim):
            buff = matrix[indexmax][i]
            matrix[indexmax][i] = matrix[k][i]
            matrix[k][i] = buff
        buff = freeterm[indexmax][0]
        freeterm[indexmax][0] = freeterm[k][0]
        freeterm[k][0] = buff
        #
        leader = matrix[k][k]
        for i in range(k, dim):
            matrix[k][i] = matrix[k][i] / leader
        freeterm[k][0] = freeterm[k][0] / leader
        for i in range(k + 1, dim):
            leader = matrix[i][k]
            freeterm[i][0] = freeterm[i][0] - freeterm[k][0] * leader
            for j in range(k, dim):
                matrix[i][j] = matrix[i][j] - matrix[k][j] * leader
        #print "k=", k
        #print matrixtostring(matrix, field=6)
    #Подготовка векторов решений
    if (rank == dim):
        xvec = [[[0.0] for i in range(dim)]]
        solvq=1 #Количество решений
    elif (rank < dim): #Если ранг матрицы меньше размерности, то подготовить dim-rank векторов решений
        solvq=dim-rank
        xvec = [[[0.0] for i in range(dim)] for k in range(solvq)]
        for k in range(solvq):
            for i in range(dim-1,rank-1,-1):
                if (precise == 1):
                    if (i == k+rank):
                        xvec[k][i][0]=decimal.Decimal(1.0)
                    else:
                        xvec[k][i][0]=decimal.Decimal(0.0)
                else:
                    if (i == k+rank):
                        xvec[k][i][0]=1.0

    #Для каждого решения
    for k in range(solvq):
        #Обратный ход
        for i in range(rank-1, -1, -1):
            xvec[k][i][0] = freeterm[i][0]
            for j in range(i + 1, dim):
                xvec[k][i][0] = xvec[k][i][0] - matrix[i][j] * xvec[k][j][0]
    return xvec


def determinant(inmatrix, copy=0, precise=0,size=0):
    if (size == 0 ):
        dim = len(inmatrix)
    else:
        dim=size
	matrix=[[inmatrix[i][j] for j in range(dim)] for i in range(dim)]
    if (precise == 1):
        matrix = [[decimal.Decimal(inmatrix[i][j]) for j in range(dim)] for i in range(dim)]
    else:
        if (copy == 1):
            matrix = deepcopy(inmatrix)
        else:
	    if (size == 0):
            	matrix = inmatrix
    permutation=0.0
    #Прямой ход
    for k in range(dim):
        #Поиск лидирующего элемента в k-том столбце
        s = [abs(matrix[i][k]) for i in range(k, dim)]
        indexmax = s.index(max(s)) + k #Определение индекса строки, содержащей лидирующий элемент
        if (indexmax != k): #Если это не первая строка в анализируемом диапазоне
            for i in range(dim): #Перестановка первой строки и строки с лидирующим элементом местами
                buff = matrix[indexmax][i]
                matrix[indexmax][i] = matrix[k][i]
                matrix[k][i] = buff
            permutation+=1 #Увеличение на единицу счетчика перестановок
	#���� ������� ������� ������� (��� ������� ������� � ������� ��� ���� ������� ��������� ������),
        #�� ��������� �������� � ���������� ���� - ������������ � ������ ���� �������
        if (s[indexmax-k] == 0.0):
            return 0.0
        ########
        for i in range(k + 1, dim):
            leader = matrix[i][k]
            for j in range(k, dim):
                matrix[i][j] = matrix[i][j] - (matrix[k][j] / matrix[k][k]) * leader
#        print "k=", k
#        print matrixtostring(matrix, field=6)
    #Вычисление определителя
    result=1.0
    for k in range(dim):
        result=result*matrix[k][k]
    result=result*((-1.0)**permutation)
    return result


def matrixtostring(matrix, prec=2, field=8,borders=""):
    string = ""
    dim = len(matrix)
    for i in range(dim):
        string=string+borders
        for koef in matrix[i]:
            string = string + " {0:{1}.{2}f}".format(koef, field, prec)
        string = string + " "+borders+"\n"
    return string


def frobeniusnorm(matrix):
    dim = len(matrix)
    result = 0.0
    for i in range(dim):
        for j in range(dim):
            result = result + abs(matrix[i][j]) ** 2
    return sqrt(result)


def euclidnorm(vector, precise=0):
    dim = len(vector)
    if (precise == 1):
        result = decimal.Decimal(0.0)
    else:
        result = 0.0
    for i in range(dim):
        result = result + abs(vector[i][0]) ** 2
    return sqrt(result)


def scalarmult(vec1, vec2):
    dim = len(vec1)
    result = 0.0
    for i in range(dim):
        result = result + vec1[i][0] * vec2[i][0]
    return result

def transpose(matrix):
    vdim=len(matrix)
    hdim=len(matrix[0])
    return [[matrix[j][i] for j in range(vdim)] for i in range(hdim)]