#! /usr/bin/python -u
#-*- coding: utf-8 -*-
#usage test.py -i input.xvg -n input.ndx -o output.ndx -dih anglename -out {plot,sort}
import sys
from math import floor, ceil, trunc, log
from copy import copy, deepcopy

#Строка справки
helpline = "\
Использование: -i dihedrals.xvg -n dihedrals.ndx -o dihedrals\n\
Опции:\n\
    -diststep <значение> - шаг при вычислении распределения углов, по умолчанию 0.1.\n\
В  файле  dihedrals.xvg первый столбец  -  время, второй  -  среднее значение по группе,\n\
остальные  -  величины двугранных  углов в  том же порядке,  в котором они перечислены в\n\
файле dihedrals.ndx\n\
"
#
argumentum = list()
index = []
dihedrals = list()
covariance_matrix = []
cutmatrix = []
dihedral_quantity = 0
value_quantity = 0
binid = 0
distribution_wide = 0
current_covariance = 0.0
current_pirson = 0.0
counter = 0
R=8.314459848484848


def parse_args():#Разбор аргументов командной строки
    keys = ["-i", "-n", "-o", "-diststep"]
    arguments = [["-i", "dihedrals.xvg"],\
                 ["-n", "dihedrals.ndx"],\
                 ["-o", "analysis"],\
                 ["-diststep", 1.0],\
                 ["-entropy",False],\
                 ["-dimension",False]]
    arguments_number = len(sys.argv)
    for i in xrange(0, arguments_number):
        try:
            if (sys.argv[i] == "-i"):
                arguments[0][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-n"):
                arguments[1][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-o"):
                arguments[2][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-diststep"):
                try:
                    arguments[3][1] = float(sys.argv[i + 1])
                except ValueError:
                    print("Некорректное значение при аргументе -diststep. Принято значение по умолчанию")
            elif(sys.argv[i] == "-entropy"):
                arguments[4][1]=True
            elif(sys.argv[i] == "-dimension"):
                arguments[5][1]=True
            elif (sys.argv[i] == "-h"):
                print(helpline)
                raise SystemExit(0)
        except IndexError:
            print("Ошибка! Не хватает аргументов командной строки для корректного выполнения программы.")
    return arguments


def init_ndx(filename): #Здесь генерируется список записей о двугранных углах
    ndxin = list()
    dihedral = []
    indices = list()
    ndxlength = 0
    try:
        inf_ndx = open(filename, 'r')
        ndxin = inf_ndx.readlines()
        inf_ndx.close()
    except IOError:
        print("Ошибка чтения из файла {0:s}".format(filename))
        raise SystemExit(1)
    ndxin.pop(0)
    ndxlength = len(ndxin)
    #Разбор строк на записи
    for i in xrange(0, ndxlength):
        dihedral = ndxin[i].split(' ')
        dihedral = filter(lambda val: val != '', dihedral)
        indices.append([int(dihedral[6].strip()[2:]), \
                        dihedral[5].strip(), \
                        [dihedral[0].strip(), dihedral[1].strip(), dihedral[2].strip(), dihedral[3].strip()], \
                        0.0, 0.0, []][:])
    return indices #Возврат списка с индексами


def get_xvg(filename, dih_quant):#Извлечение значений времен и соответствующих им значений двугранных углов
    result = ["", "", "", []]
    xvg = []
    xvg_property = []
    xvg_value = []
    xvg_length = 0
    try:
        inf_xvg = open(filename, 'r')
        xvg = inf_xvg.readlines()
        inf_xvg.close()
    except IOError:
        print("Ошибка при чтении из файла {0:s}".format(filename))
    result[3] = [[] for i in xrange(0, dih_quant)]
    #Разбор *.xvg
    xvg_length = len(xvg)
    for i in xrange(0, xvg_length):
        if (xvg[i][0] != '#'): #Если не комментарий
            if (xvg[i][0] == '@'): #Если данные о *.xvg, описываемых величинах, размерностях
                xvg_property = xvg[i].split()
                if (xvg_property[1] == 'title'):#Заголовок
                    result[0] = ' '.join(xvg_property[2:])
                elif (xvg_property[1] == 'xaxis'):#Величина по абсциссе
                    result[1] = ' '.join(xvg_property[3:])
                elif (xvg_property[1] == 'yaxis'):#Величина по ординате
                    result[2] = ' '.join(xvg_property[3:])
            else:#Значения
                try:
                    xvg_value = xvg[i].split()
                    xvg_value = filter(lambda val: val != '', xvg_value)
                    #					print(len(xvg_value),xvg_value)
                    for j in xrange(1, dih_quant + 1):
                        result[3][j - 1].append(float(xvg_value[j].strip()))
                except IndexError:
                    print("Ошибка!Нарушена целостность данных в файле {0:s} строка {1:d}".format(filename, i + 1))
                    raise SystemExit(1)
    return result


def distrib_output(index, filename, diststep,min_value):
    dist_length = len(index[0][5])
    ind_length = len(index)
    sys.stdout.write("Вывод распределений...       ")
    try:
        outf = open(filename + "_distributions.dat", 'w')
        outf.write("#Величина  {0:s}\n".format(" ".join([index[k][1] for k in xrange(0, ind_length)])))
        for k in xrange(0, dist_length):
            outf.write("{0:8.3f} {1:s}\n".format(diststep * (k + 0.5) +min_value, " ".join(
                ["{0:7.4f}".format(index[i][5][k]) for i in xrange(0, ind_length)])))
        outf.close()
    except IOError:
        print("Ошибка записи в файл {0:s}".format(filename + "_distributions.dat"))
    sys.stdout.write("  Готово\n")

def entropy_output(entropy, filename):
    sys.stdout.write("Вывод энтропий...       ")
    try:
        outf = open(filename + "_entropy.dat", 'w')
        if argumentum[5][1]==True:
            outf.write("{0:20s} {1:20s}\n".format("#Value","Entropy,kJ/mol*K"))
        else:
            outf.write("{0:20s} {1:8s}\n".format("#Value","Entropy"))
        for s in entropy:
            outf.write("{0:20s} {1:8.4f}\n".format(s[0],s[1]))
        outf.close()
    except IOError:
        print("Ошибка записи в файл {0:s}".format(filename + "_entropy.dat"))
    sys.stdout.write("  Готово\n")


######################################################################
#Инициализация
argumentum = parse_args()
index = init_ndx(argumentum[1][1])
data_quantity = len(index)
data = get_xvg(argumentum[0][1], data_quantity)
print("Данные из файла {0:s} загружены".format(argumentum[0][1]))
value_quantity = len(data[3][0])
#Ищем размах величин
min_value=floor(\
    min(\
        map(\
            lambda s: min(s),\
            data[3])))
max_value=ceil(\
    max(\
        map(\
            lambda s: max(s),\
            data[3])))
#Вычисление распределения
distribution_wide = int((max_value-min_value )/ argumentum[3][1]) + 1 #Во сколько ячеек будем собирать значения
bins_basic=[0.0 for k in xrange(0, distribution_wide)] #формируем набор ячеек ширины -diststep
for i in range(distribution_wide):
    print '{0:<3d} {1:4.2f}'.format(i,min_value+argumentum[3][1]*i)
sys.stdout.write("Вычисление распределений {0:3d}%".format(0))
for i in xrange(0, data_quantity):#Для каждой величины
    index[i][5] = copy(bins_basic)
    for j in xrange(0, value_quantity):#Для каждого значения
        binid = int(
            trunc((data[3][i][j]- min_value) / argumentum[3][1])) #вычисляем номер ячейки, в которой оно лежит
        index[i][5][binid] += 1 #и увеличиваем значение в ячейке на единицу
    sys.stdout.write(" \b\b\b\b\b{0:3d}%".format(100 * (i + 1) / data_quantity))
print
#Нормируем значение в ней на общее число значений величин
index=map(\
    lambda s: s[:5]+[map(\
        lambda x: x/value_quantity,\
        s[5])],\
    index)

#Расчет энтропий распределений, если включен
if argumentum[4][1]==True:
    entropies=[]
    #По каждому распределению
    for distribution in index:
        #По всем ненулевым вероятностям
        entropy=0.0
        for probability in distribution[5]:
            if probability!=0.0:
                entropy+=probability*log(probability)
        if argumentum[5][1]==True:
            entropy=-R*entropy
        else:
            entropy=-entropy
        entropies.append((distribution[1],entropy))
    entropy_output(entropies,argumentum[2][1])

#Вывод распределений по основаниям
distrib_output(index, argumentum[2][1], argumentum[3][1],min_value)
sys.stdout.write("  Готово\n")