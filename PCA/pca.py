#! /usr/bin/python
#-*- coding: utf-8 -*-
# coding=utf-8
import sys
from matrixops import *
from pca_lib   import *
from pickle    import *

def parse_args():#Разбор аргументов командной строки
    arguments = [["-dim", 0],\
                 ["-firstcolumn", 0],\
                 ["-id", "function.dat"],\
                 ["-splitter",' '],\
                 ["-restart","none"],\
                 ["-backupfreq",10],\
                 ["-resdisp",0.0],\
                 ["-finaldim",0],\
                 ["-paramcolumn",None],
                 ["-normalize",None]]
    helpline="\
Использование: -id function.dat -dim <число переменных> -firstcolumn <первый столбец данных> ...\n\
Программа проводит анализ основных компонент подаваемого ей в файле function.dat набора данных, вычисляя матрицу\n\
ковариаций, её собственные значения и собственные вектора. Опции командной строки:\n\
    -id function.dat        - входные данные.  Допускается комментировать строки символами  '#',';','%','@','$'.\n\
    -dim <значение>         - размерность входных данных.\n\
    -firstcolumn <значение> - столбец, начиная с которого идут анализируемые данные,  по умолчанию 0.  Программа\n\
                              считывает dim столбцов начиная с firstcolumn.\n\
    -paramcolumn <значение> - номер столбца,  содержащего не  анализируемый  параметр,  например,  время.  Опция\n\
                              необязательна.Если её указать, то в обработанные данные этот столбец будет дописан\n\
                              в неизменном виде.\n\
    -splitter <значение>    - разделитель столбцов, по умолчанию пробел.\n\
    -backupfreq <значение>  - частота,  с  которой  сохраняются  промежуточные данные при вычислении собственных\n\
                              значений и собственных векторов, по умолчанию каждые 10 шагов.\n\
    -resdisp <значение>     - доля необъясненной дисперсии, проценты, по умолчанию 0.0. Чем больше это значение,\n\
                              тем меньше главных компонент оставит программа.\n\
    -finaldim <значение>    - количество главных компонент,  которые нужно оставить.  По умолчанию выводятся все\n\
                              главные компоненты.\n\
    -restart <значение>     - перезапуск программы с частично проведенным анализом. Допустимые значения:\n\
             * none         - по умолчанию. Означает выполнение анализа с самого начала, с загрузки данных и\n\
                              расчета матрицы ковариаций\n\
             * covariance   - запуск  анализа с ранее  полученной матрицей  ковариаций, сохраненной  в файле\n\
                              function.dat.covariance\n\
             * qrdecomp     - запуск анализа с незавершенным  расчетом  собственных значений.  Промежуточные\n\
                              результаты подгружаются из файла function.dat.qrdecomp\n\
             * eigenvalues  - запуск анализа с ранее полученными  собственными значениями и векторами, ранее\n\
                              сохраненными в файлах function.dat.eigenvalues и function.dat.eigenvectors\n\
             * filter       - программа  только  подгружает  исходные  данные и проецирует их в пространство\n\
                              главных компонент меньшей размерности, используя ранее полученную проецирующую\n\
                              матрицу из файла function.dat.filtermatrix. Сам анализ при этом не проводится,\n\
                              предполагается, что его провели ранее\n\
    -normalize <значение>   - файл, содержащий список границ диапазонов для нормировки входных данных.Формат\n\
                              файла таков: список строк вида   \"min max \\n \", строк столько же, сколько и\n\
                              переменных в файле данных. При нормировке min отображается в -1, max в +1.  По\n\
                              умолчанию нормировка не проводится."
    restartkeys=set(["none","covariance","qrdecomp","eigenvalues","eigenvectors","filter"])
    arguments_number = len(sys.argv)
    for i in xrange(0, arguments_number):
        try:
            if (sys.argv[i] == "-dim"):
                arguments[0][1] = int(sys.argv[i + 1])
            elif (sys.argv[i] == "-firstcolumn"):
                arguments[1][1] = int(sys.argv[i + 1])-1
            elif (sys.argv[i] == "-id"):
                arguments[2][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-splitter"):
                if (sys.argv[i+1] == "tab"):
                    arguments[3][1]='\t'
                elif (sys.argv[i+1] == "space"):
                    arguments[3][1]=' '
            elif (sys.argv[i] == "-restart"):
                if (sys.argv[i+1] in restartkeys):
                    arguments[4][1] = sys.argv[i + 1]
                else:
                    print("Ошибка! Неизвестная опция перезапуска.")
                    raise SystemExit(1)
            elif (sys.argv[i] == "-backupfreq"):
                arguments[5][1] = int(sys.argv[i + 1])
            elif (sys.argv[i] == "-resdisp"):
                arguments[6][1] = float(sys.argv[i + 1])/100
            elif (sys.argv[i] == "-finaldim"):
                arguments[7][1] = int(sys.argv[i + 1])
            elif (sys.argv[i] == "-paramcolumn"):
                arguments[8][1] = int(sys.argv[i + 1])
            elif (sys.argv[i] == "-normalize"):
                arguments[9][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-h"):
                print(helpline)
                raise SystemExit(0)
        except IndexError:
            print("Ошибка! Не хватает аргументов командной строки для корректного выполнения программы.")
            raise SystemExit(1)
        except ValueError:
            print("Ошибка! Неверный тип аргумента командной строки - {0:s}".format(sys.argv[i+1]))
            raise SystemExit(1)
    if arguments[7][1] == 0: arguments[7][1]=arguments[0][1]
    if (arguments[0][1]==1):
        print("Программа прекращает работу, ибо тщетна она, ибо подают ей на вход\n"+\
              "единственный столбец данных,  который ни от чего не зависит,  и от\n"+\
              "которого ничто не зависит")
        raise SystemExit(1)
    return arguments

def load_values(arguments):
    def empty(point):
        return point != ' ' and point != '\t' and point != ''

    comments = ['#', '%', '$', '@', ';']
    function = []
    parameters=[]
    inf1 = open(arguments[2][1], 'r')
    lines = inf1.readlines()
    inf1.close()
    for s in lines:
        if ((s[0] not in comments) and (s != '')):
            s = s.strip()
            line = s.split(' ')
            line = filter(empty, line)
            try:
                function.append([float(line[i].strip()) for i in range(arguments[1][1],arguments[1][1]+arguments[0][1])])
                if (arguments[8][1] != None): parameters.append(float(line[arguments[8][1]].strip(arguments[3][1])))
            except IndexError:
                print("Во входном файле нет нужного количества столбцов данных.\n"+\
                      "Нехорошо издеваться над программой, она и так написана на Python")
                raise SystemExit(1)
    return (parameters,function)

restart_eigenval=set(["none","covariance","qrdecomp","filter"])
restart_qrdecomp=set(["none","covariance","qrdecomp"])

argumentum=parse_args()

if ((argumentum[4][1] == "none") or (argumentum[4][1] == "filter")):
    parameters,data = load_values(argumentum)
    print("Загружены данные из файла {0:s}, начиная с {1:d} по {2:d} столбец".\
              format(argumentum[2][1],argumentum[1][1]+1,argumentum[1][1]+argumentum[0][1]))
    #Нормировка данных
    if (argumentum[9][1] != None):
        data=normalize_data(data,argumentum[9][1])

if (argumentum[4][1] == "none"):
    #Вычисление средних
    averages=calcaverages(data,argumentum[0][1])
    print("Вычислены средние значения")
    #Вычисление матрицы ковариаций
    covmatrix=covariance_matrix(data,averages,argumentum[0][1],argumentum[2][1])
    print("Рассчитана матрица ковариаций")

if (argumentum[4][1] != "none") and (argumentum[4][1] != "filter"):#Если перезапуск с уже посчитанной матрицей ковариаций
    print("Перезапуск с ранее полученной матрицы ковариаций, файл {0:s}".format(argumentum[2][1]+".covariance"))
    covmatrix=load_covmatrix(argumentum[2][1]+".covariance")

if (argumentum[4][1] in restart_qrdecomp):#Если перезапуск с незавершенного вычисления собственных чисел
    if (argumentum[4][1] == "qrdecomp"):
        print("Перезапуск с незавершенного расчета собственных значений, файл {0:s}".format(argumentum[2][1]+".qrdecomp"))
        qr=load_qrdecomp(argumentum[2][1]+".qrdecomp")
    else:
        qr=None
    #Вычисление собственных значений
    eigenvectors,eigenvalues = eigenvalue_sim(covmatrix, prec=16,initqr=qr,backupfile=argumentum[2][1]+".qrdecomp",backup_freq=argumentum[5][1])
    print("Рассчитаны собственные значения и собственные вектора")
    #Запись в файл
    eigf=open(argumentum[2][1]+".eigenvalues",'w')
    eigf.write(matrixtostring(eigenvalues,prec=16,field=20))
    eigf.close()
    eigf=open(argumentum[2][1]+".eigenvectors",'w')
    eigf.write(matrixtostring(eigenvectors,prec=16,field=20))
    eigf.close()

if (argumentum[4][1] not in restart_eigenval): #Если перезапуск с известными собственными значениями
    print("Перезапуск с ранее полученных собственных значений, файл {0:s}".format(argumentum[2][1]+".eigenvalues"))
    eigenvalues=load_eigenvalues(argumentum[2][1]+".eigenvalues")
    eigenvectors=load_eigenvectors(argumentum[2][1]+".eigenvectors")

if (argumentum[4][1] != "filter"):
    #Вычисление долей суммарной дисперсии, обуславливаемых главными компонентами
    totaldisp=0.0
    for i in range(argumentum[0][1]):
        totaldisp=totaldisp+eigenvalues[i][0]
    shares=[eigenvalues[i][0]/totaldisp for i in range(argumentum[0][1])]

    report=open(argumentum[2][1]+".report",'w')
    report.write("               + Анализ основных компонент в файле {0:s}, начиная с {1:d} по {2:d} столбец +\n".\
        format(argumentum[2][1],argumentum[1][1]+1,argumentum[1][1]+argumentum[0][1]))

    explained_disp=0.0
    i=0
    while (explained_disp <= (1.0-argumentum[6][1])) and (i<=argumentum[0][1]-1) and(i<=argumentum[7][1]-1):
        report.write("  Cобственное значение {0:.4f}\n  Доля суммарной дисперсии, вносимой компонентой {2:.2%}\n  Собственный вектор\n{1:s}\n".\
            format(eigenvalues[i][0], matrixtostring([[eigenvectors[j][i]] for j in range(argumentum[0][1])],prec=16,field=20,borders="|"),shares[i]))
        report.write("==========================================================================================\n")
        explained_disp=explained_disp+shares[i]
        i+=1
    report.close()
    #Подготовка матрицы, проецирующей вектора данных в пространство основных компонент
    counter=i
    filtermatrix=[[eigenvectors[i][j] for j in range(counter)] for i in range(argumentum[0][1])]
    filterf=open(argumentum[2][1]+".filtermatrix",'w')
    filterf.write(matrixtostring(filtermatrix,prec=16,field=20))
    filterf.close()

if (argumentum[4][1] == "filter"):
    #Загрузка прецирующей матрицы
    filterf=open(argumentum[2][1]+".filtermatrix",'r')
    filterf_lines=filterf.readlines()
    filterf.close()
    filtermatrix=readmatrix(filterf_lines)
    filterf_lines=[]
    #Проецирование
    data_filtered=[]
    value_quantity=len(data)
    for i in range(value_quantity):
        filtered_vec=multiplymatrix([data[i]],filtermatrix)
        if (argumentum[8][1] != None):
            data_filtered.append(parameters[i]+filtered_vec[0])
        else:
            data_filtered.append(filtered_vec[0])
    #Сохранение проецированных данных в файл
    fdf=open(argumentum[2][1]+".filtered",'w')
    fdf.write(matrixtostring(data_filtered,prec=4,field=10))
    fdf.close()

