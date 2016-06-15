#! /usr/bin/python -u
#-*- coding: utf-8 -*-
__author__ = 'MGI'

import sys
import math
from matrixops import plane_3points, matrixtostring, scalarmult, euclidnorm, normalize, gaussmethod
from copy import deepcopy
from xdrfile import *


helpline = ''

plane_atoms = \
    {"U"  : ('C5','N1','N3'),\
     "C"  : ('C5','N1','N3'),\
     "DT" : ('C5','N1','N3'),\
     "DC" : ('C5','N1','N3'),\
     "Cm" : ('C5','N1','N3'),\
     "m5C": ('C5','N1','N3'),\
     "Um" : ('C5','N1','N3'),\
     "m5U": ('C5','N1','N3'),\
     "D"  : ('C5','N1','N3'),\
     "Y"  : ('C5','N1','N3'),\
     "m3Y": ('C5','N1','N3'),\
     "A"  : ('C6','N3','C8'),\
     "DA" : ('C6','N3','C8'),\
     "G"  : ('C6','N3','C8'),\
     "DG" : ('C6','N3','C8'),\
     "Gm" : ('C6','N3','C8'),\
     "m1G": ('C6','N3','C8'),\
     "m2G": ('C6','N3','C8'),\
     "m7G": ('C6','N3','C8'),\
     "m2A": ('C6','N3','C8'),\
     "m6A": ('C6','N3','C8'),\
     "28A": ('C6','N3','C8')}

def parse_args(): #Разбор аргументов командной строки
    arguments = [["-top", "data.pdb"],\
                 ["-trj", "data.xtc"],\
                 ["-angle", 30.0],\
                 ['-out', 'angles.dat'],\
                 ['-residues',None],\
                 ['-pairs',None],\
                 ['-mindist',0.28],\
                 ['-maxdist',0.40],
                 ['-timestep',1.0],
                 ['-pair_maxdist',0.1],
                 ['-pair_cutoff',1.0],
                 ['-stack_cutoff',0.5]]
    arguments_number = len(sys.argv)
    for i in xrange(0, arguments_number):
        try:
            if (sys.argv[i] == "-top"):
                arguments[0][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-trj"):
                arguments[1][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-angle"):
                arguments[2][1] = float(sys.argv[i + 1])
            elif (sys.argv[i] == "-out"):
                arguments[3][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-residues"):
                arguments[4][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-pairs"):
                arguments[5][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-mindist"):
                arguments[6][1] = float(sys.argv[i + 1])/10.0
            elif (sys.argv[i] == "-maxdist"):
                arguments[7][1] = float(sys.argv[i + 1])/10.0
            elif (sys.argv[i] == "-timestep"):
                arguments[8][1] = float(sys.argv[i + 1])
            elif (sys.argv[i] == "-pair_maxdist"):
                arguments[9][1] = float(sys.argv[i + 1])
            elif (sys.argv[i] == "-pair_cutoff"):
                arguments[10][1] = float(sys.argv[i + 1])
            elif (sys.argv[i] == "-stack_cutoff"):
                arguments[11][1] = float(sys.argv[i + 1])
            elif (sys.argv[i] == "-h"):
                print helpline
        except IndexError:
            print("Ошибка! Не хватает аргументов командной строки для корректного выполнения программы.")
    return arguments

def angle_planes(normal1,normal2):
    return math.acos(scalarmult(normal1,normal2)/(euclidnorm(normal1)*euclidnorm(normal2)))

class Index():
    def __init__(self, atoms):
        self.__tree__ = dict()
        i=0
        for atom in atoms:
            if atom[3] not in self.__tree__.keys():
                self.__tree__.update({atom[3]: {}})
            if atom[5] not in self.__tree__[atom[3]].keys():
                self.__tree__[atom[3]].update({atom[5]: {}})
            self.__tree__[atom[3]][atom[5]].update({atom[1]: i})
            i+=1
    def __getitem__(self, item):
        try:
            return self.__tree__[item]
        except KeyError:
            return None
    def newindex(self,item):
        self.__tree__.update(item)

    def __contains__(self, item):
        if item in self.__tree__:
            return True
        else:
            return False


class PDB():
    def __init__(self, file=''):
        self.atoms = []
        self.bonds = []
        self.atom_counter = 0
        self.HB_counter = 0

        if file != '':
            pdbin_file = open(file, 'r')
            pdbin = pdbin_file.readlines()
            pdbin_file.close()
            #Выделение кристаллографических данных
            cryst = ''
            for line in pdbin:
                if line[:6] == 'CRYST1':
                    cryst = line
                elif line[:6] == 'ATOM  ' or line[:6] == 'HETATM':
                    self.atoms.append(line)
                elif line[:6] == 'CONECT':
                    self.bonds.append(line)
                #Разбиение на поля
            self.atoms = map( \
                lambda line: [int(line[6:11].strip()), \
                              line[12:16].strip(), \
                              line[16], \
                              line[17:20].strip(), \
                              line[21], \
                              int(line[22:26].strip()), \
                              line[26], \
                              float(line[30:38].strip()), \
                              float(line[38:46].strip()), \
                              float(line[46:54].strip()), \
                              float(line[54:60].strip()), \
                              float(line[60:66].strip()), \
                              line[76:78].strip(), \
                              line[78:80].strip()],
                #Номер атома
                #Имя атома
                #Альтернативное положение
                #Остаток
                #Цепь
                #Номер остатка
                #Код вставки остатка
                #X
                #Y
                #Z
                #Занятость
                #B-фактор
                #Элемент
                #Заряд
                self.atoms)

            self.bonds = map( \
                lambda line: map( \
                    lambda x: int(x), \
                    filter( \
                        lambda x: x != '', \
                        line.split())[1:]), \
                self.bonds)
            #Преобразование связей в плоский список пар нумеров атомов
            bonds_back = list( \
                set( \
                    reduce( \
                        lambda x, y: x + y, \
                        map( \
                            lambda x: map( \
                                lambda y: tuple(sorted([x[0], y])), \
                                x[1:]), \
                            self.bonds))))
            self.bonds = map(lambda x: set(x), bonds_back)
            #Индексация остатков
            self.index = Index(self.atoms)
        else:
            self.index=Index([])

    def reindex(self):
        self.index = Index(self.atoms)

    def add_atom(self, name, resname, resnumber, chain, elem, x, y, z, altloc='', insertion='', occupancy=1.0,
                 bfact=0.0, charge=0.0):
        if len(self.atoms)==0:
            lastindex=0
        else:
            lastindex=self.atoms[-1][0]
        #Добавление атома
        self.atoms.append([lastindex+1,name,altloc,resname,chain,resnumber,insertion,x,y,z,occupancy,bfact,elem,charge])
        #Перестройка индекса
        if resname not in self.index:
            self.index.newindex({resname:{}})
        if resnumber not in self.index[resname]:
            self.index[resname].update({resnumber:{}})
        self.index[resname][resnumber].update({name:len(self.atoms)-1})
        return lastindex+1

    def add_bond(self,number1,number2):
        self.bonds.append([number1,number2])

    def droptofile(self, file):
        text = ''
        for atom in self.atoms:
            text += "{0:<6s}{1:>5d}{2:>5s}{3:>1s}{4:>3s}{5:>2s}{6:>4d}{7:>1s}{8:>11.3f}{9:>8.3f}{10:>8.3f}{11:>6.2f}{12:>6.2f}{13:>12s}{14:>2.0f}\n".format(
                'ATOM', \
                atom[0],  atom[1],  atom[2],  atom[3], atom[4], \
                atom[5],  atom[6],  atom[7],  atom[8], atom[9], \
                atom[10], atom[11], atom[12], atom[13])
        for bond in self.bonds:
            text += "{0:<6s}{1:>5d}{2:>5d}\n".format('CONECT', bond[0], bond[1])
        pdbout_file = open(file, 'w')
        pdbout_file.write(text)
        pdbout_file.close()
    pass


class XPM():
    def __init__(self,table,discrete=True):

        self.header='static char *gromacs_xpm[] = {'
        self.tail='};'
        try:
            self.height=len(table[0])
            self.width=len(table)
        except IndexError:
            print("Ошибка! Таблица не содержит данных.")
            raise SystemExit(1)
        self.palette=dict()
        self.data=deepcopy(table)
        self.colors=dict()
        self.xtics=[]
        self.ytics=[]
        self.xlabel=''
        self.ylabel=''
        self.title=''
        self.legend=''
        self.isdiscrete=discrete
        self.codechars=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', \
                        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',\
                        'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', \
                         'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',\
                         'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',\
                         '!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', \
                         ':', ';', '<', '=', '>', '?', '@',\
                        '[', '\\', ']', '^', '_', '`', '{', '|', '}', '~']

    def set_colors(self,colors,char_per_color=1,nsteps=0):
        def encode_by_char(index,field):
            code=''
            codepage_length=len(self.codechars)
            while index>=codepage_length:
                code=self.codechars[index % codepage_length]+code
                index=index / codepage_length
            else:
                code=self.codechars[index]+code
            while len(code)<field:
                code=self.codechars[0]+code
            return code

        #Пример colors={1.0:('N','#ff0000')}
        #Дискретные данные
        if self.isdiscrete==True:
            self.colors=colors
            self.char_per_color=char_per_color
        #Если данные не дискретного характера, то преобразовываем их согласно подгруженной сетке
        else:
            #Сколько символов нужно на цвет?
            self.char_per_color=int(\
                math.ceil(\
                    math.log(\
                        nsteps,\
                        len(self.codechars))))
            #Построение сетки
            minvalue,maxvalue=tuple(colors.keys())
            self.minvalue=minvalue
            mincolor=colors[minvalue][1]
            maxcolor=colors[maxvalue][1]
            mincolor=[int('0x'+mincolor[1:3],16),int('0x'+mincolor[3:5],16),int('0x'+mincolor[5:7],16)]
            maxcolor=[int('0x'+maxcolor[1:3],16),int('0x'+maxcolor[3:5],16),int('0x'+maxcolor[5:7],16)]
            delta_color=[maxcolor[0]-mincolor[0],maxcolor[1]-mincolor[1],maxcolor[2]-mincolor[2]]
            delta_value=maxvalue-minvalue
            self.step=delta_value/nsteps
            for i in xrange(nsteps+1):
                R=int(mincolor[0]+(i*1.0/nsteps)*delta_color[0])
                G=int(mincolor[1]+(i*1.0/nsteps)*delta_color[1])
                B=int(mincolor[2]+(i*1.0/nsteps)*delta_color[2])
                value=round(minvalue+(i*1.0/nsteps)*delta_value,8)
                self.colors.update({value:(encode_by_char(i,self.char_per_color),'#{0:02x}{1:02x}{2:02x}'.format(R,G,B),value)})

    def droptofile(self,filename):
        #Параметры картинки
        values='"{0:>8d} {1:>8d} {2:>5d} {3:>2d}",\n'.format(self.width,self.height,len(self.colors),self.char_per_color)
        #Цвета
        colors=reduce(\
            lambda x,y:x+y,\
            map(\
                lambda s:'"{0:<4s} c {1:>8s} " /* "{2:s}" */,\n'.format(self.colors[s][0],self.colors[s][1],str(self.colors[s][2])),\
                sorted(\
                    self.colors)),\
            '')
        #Метки на осях
        #Ось X
        xlabels=''
        i=0
        length=len(self.xtics)
        while i<length:
            if i % 20 ==0:
                xlabels+='/* x-axis:'
            xlabels+=' '+self.xtics[i]+' '
            i+=1
            if i % 20 ==0:
                xlabels+='*/\n'
        else:
            if xlabels[-1]!='\n':
                xlabels+='*/\n'
        #Ось Y
        ylabels=''
        i=0
        length=len(self.ytics)
        while i<length:
            if i % 20 ==0:
                ylabels+='/* y-axis:'
            ylabels+=' '+self.ytics[i]+' '
            i+=1
            if i % 20 ==0:
                ylabels+='*/\n'
        else:
            if ylabels[-1]!='\n':
                ylabels+='*/\n'
        #Таблица
        table=''
        string=''
        minvalue=self.colors.keys()[0]
        for i in xrange(self.height-1,-1,-1):
            string='"'
            for j in xrange(self.width):
                if self.isdiscrete==True:
                    string+=self.colors[self.data[j][i]][0]
                else:
                    k=math.floor((self.data[j][i]-self.minvalue)/self.step)
                    index=round(self.minvalue+k*self.step,8)
                    string+=self.colors[index][0]
            string+='",\n'
            table+=string
        table=table[:-2]+'\n'
        #
        if self.isdiscrete==False:
            type='Continuous'
        else:
            type='Discrete'
        #Вывод в файл
        xpmout=open(filename,'w')
        #Заголовок
        xpmout.write('/* XPM */\n/* Generated by g_rms */\n/* This file can be converted to EPS by the GROMACS program xpm2ps */\n')
        xpmout.write('/* title: "{0:s}" */\n/* legend: "{1:s}" */\n/* x-label: "{2:s}" */\n/* y-label: "{3:s}" */\n/* type:      "{4:s}" */\n'.format(self.title,self.legend,self.xlabel,self.ylabel,type))
        xpmout.write('static char *gromacs_xpm[] = {\n')
        #Параметры картинки и цвета
        xpmout.write('{0:s}{1:s}'.format(values,colors))
        #Подписи осей
        xpmout.write('{0:s}{1:s}'.format(xlabels,ylabels))
        #Данные
        xpmout.write(table)
        xpmout.close()

    pass







#Аргументы командной строки
arguments = parse_args()
#Чтение и индексация PDB
topology = PDB(file=arguments[0][1])
#Чтение выбранных остатков или пар
if (arguments[4][1]!=None) and (arguments[5][1]==None):
    #Чтение остатков
    residues_file=open(arguments[4][1],'r')
    residues=residues_file.readlines()
    residues_file.close()
    residues=map(\
        lambda res: (res[0].strip(),int(res[1])),\
        map(\
            lambda res: res.split(),\
            residues))
    #Изготовление пар
    length=len(residues)
    pairs=reduce(\
        lambda x,y: x+y,\
        [[(i,j) for j in range(i+1,length)] for i in range(length)],\
        [])
elif (arguments[4][1]==None) and (arguments[5][1]!=None):
    #Чтение готовых пар
    residues_file=open(arguments[5][1],'r')
    residues=residues_file.readlines()
    residues_file.close()
    pairs=map(\
        lambda res: [(res[0].strip(),int(res[1])),(res[2].strip(),int(res[3]))],\
        map(\
            lambda res: res.split(),\
            residues))
    print pairs
    #Изготовление списка остатков
    residues=list(\
        set(\
            reduce(\
                lambda x,y:x+y,\
                pairs,\
                [])))
    print residues
    #Изготовление списка пар, соответствующих списку остатков
    pairs=map(\
        lambda pair:(residues.index(pair[0]),residues.index(pair[1])),\
        pairs)
    print pairs
#Индексация атомов выбранных плоскостей
planes=[]
for residue in residues:
    planes.append((\
        topology.index[residue[0]][residue[1]][plane_atoms[residue[0]][0]],\
        topology.index[residue[0]][residue[1]][plane_atoms[residue[0]][1]],\
        topology.index[residue[0]][residue[1]][plane_atoms[residue[0]][2]]))


#Чтение траектории
xtc=xdrfile(arguments[1][1])
#Анализ стэкинга
stacking=[]
equations=[[] for plane in planes]
residue_quantity=len(residues)
times=[]
#По всем состояниям
count=0
for frame in xtc:
    stacking+=[[]]
    #stacking[-1]+=[frame.step]
    times+=[str(frame.time/1000)]
    #По всем остаткам считаем уравнения плоскостей
    for i in range(residue_quantity):
        try:
            atom1=frame.x[planes[i][0]]
            atom2=frame.x[planes[i][1]]
            atom3=frame.x[planes[i][2]]
            p1=[[float(atom1[0])],[float(atom1[1])],[float(atom1[2])]]
            p2=[[float(atom2[0])],[float(atom2[1])],[float(atom2[2])]]
            p3=[[float(atom3[0])],[float(atom3[1])],[float(atom3[2])]]
            equations[i]=plane_3points(p1, p2, p3)
        except ZeroDivisionError:
            print atom1,atom2,atom3
            print atom1[1],atom2[1]

    #По всем парам
    for pair in pairs:
        i=pair[0]
        j=pair[1]
        #Вычисляем межплоскостной угол
        normal_1=equations[i][:-1]
        normal_2=equations[j][:-1]
        try:
            angle=angle_planes(normal_1,normal_2)*180.0/math.pi
        except ZeroDivisionError:
            print normal_1,normal_2
        if angle > arguments[2][1]:
            stacking[-1]+=[0]#Нет стэкинга
        else:
            #Плоскости оснований
            plane_1=equations[i]
            plane_2=equations[j]
            #Вычисляем межплоскостные расстояния
            #Спорный способ, поэтому закомментируем
            #distance=0.0
            #for k in range(3):
            #    atom_1=frame.x[planes[i][k]]
            #    t=-(float(atom_1[0])+float(atom_1[1])*plane_2[1][0]+float(atom_1[2])*plane_2[2][0]+plane_2[3][0])\
            #      /(1.0+plane_1[1][0]*plane_2[1][0]+plane_1[2][0]*plane_2[2][0])
            #    distance_i=t*math.sqrt(1+plane_1[1][0]**2+plane_1[2][0]**2)
            #    distance+=distance_i
            #distance=distance/3.0
            #Вычисляем расстояние между центром одного основания и пересечением нормали из центра другого
            #Центры оснований
            base_center_1=[(float(frame.x[planes[i][0]][0])+float(frame.x[planes[i][1]][0])+float(frame.x[planes[i][2]][0]))/3,\
                           (float(frame.x[planes[i][0]][1])+float(frame.x[planes[i][1]][1])+float(frame.x[planes[i][2]][1]))/3,\
                           (float(frame.x[planes[i][0]][2])+float(frame.x[planes[i][1]][2])+float(frame.x[planes[i][2]][2]))/3]
            base_center_2=[(float(frame.x[planes[j][0]][0])+float(frame.x[planes[j][1]][0])+float(frame.x[planes[j][2]][0]))/3,\
                           (float(frame.x[planes[j][0]][1])+float(frame.x[planes[j][1]][1])+float(frame.x[planes[j][2]][1]))/3,\
                           (float(frame.x[planes[j][0]][2])+float(frame.x[planes[j][1]][2])+float(frame.x[planes[j][2]][2]))/3]
            #Точка проекции
            t=-(base_center_1[0]+base_center_1[1]*plane_2[1][0]+base_center_1[2]*plane_2[2][0]+plane_2[3][0])\
                  /(1.0+plane_1[1][0]*plane_2[1][0]+plane_1[2][0]*plane_2[2][0])
            projection=[base_center_1[0]+plane_1[0][0]*t,\
                        base_center_1[1]+plane_1[1][0]*t,\
                        base_center_1[2]+plane_1[2][0]*t]
            #Искомое расстояние
            base_centers_distance=math.sqrt((base_center_2[0]-projection[0])**2+(base_center_2[1]-projection[1])**2+(base_center_2[2]-projection[2])**2)
            #Межплоскостное расстояние по нормали, исходящей из центра первого основания
            distance=t*math.sqrt(1+plane_1[1][0]**2+plane_1[2][0]**2)

            if base_centers_distance > arguments[10][1]:
                stacking[-1]+=[0] #Стэкинга нет
            elif arguments[11][1] <= base_centers_distance <= arguments[10][1]:
                if distance <= arguments[8][1]:
                    stacking[-1]+=[-1] #Стэкинга нет, но основания, скорее всего, образовали пару
                else:
                    stacking[-1]+=[0] #Стэкинга нет
            else:
                if (arguments[6][1] <= distance <= arguments[7][1]):
                    stacking[-1]+=[1] #Стэкинг есть
                else:
                    stacking[-1]+=[0] #Стэкинга нет

    count+=1
    if count % 100 == 0: print 'Кадр {0:s} нс'.format(times[-1])

#Вывод результата
#Вывод в виде простой таблицы. Работает. Закомментирован в пользу вывода в формат *.xpm и формат таблицы встречаемости
#for s in stacking:
#    s[0]=str(s[0])
# title=reduce(\
#     lambda x,y: x+y,\
#     map(\
#         lambda s: "{0:>15s}".format(s),\
#         map(\
#             lambda pair: residues[pair[0]][0]+'-'+str(residues[pair[0]][1])+'|'+residues[pair[1]][0]+'-'+str(residues[pair[1]][1]),\
#             pairs)),\
#     '#Step          ')
# output=reduce(\
#     lambda x,y:x+y+'\n',\
#     map(\
#         lambda s:reduce(\
#             lambda i,j:i+j,\
#             map(\
#                 lambda p: '{0:>15s}'.format(p),\
#                 s),\
#             ''),\
#         stacking),\
#     '')
#outfile=open(arguments[3][1],'w')
#outfile.write(title+'\n')
#outfile.write(output)
#Вывод таблицы встречаемости
#Таблица для сбора встречаемости
occurency=map(\
    lambda pair: [residues[pair[0]][0]+'-'+str(residues[pair[0]][1]),\
                  residues[pair[1]][0]+'-'+str(residues[pair[1]][1]),\
                  pair[0],\
                  pair[1],\
                  0,0,0],\
    pairs)
#xpm_pairs=map(lambda pair: pair[0]+'/'+pair[1])
#По всем состояниям считаем события
length=len(pairs)
for frame in stacking:
    for i in xrange(length):
        if frame[i]==1:
            occurency[i][4]+=1
        elif frame[i]==-1:
            occurency[i][5]+=1
        elif frame[i]==0:
            occurency[i][6]+=1

frame_quant=len(stacking)
occurency_out=reduce(\
    lambda x,y: x+y,\
    map(\
        lambda s: '{0:<10s} {1:<10s} {2:>09d} {3:>09d} {4:9.4f} {5:9.4f} {6:9.4f}\n'.format(s[0],s[1],s[2],s[3],s[4],s[5],s[6]),\
        map(\
            lambda pair: pair[:4]+[(100.0*pair[4])/frame_quant,(100.0*pair[5])/frame_quant,(100.0*pair[6])/frame_quant],\
            occurency)),\
    '{0:<10s} {1:<10s} {2:>9s} {3:>9s} {4:9s} {5:9s} {6:9s} Состояний: {7:d}\n'.format('#Остаток 1','Остаток 2','i','j','Стэкинг,%','Пара,%','Нет,%',frame_quant))
outfile=open(arguments[3][1]+'.dat','w')
outfile.write(occurency_out+'\n')
occurency_out=reduce(\
    lambda x,y: x+y,\
    map(\
        lambda s: '{0:<10s} {1:<10s} {2:>09d} {3:>09d} {4:9.4f} {5:9.4f} {6:9.4f}\n'.format(s[0],s[1],s[2],s[3],s[4],s[5],s[6]),\
        map(\
            lambda pair: [pair[1],pair[0],pair[3],pair[2],(100.0*pair[4])/frame_quant,(100.0*pair[5])/frame_quant,(100.0*pair[6])/frame_quant],\
            occurency)),\
    '')
outfile.write(occurency_out)
occurency_out=reduce(\
    lambda x,y: x+y,\
    map(\
        lambda s: '{0:<10s} {1:<10s} {2:>09d} {3:>09d} {4:9.4f} {5:9.4f} {6:9.4f}\n'.format(s[0],s[1],s[2],s[3],s[4],s[5],s[6]),\
        [[residues[i][0]+'-'+str(residues[i][1]),residues[i][0]+'-'+str(residues[i][1]),i,i,0.0,0.0,0.0] for i in xrange(len(residues))]),\
    '')
outfile.write(occurency_out)
outfile.close()

#Вывод *.xpm
xpm=XPM(stacking)
xpm.xlabel='Time, ns'
xpm.ylabel='Residues'
xpm.title='Stacking interactions'
xpm.legend=''
xpm.xtics=times
#xpm2ps не умеет принимать текстовые подписи
#xpm.ytics=xpm_pairs
xpm.ytics=[str(i+1) for i in xrange(len(pairs))]
xpm.set_colors({0:('N','#ffffff','None'),1:('S','#ff0000','Stacking'),-1:('P','#00ff00','Pair')})
xpm.droptofile(arguments[3][1]+'.xpm')

print 'Готово'
