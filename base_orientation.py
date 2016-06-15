#! /usr/bin/python -u
#-*- coding: utf-8 -*-
__author__ = 'MGI'

import sys
import math
from matrixops import plane_3points, matrixtostring, scalarmult, euclidnorm, normalize, gaussmethod
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
    arguments = [["-top", "data.pdb"], ["-trj", "data.xtc"], ["-base", ['', 0]], ['-out', 'angles.dat'],['-residues','residues.dat']]
    arguments_number = len(sys.argv)
    for i in xrange(0, arguments_number):
        try:
            if (sys.argv[i] == "-top"):
                arguments[0][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-trj"):
                arguments[1][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-base"):
                arguments[2][1][0] = sys.argv[i + 1]
                arguments[2][1][1] = int(sys.argv[i + 2])
            elif (sys.argv[i] == "-out"):
                arguments[3][1] = sys.argv[i + 1]
            elif (sys.argv[i] == "-residues"):
                arguments[4][1] = sys.argv[i + 1]
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
            self.height=len(table)
            self.width=len(table[0])
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
                self.colors.update({value:(encode_by_char(i,self.char_per_color),'#{0:02x}{1:02x}{2:02x}'.format(R,G,B))})

    def droptofile(self,filename):
        #Параметры картинки
        values='"{0:>8d} {1:>8d} {2:>5d} {3:>2d}",\n'.format(self.width,self.height,len(self.colors),self.char_per_color)
        #Цвета
        colors=reduce(\
            lambda x,y:x+y,\
            map(\
                lambda s:'"{0:<4s} c {1:>8s} " /* "{2:s}" */,\n'.format(self.colors[s][0],self.colors[s][1],str(s)),\
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
base_residue = arguments[2][1]
base_atoms = plane_atoms[base_residue[0]]
#Чтение и индексация PDB
topology = PDB(file=arguments[0][1])
#Чтение выбранных остатков
residues_file=open(arguments[4][1],'r')
residues=residues_file.readlines()
residues_file.close()
residues=map(\
    lambda res: (res[0].strip(),int(res[1])),\
    map(\
        lambda res: res.split(),\
        residues))
#Индексация атомов выбранных плоскостей
planes=[]
for residue in residues:
    planes.append((\
        topology.index[residue[0]][residue[1]][plane_atoms[residue[0]][0]],\
        topology.index[residue[0]][residue[1]][plane_atoms[residue[0]][1]],\
        topology.index[residue[0]][residue[1]][plane_atoms[residue[0]][2]]))


#Задание базисных плоскостей

#Плоскость выбранного основания
b_p1 = [[topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[0]]][7]], \
        [topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[0]]][8]], \
        [topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[0]]][9]]]
b_p2 = [[topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][7]], \
        [topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][8]], \
        [topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][9]]]
b_p3 = [[topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[2]]][7]], \
        [topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[2]]][8]], \
        [topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[2]]][9]]]
base_plane_1 = plane_3points(b_p1, b_p2, b_p3)

#Построение перпендикулярных ей и друг другу плоскостей
base_plane_1_normal = normalize(base_plane_1[:-1])
base_plane_2_normal = normalize([[topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][7] -
                                  topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[0]]][7]], \
                                 [topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][8] -
                                  topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[0]]][8]], \
                                 [topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][9] -
                                  topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[0]]][9]]])
matrix = [[base_plane_1_normal[1][0], base_plane_1_normal[2][0]], \
          [base_plane_2_normal[1][0], base_plane_2_normal[2][0]]]
term = [[-base_plane_1_normal[0][0]], \
        [-base_plane_2_normal[0][0]]]
base_plane_3_normal = normalize([[1.0]] + gaussmethod(matrix, term))


out=PDB()
ind1=out.add_atom('CEN','Bas',1,'','C',\
                   topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][7],\
                   topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][8],\
                   topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][9])
ind2=out.add_atom('1Ort','Bas',1,'','C',\
                   topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][7]+base_plane_1_normal[0][0],\
                   topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][8]+base_plane_1_normal[1][0],\
                   topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][9]+base_plane_1_normal[2][0])
ind3=out.add_atom('2Ort','Bas',1,'','C',\
                   topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][7]+base_plane_2_normal[0][0],\
                   topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][8]+base_plane_2_normal[1][0],\
                   topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][9]+base_plane_2_normal[2][0])
ind4=out.add_atom('3Ort','Bas',1,'','C',\
                   topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][7]+base_plane_3_normal[0][0],\
                   topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][8]+base_plane_3_normal[1][0],\
                   topology.atoms[topology.index[base_residue[0]][base_residue[1]][base_atoms[1]]][9]+base_plane_3_normal[2][0])

out.add_bond(ind1,ind2)
out.add_bond(ind1,ind3)
out.add_bond(ind1,ind4)
out.droptofile(arguments[0][1][:-4]+'.basis.pdb')

#Чтение траектории
xtc=xdrfile(arguments[1][1])
print sys.getsizeof(xtc)
print sys.getsizeof(topology)
#Вычисление углов
angles=[]
#По всем состояниям
for frame in xtc:
    angles+=[[]]
    angles[-1]+=[frame.step]
    #По всем углам
    for plane in planes:
	try:
            atom1=frame.x[plane[0]]
            atom2=frame.x[plane[1]]
            atom3=frame.x[plane[2]]
            p1=[[float(atom1[0])],[float(atom1[1])],[float(atom1[2])]]
            p2=[[float(atom2[0])],[float(atom2[1])],[float(atom2[2])]]
            p3=[[float(atom3[0])],[float(atom3[1])],[float(atom3[2])]]
            #c=((atom1[0]-atom3[0])*(atom1[1]-atom2[1])+(atom1[0]-atom2[0])*(atom1[1]-atom3[1]))/\
            #  ((atom1[2]-atom2[2])*(atom1[1]-atom3[1])-(atom1[2]-atom3[2])*(atom1[1]-atom2[1]))
            #b=(c*(atom1[2]-atom2[2])-atom1[0]+atom2[0])/(atom1[1]-atom2[1])
            #a=1.0
            normal=plane_3points(p1, p2, p3)[:-1]
            angle_1=angle_planes(normal,base_plane_1_normal)*180.0/math.pi
            angle_2=angle_planes(normal,base_plane_2_normal)*180.0/math.pi
            angle_3=angle_planes(normal,base_plane_3_normal)*180.0/math.pi
            angles[-1]+=[angle_1,angle_2,angle_3]
        except ZeroDivisionError:
            print atom1,atom2,atom3
            print atom1[1],atom2[1]

#Вывод результата
outfile=open(arguments[3][1],'w')
title=reduce(\
    lambda x,y: x+y,\
    map(\
        lambda s: "{0:>15s}".format(s),\
        map(\
            lambda s: s[0]+'-'+str(s[1])+'_a1',\
            residues)),\
    '#Step          ')

outfile.write(title+'\n')
outfile.write(matrixtostring(angles,prec=4,field=15))
outfile.close()
