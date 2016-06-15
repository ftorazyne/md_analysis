#! /usr/bin/python -u
#-*- coding: utf-8 -*-
__author__ = 'MGI'

import sys
from copy import copy

flat_groups = {
    #Нуклеотиды
    'A'   : ('N9','C8','H8','N7','C5','C6','N6','H61','H62','N1','C2','H2','N3','C4'), \
    'U'   : ('N1','C6','H6','C5','H5','C4','O4','N3','H3','C2','O2'), \
    'G'   : ('N9','C8','H8','N7','C5','C6','O6','N1','H1','C2','N2','H21','H22','N3','C4'), \
    'C'   : ('N1','C6','H6','C5','H5','C4','N4','H41','H42','N3','C2','O2'), \
    'Cm'  : ('N1','C6','H6','C5','H5','C4','N4','H41','H42','N3','C2','O2'), \
    'm5C' : ('N1','C6','H6','C5','CM','HM1','HM2','HM3','C4','N4','H41','H42','N3','C2','O2'), \
    'Gm'  : ('N9','C8','H8','N7','C5','C6','O6','N1','H1','C2','N2','H21','H22','N3','C4'), \
    'm1G' : ('CM','HM1','HM2','HM3','N1','C2','N2','H21','H22','N3','C4','C5','C6','O6','N7','C8','H8','N9'), \
    'm2G' : ('CM','HM1','HM2','HM3','N1','H1','C2','N2','H2','N3','C4','C5','C6','O6','N7','C8','H8','N9'), \
    'm7G' : ('CM','HM1','HM2','HM3','N1','H1','C2','N2','H21','H22','N3','C4','C5','C6','O6','N7','C8','H8','N9'), \
    'm2A' : ('N9','C8','H8','N7','C5','C6','N6','H61','H62','N1','C2','CM','HM1','HM2','HM3','N3','C4'), \
    'm6A' : ('N9','C8','H8','N7','C5','C6','N6','H6','CM','HM1','HM2','HM3','N1','C2','H2','N3','C4'), \
    'dmA' : ('N1','C2','H2','N3','C4','C5','C6','N6','N7','C8','H8','N9','CM1','HM1','HM2','HM3','CM2','HM4','HM5','HM6'), \
    'Um'  : ('N1','C6','H6','C5','H5','C4','O4','N3','H3','C2','O2'), \
    'm5U' : ('N1','C6','H6','C5','CM','HM1','HM2','HM3','C4','O4','N3','H3','C2','O2'), \
    'Y'   : ('N1','H1','C2','O2','N3','H3','C4','O4','C5','C6','H6'), \
    'D'   : ('N1','C2','O2','N3','H3','C4','O4','C5','H51','H52','C6','H61','H62'), \
    'm3Y' : ('N1','H1','C2','O2','N3','CM','HM1','HM2','HM3','C4','O4','C5','C6','H6'), \
    'm1A' : ('N1','CM','HM1','HM2','HM3','C2','H2','N3','C4','C5','C6','N6','H61','H62','N7','C8','H8','N9'), \
    'm3U' : ('N1','C2','O2','N3','CM','HM1','HM2','HM3','C4','O4','C5','H5','C6','H6'), \
    'mCm' : ('CM1','HM1','HM2','HM3','N1','C2','O2','N3','C4','N4','H41','C5','H5','C6','H6'), \
    '28A' : ('CM2','HM1','HM2','HM3','CM8','HM4','HM5','HM6','N1','C2','N3','C4','C5','C6','N6','H61','H62','N7','C8','N9'), \
    #Аминокислоты
    'PHE' : ('CG','CD1','HD1','CE1','HE1','CZ','HZ','CE2','HE2','CD2','HD2'), \
    'TYR' : ('CG','CD1','HD1','CE1','HE1','CZ','OH','HH','CE2','HE2','CD2','HD2'), \
    'TRP' : ('CG','CD1','HD1','NE1','HE1','CE2','CZ2','HZ2','CH2','HH2','CZ3','HZ3','CE3','HE3','CD2'), \
    'HIS' : ('CG','ND1','HD1','CE1','HE1','NE2','CD2','HD2','HE2'), \
    'ARG' : ('NE','HE','CZ','NH1','HH11','HH12','NH2','HH21','HH22'), \
    #Лиганды
    'CLM' : ( 'C6','C7','H7','C8','H8','C9','N9','O9A','O9B','C10','H10','C11','H11')}

class Arguments():
    def __init__(self):
        self.helpline='multijoin.py -files files.dat -joined data.dat -intersection \n\
                      Программа объединяет данные из файлов, перечисленных в files.dat.\n\
                      Формат записей: <ИмяФайла> <Номер ключевого поля> <Поле 1> ... <Поле N>\n'

        arguments_number = len(sys.argv)
        self.topology='topology.pdb'
        self.bases=None
        self.pairs=None
        for i in xrange(0, arguments_number):
            try:
                if (sys.argv[i] == "-top"):
                    self.topology = sys.argv[i + 1]
                if (sys.argv[i] == "-bases"):
                    self.bases = sys.argv[i + 1]
                if (sys.argv[i] == "-pairs"):
                    self.pairs = sys.argv[i + 1]
                if (sys.argv[i] == "-o"):
                    self.ndx=sys.argv[i+1]
                elif (sys.argv[i] == "-h"):
                    print self.helpline
            except IndexError:
                print("Ошибка! Не хватает аргументов командной строки для корректного выполнения программы.")
    pass

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

#Аргументы командной строки
arguments = Arguments()
#Чтение и индексация PDB
topology = PDB(file=arguments.topology)
#Чтение выбранных остатков или пар
#Чтение выбранных остатков или пар
if (arguments.bases!=None) and (arguments.pairs==None):
    #Чтение остатков
    residues_file=open(arguments.bases,'r')
    residues=residues_file.readlines()
    residues_file.close()
    residues=map(\
        lambda res: (res[0].strip(),int(res[1])),\
        map(\
            lambda res: res.split(),\
            residues))
elif (arguments.bases==None) and (arguments.pairs!=None):
    #Чтение готовых пар
    residues_file=open(arguments.pairs,'r')
    residues=residues_file.readlines()
    residues_file.close()
    residues=filter(lambda s: s!='\n',residues)
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
#Индексация атомов выбранных плоскостей
indices=map(lambda residue: [residue[0],residue[1],map(lambda atom: topology.index[residue[0]][residue[1]][atom]+1,flat_groups[residue[0]])],residues)
print
#Вывод
output=reduce(lambda x,y: x+y+'\n',map(lambda residue: reduce(lambda x,y: x+y+' ',map(str,residue[2]),'[ {0:s}{1:d} ]\n'.format(residue[0],residue[1])),indices),'')
output_f=open(arguments.ndx,'w')
output_f.write(output)
output_f.close()
