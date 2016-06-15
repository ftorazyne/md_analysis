#! /usr/bin/python -u
#-*- coding: utf-8 -*-
__author__ = 'MGI'

import sys,math
from copy import deepcopy

class Arguments():
    def __init__(self):
        self.helpline='hbond_analysis.py -xpm hbonds.xpm -ndx hbonds.ndx -gro system.gro -o output \n'
        arguments_number = len(sys.argv)
        self.xpm='hbonds.xpm'
        self.ndx='hbonds.ndx'
        self.gro='system.gro'
        self.out='output'
        self.begin=None
        self.end=None
        self.cutoff=None
        for i in xrange(0, arguments_number):
            try:
                if (sys.argv[i] == "-xpm"):
                    self.xpm = sys.argv[i + 1]
                if (sys.argv[i] == "-ndx"):
                    self.ndx = sys.argv[i + 1]
                if (sys.argv[i] == "-gro"):
                    self.gro = sys.argv[i + 1]
                if (sys.argv[i] == "-o"):
                    self.out = sys.argv[i + 1]
                if (sys.argv[i] == "-b"):
                    self.begin = int(sys.argv[i + 1])-1
                if (sys.argv[i] == "-e"):
                    self.end = int(sys.argv[i + 1])
                if (sys.argv[i] == "-cutoff"):
                    self.cutoff = float(sys.argv[i + 1])/100.0
                elif (sys.argv[i] == "-h"):
                    print self.helpline
                    raise SystemExit(0)
            except IndexError:
                print("Ошибка! Не хватает аргументов командной строки для корректного выполнения программы.")
    pass

class Atom():
    def __init__(self,*args):
        #Если type='gmx', то atom_string интерпретируется как описание из файлов топологии GROMACS
        #Если type='pdb', то atom_string интерпретируется как описание из файла координат PDB
        #Если type='gro', то atom_string интерпретируется как описание из файлов координат GROMACS
        self.index=0
        self.type=None
        self.resname=None
        self.resnumber=0
        self.name=None
        self.chargegroup=0
        self.charge=0.0
        self.mass=0.0
        if len(args)>=2:
            if args[1]=='gmx': self.load_gmx_atom(args[0])
            elif args[1]=='pdb': self.load_pdb_atom(args[0])
            elif args[1]=='gro': self.load_gro_atom(args[0])

    def load_gmx_atom(self,gmx_top_atom):
        try:
            atom_string=filter(lambda s: s!='',gmx_top_atom.split())
            self.index=int(atom_string[0])
            self.type=atom_string[1]
            self.resname=atom_string[3]
            self.resnumber=int(atom_string[2])
            self.name=atom_string[4]
            self.chargegroup=int(atom_string[5])
            self.charge=float(atom_string[6])
            self.mass=float(atom_string[7])
        except IndexError:
            print "Не хватает полей в описании атома"
            raise SystemExit(1)
        except ValueError:
            print "Неверный формат описания атома"
            raise SystemExit(1)

    def load_pdb_atom(self,line):
        try:
            self.index=int(line[6:11].strip())
            self.type=line[76:78].strip()
            self.resname=line[17:20].strip()
            self.resnumber=int(line[22:26].strip())
            self.name=line[12:16].strip()
            self.charge=float(line[78:80].strip())
            self.alter_position=line[16]
            self.chain=line[21]
            self.B_factor=float(line[60:66].strip())
            self.occupancy=float(line[54:60].strip())
            self.x=float(line[30:38].strip())
            self.y=float(line[38:46].strip())
            self.z=float(line[46:54].strip())
        except IndexError:
            print "Не хватает полей в описании атома"
            raise SystemExit(1)
        except TypeError:
            print "Неверный формат описания атома"
            raise SystemExit(1)

    def load_gro_atom(self,line):
        try:
            self.resnumber=int(line[0:5].strip())
            self.resname=line[5:10].strip()
            self.name=line[10:15].strip()
            self.index=int(line[15:20].strip())
            self.x=float(line[20:28].strip())
            self.y=float(line[28:36].strip())
            self.z=float(line[36:44].strip())
            self.v_x=float(line[44:52].strip())
            self.v_y=float(line[52:60].strip())
            self.v_z=float(line[60:68].strip())
        except IndexError:
            print "Не хватает полей в описании атома"
            raise SystemExit(1)
        except TypeError:
            print "Неверный формат описания атома"
            raise SystemExit(1)

    pass

class Coordinates():
    def __init__(self,*args):
        self.atoms_quantity=0
        self.atoms=[]
        self.name=''
        self.box_parameters=[]
        if len(args)>=2:
            topology_f=open(args[0],'r')
            topology_lines=topology_f.readlines()
            topology_f.close()
            if args[1]=='gro':
                self.load_gro_coordinates(topology_lines)

    def load_gro_coordinates(self,coordinates):
        self.name=coordinates[0]
        self.atoms_quantity=int(coordinates[1])
        self.box_parameters=coordinates[-1].split()
        coordinates=coordinates[2:-1]
        for line in coordinates:
            self.atoms.append(Atom(line,'gro'))

    pass

class Ndx():
    def __init__(self,*args):
        class Group():
            def __init__(self):
                self.group=set()
                self.groupname=''
            def add(self,elements):
                self.group.update(map(int,elements))
        class Hbonds():
            def __init__(self):
                self.hbonds=[]
            def add(self,elements):
                self.hbonds.append((int(elements[0]),int(elements[1]),int(elements[2])))
        class Donors():
            def __init__(self):
                self.donors=[]
            def add(self,elements):
                self.donors.append(elements)
        self.groups=[Group(),Group()]
        self.donors_H=[Donors(),Donors()]
        self.acceptors=[Group(),Group()]
        self.hbonds=Hbonds()
        if len(args)>=1:
            ndx_f=open(args[0],'r')
            ndx=ndx_f.readlines()
            ndx_f.close()
            current_group=None
            group_selector=-1
            group_name=''
            for line in ndx:
                line=line.strip()
                if line[0]=='[':
                    group_name=line[1:-1].strip()
                    if group_name.split('_')[0]=='donors':
                        current_group=self.donors_H[group_selector]
                    elif group_name.split('_')[0]=='acceptors':
                        current_group=self.acceptors[group_selector]
                    elif group_name.split('_')[0]=='hbonds':
                        current_group=self.hbonds
                    else:
                        group_selector+=1
                        self.groups[group_selector].groupname=group_name
                        current_group=self.groups[group_selector]
                else:
                    current_group.add(line.split())
    pass

class XPM():
    def __init__(self,discrete=True,table=[[0]]):

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

    def loadfromfile(self,filename):
        def block_selection(sequence,blocksize):
            length=len(sequence)
            for i in xrange(0,length,blocksize):
                yield sequence[i:i+blocksize]

        xpm_f=open(filename,'r')
        xpm=xpm_f.readlines()
        xpm_f.close()
        field_counter=0
        color_quantity=0
        raw_table=[]
        #Разбор файла
        for line in xpm:
            selector=line[0]
            if selector=='/':
                keyword=line[2:].split(':')[0].strip()
                if keyword=='title':
                    self.title=line.split('\"')[1]
                elif keyword=='legend':
                    self.legend=line.split('\"')[1]
                elif keyword=='x-label':
                    self.xlabel=line.split('\"')[1]
                elif keyword=='y-label':
                    self.ylabel=line.split('\"')[1]
                elif keyword=='type':
                    type=line.split('\"')[1]
                    if type=='Discrete':
                        self.isdiscrete=True
                    else:
                        self.isdiscrete=False
                elif keyword=='x-axis':
                    self.xtics+=line.split()[2:-1]
                elif keyword=='y-axis':
                    self.ytics+=line.split()[2:-1]
            elif selector=='\"':
                field_counter+=1
                if field_counter==1:
                    line=map(int,line[1:-3].split())
                    self.width=line[0]
                    self.height=line[1]
                    self.char_per_color=line[3]
                    color_quantity=line[2]
                elif 1<field_counter<=1+color_quantity:
                    value=line[1:].split('\"')[2]
                    self.colors.update({value:(line[1:self.char_per_color+1],line[self.char_per_color+5:self.char_per_color+5+7])})
                else:
                    raw_table+=[line.split('\"')[1]]
        #Преобразование данных из таблицы
        raw_table=[raw_table[i] for i in range(self.height-1,-1,-1)]
        self.data=[]
        reverse_transform=dict()
        for value,code in self.colors.iteritems():
            reverse_transform.update({code[0]:value})
        for raw_line in raw_table:
            self.data+=[[reverse_transform[code] for code in block_selection(raw_line,self.char_per_color)]]
        self.data=[[self.data[j][i] for j in xrange(self.height)] for i in xrange(self.width)]

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

class Hbond():
    def __init__(self,coordinates,groups,number):
        self.quantity=0
        self.frequency=0.0
        triade=index.hbonds.hbonds[number]
        #Отнесение атомов
        self.donor=coordinates.atoms[triade[0]-1].name
        self.hydrogen=coordinates.atoms[triade[1]-1].name
        self.acceptor=coordinates.atoms[triade[2]-1].name
        self.donor_residue=coordinates.atoms[triade[0]-1].resname+'-'+str(coordinates.atoms[triade[0]-1].resnumber)
        self.acceptor_residue=coordinates.atoms[triade[2]-1].resname+'-'+str(coordinates.atoms[triade[2]-1].resnumber)
        if triade[0] in groups.groups[0].group:
            self.donor_group=groups.groups[0].groupname
            self.acceptor_group=groups.groups[1].groupname
        else:
            self.donor_group=groups.groups[1].groupname
            self.acceptor_group=groups.groups[0].groupname


arguments=Arguments()
#Загрузка необходимых данных
atoms=Coordinates(arguments.gro,'gro')
index=Ndx(arguments.ndx)
xpm=XPM()
xpm.loadfromfile(arguments.xpm)
if arguments.begin==None:
    arguments.begin=0
if arguments.end==None:
    arguments.end=xpm.width
else:
    if arguments.end>xpm.width:
        arguments.end=xpm.width
        print 'Заданный номер конечного состояния превышает количество доступных состояний. Расчет производится до последнего доступного состояния'

#Именование водородных связей и подготовка данных для анализа
average_hbond_quantity=0.0
standard_deviation=0.0
#Частоты по каждой связи
hbonds=[]
used_frame_quantity=arguments.end-arguments.begin
ret=lambda x: x
selector=lambda s:(s=='Present' and ret(1)) or (s=='None' and ret(0))
for i in xrange(xpm.height):
    hbonds.append(Hbond(atoms,index,i))
    hbonds[-1].quantity=sum([selector(xpm.data[j][i]) for j in xrange(arguments.begin,arguments.end)])
    hbonds[-1].frequency=1.0*hbonds[-1].quantity/used_frame_quantity
#Средняя величина и дисперсия
bonds_per_state=[float(sum(map(selector,xpm.data[i]))) for i in xrange(arguments.begin,arguments.end,1)]
average_hbond_quantity=sum(bonds_per_state)/used_frame_quantity
standard_deviation=math.sqrt(sum(map(lambda x: (x-average_hbond_quantity)**2,bonds_per_state))/(used_frame_quantity-1))
#Фильтрация, если включено
if arguments.cutoff!=None:
    filtered_xpm=map(lambda s: [], xpm.data)
    filtered_hbonds=[]
    for i in xrange(xpm.height):
        if hbonds[i].frequency>=arguments.cutoff:
            for j in xrange(xpm.width):
                filtered_xpm[j].append(xpm.data[j][i])
            filtered_hbonds.append(hbonds[i])
    filtered_xpm=XPM(table=filtered_xpm)
    filtered_xpm.set_colors(xpm.colors)
    filtered_xpm.title=xpm.title
    filtered_xpm.legend=xpm.legend
    filtered_xpm.xlabel=xpm.xlabel
    filtered_xpm.ylabel=xpm.ylabel
    filtered_xpm.xtics=xpm.xtics
    filtered_xpm.ytics=map(str,range(filtered_xpm.height))+['0']
    hbonds=filtered_hbonds
    filtered_xpm.droptofile(arguments.out+'.filtered.xpm')
    xpm=filtered_xpm
#Вывод
#Определение длины имени водородной связи
hbond=hbonds[0]
test_name=hbond.donor_group+'/'+hbond.donor_residue+'/'+hbond.donor+'-'+hbond.hydrogen+'...'+hbond.acceptor_group+'/'+hbond.acceptor_residue+'/'+hbond.acceptor
name_length=len(test_name)+10
#
if arguments.cutoff!=None:
    outf_xpm=open(arguments.out+'.filtered.dat','w')
else:
    outf_xpm=open(arguments.out+'.analyse.dat','w')
outf_xpm.write("#Анализ карты водородных связей {0} \n".format(arguments.xpm))
outf_xpm.write("#Состояний {0:d}\n".format(used_frame_quantity))
outf_xpm.write("#Анализируется участок с {0:d} по {1:d} состояние\n".format(arguments.begin+1,arguments.end))
outf_xpm.write("#Водородных связей {0}\n".format(xpm.height))
outf_xpm.write("#Среднее число водородных связей на состояние {0:.2f}\n".format(average_hbond_quantity))
outf_xpm.write("#Стандартное отклонение числа водородных связей на состояние {0:.2f}\n".format(standard_deviation))
outf_xpm.write("#Встречаемость водородных связей \n")
#outf_xpm.write("   +-------------+--------------------------------+------------+------------------+ \n")
outf_xpm.write("#{0:>11s}  {1:<{width}s}  {2:15s}\n".format('Номер связи','Пара донор-акцептор','Встречаемость, %',width=name_length))
#outf_xpm.write("   +-------------+--------------------------------+------------+------------------+ \n")
#Инициализация массивов строк *.gro и *.ndx
hbond_table=''
counter=0
for hbond in hbonds:
    counter+=1
    hbond_table=hbond_table+" {0:>11d}  {1:<{width}s}  {2:15.2f}\n".format(\
        counter,\
        hbond.donor_group+'/'+hbond.donor_residue+'/'+hbond.donor+'-'+hbond.hydrogen+'...'+hbond.acceptor_group+'/'+hbond.acceptor_residue+'/'+hbond.acceptor,\
        hbond.frequency*100.0,\
        width=name_length)

outf_xpm.write(hbond_table)
#outf_xpm.write("   +-------------+--------------------------------+------------+------------------+ \n")
outf_xpm.close()
