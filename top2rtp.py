#! /usr/bin/python -u
#-*- coding: utf-8 -*-
__author__ = 'MGI'

import sys

class Arguments():
    def __init__(self):
        self.helpline='top2rtp.py -top topology.top -rtp residue.rtp \n'
        arguments_number = len(sys.argv)
        self.topology='topology.pdb'
        self.rtp='residue.rtp'
        self.pairs=None
        for i in xrange(0, arguments_number):
            try:
                if (sys.argv[i] == "-top"):
                    self.topology = sys.argv[i + 1]
                if (sys.argv[i] == "-rtp"):
                    self.rtp = sys.argv[i + 1]
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

    pass

class Bond():
    def __init__(self,*args):
        bond=[None,None,0,0.0,0.0]
        if len(args)>=1:
            try:
                bond_raw=filter(lambda s: s!='',args[0].split())
                for i in range(5):
                    if i+1<=len(bond_raw):
                        bond[i]=bond_raw[i]
                self.atom1=int(bond[0])
                self.atom2=int(bond[1])
                self.type=int(bond[2])
                self.length=float(bond[3])
                self.energy=float(bond[4])
            except ValueError:
                print "Неверный формат записи связи"
                raise SystemExit(1)
        else:
            self.atom1=None
            self.atom2=None
            self.type=None
            self.length=0.0
            self.energy=0.0

    pass

class Angle():
    def __init__(self,*args):
        bond=[None,None,None,0,0.0,0.0]
        if len(args)>=1:
            try:
                bond_raw=filter(lambda s: s!='',args[0].split())
                for i in range(6):
                    if i+1<=len(bond_raw):
                        bond[i]=bond_raw[i]
                self.atom1=int(bond[0])
                self.atom2=int(bond[1])
                self.atom3=int(bond[2])
                self.type=int(bond[3])
                self.angle=float(bond[4])
                self.energy=float(bond[5])
            except TypeError:
                print "Неверный формат записи валентного угла"
        else:
            self.atom1=None
            self.atom2=None
            self.atom3=None
            self.type=None
            self.angle=0.0
            self.energy=0.0

    pass

class Dihedral():
    def __init__(self,*args):
        bond=[None,None,None,None,0,0.0,0.0,0]
        if len(args)>=1:
            try:
                bond_raw=filter(lambda s: s!='',args[0].split())
                for i in range(8):
                    if i+1<=len(bond_raw):
                        bond[i]=bond_raw[i]
                self.atom1=int(bond[0])
                self.atom2=int(bond[1])
                self.atom3=int(bond[2])
                self.atom4=int(bond[3])
                self.type=int(bond[4])
                self.angle=float(bond[5])
                self.energy=float(bond[6])
                self.periodicity=int(bond[7])
            except TypeError:
                print "Неверный формат записи торсионного угла"
        else:
            self.atom1=None
            self.atom2=None
            self.atom3=None
            self.atom4=None
            self.type=None
            self.angle=0.0
            self.energy=0.0

    pass

class Atomtypes():
    pass

class Moleculetypes():
    pass

class Pairs():
    pass

def empty_function(*args):
    return None

def extract_section_name(line):
    result=line.split(']')
    result=result[0].split('[')
    result=result[1].strip()
    return result

class Topology:
    def __init__(self,*args):
        self.atoms=[]
        self.atomtypes=[]
        self.molecules=[]
        self.moleculetypes=[]
        self.bonds=[]
        self.pairs=[]
        self.angles=[]
        self.dihedrals=[]
        self.impropers=[]
        if len(args)>=2:
            topology_f=open(args[0],'r')
            topology_lines=topology_f.readlines()
            topology_f.close()
            if args[1]=='gmx':
                self.load_gmx_topology(topology_lines)
            elif args[1]=='pdb':
                self.load_pdb_topology(topology_lines)

    def load_pdb_topology(self,topology):
        bonds_raw=[]
        for line in topology:
            if line[:6] == 'CRYST1':
                self.cryst = line
            elif line[:6] == 'ATOM  ' or line[:6] == 'HETATM':
                self.atoms.append(Atom(line,'pdb'))
            elif line[:6] == 'CONECT':
                bonds_raw.append(line)
        bonds_raw = map(\
            lambda line: map(\
                lambda x: int(x), \
                filter(\
                    lambda x: x != '', \
                    line.split())[1:]), \
            self.bonds)
        #Преобразование связей в плоский список пар нумеров атомов
        bonds_raw=list(\
            set(\
                reduce(\
                    lambda x, y: x + y, \
                    map( \
                        lambda x: map( \
                            lambda y: tuple(sorted([x[0], y])), \
                            x[1:]), \
                        self.bonds))))
        self.bonds = map(lambda x: set(x), bonds_raw)

    def load_gmx_topology(self,topology):
        section_handlers={'atoms': (self.atoms,Atom),\
                         'bonds': (self.bonds,Bond),\
                         'angles': (self.angles,Angle), \
                         'dihedrals': (self.dihedrals, Dihedral)}
        sections=set(section_handlers.keys())
        section_handler=''
        for line in topology:
            if line != '\n':
                if line.strip()[0]==';':
                    pass
                elif line.strip()[0]=='[':
                    section_handler=extract_section_name(line)
                else:
                    if section_handler in sections:
                        section_handlers[section_handler][0].append(section_handlers[section_handler][1](line,'gmx'))
        #Индексация атомов
        self.atoms=dict(\
            map(\
                lambda s: (s.index,s),\
                self.atoms))
        #Индексация связей
        self.bonds=dict(\
            map(\
                lambda s: (s.atom1,dict(\
                    map(\
                        lambda x: (x.atom2,x),\
                        filter(\
                            lambda y: y.atom1==s.atom1,self.bonds)))),\
                self.bonds))
        #Индексация валентных углов
        # self.angles=dict(\
        #     map(\
        #         lambda s: (s.atom2,dict(\
        #             map(\
        #                 lambda x: (x.atom1,dict(\
        #                     map(\
        #                         lambda z: (z.atom3,z),\
        #                         filter(\
        #                             lambda q: (q.atom3==x.atom3) and (q.atom1==x.atom1),\
        #                             self.angles)))),\
        #                 filter(\
        #                     lambda y: y.atom2==s.atom2,\
        #                     self.angles)))),\
        #         self.angles))
        # #Индексация торсионных углов
        # self.dihedrals=dict(map(lambda s: (s.atom2,s),self.dihedrals))

    def make_rtp(self,filename):
        keys=self.atoms.keys()
        residue='[ {0:s} ]\n'.format(self.atoms[keys[0]].resname)
        #Преобразование атомов
        atoms=reduce(lambda a1,a2: a1+'  {0:4s}  {1:4s}  {2:8.4f}  {3:6d}\n'.format(a2.name,a2.type,a2.charge,a2.chargegroup),self.atoms.itervalues(),' [ atoms ]\n')
        #Преобразование связей
        bonds=reduce(\
            lambda a1,a2: a1+reduce(\
                lambda p,q: p+'  {0:4s}  {1:4s}  {2:8.4f}  {3:12.2f}\n'.format(self.atoms[self.bonds[a2][q].atom1].name,self.atoms[self.bonds[a2][q].atom2].name,self.bonds[a2][q].length,self.bonds[a2][q].energy),\
                sorted(self.bonds[a2].keys()),\
                ''),\
            sorted(self.bonds.keys()),\
            ' [ bonds ]\n')
        #Преобразование углов
        angles=reduce(\
            lambda a1,a2: a1+'  {0:4s}  {1:4s}  {2:4s}  {3:8.3f}  {4:12.4f}\n'.format(self.atoms[a2.atom1].name,self.atoms[a2.atom2].name,self.atoms[a2.atom3].name,a2.angle,a2.energy),\
            self.angles,\
            ' [ angles ]\n')
        #Преобразование двугранников
        dihedrals='[ dihedrals ]\n'
        impropers='[ impropers ]\n'
        for dihedral in self.dihedrals:
            if dihedral.type==9:
                dihedrals+='  {0:4s}  {1:4s}  {2:4s}  {3:4s}  {4:8.3f}  {5:12.4f} {6:1d}\n'.format(self.atoms[dihedral.atom1].name,self.atoms[dihedral.atom2].name,self.atoms[dihedral.atom3].name,self.atoms[dihedral.atom4].name,dihedral.angle,dihedral.energy,dihedral.periodicity)
            if dihedral.type==4:
                impropers+='  {0:4s}  {1:4s}  {2:4s}  {3:4s}  {4:8.3f}  {5:12.4f} {6:1d}\n'.format(self.atoms[dihedral.atom1].name,self.atoms[dihedral.atom2].name,self.atoms[dihedral.atom3].name,self.atoms[dihedral.atom4].name,dihedral.angle,dihedral.energy,dihedral.periodicity)
        rtp_f=open(filename,'w')
        rtp_f.write(residue+atoms+bonds+angles+dihedrals+impropers)
        rtp_f.close()

    pass

arguments=Arguments()
topology=Topology(arguments.topology,'gmx')
topology.make_rtp(arguments.rtp)
print