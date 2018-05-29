from ase import Atoms, Atom
from ase.data import covalent_radii, atomic_numbers
from numpy import pi, sin, cos, sqrt, deg2rad, arccos, matrix, multiply, array, rad2deg
from random import choice, randint
from ase.io import write

class Starting_population(object):
    """Generates an unrelaxed intelligent starting population"""

    def __init__(self, atom_numbers, moleculegroup):
        self.atom_numbers = atom_numbers
        self.moleculegroup = moleculegroup

    def __element_dict__(self, element, last_hyb):
        ele_list = {6: ["tetrahedral","trigonal_planar","linear"],
                    7: ["trigonal_pyramidial", "bent", "s"],
                    1: ["s"]}
        bond_choice = choice(ele_list[element])
        if last_hyb in ele_list[element]:
            chance = randint(1,100)
            if chance < 40:
                bond_choice = last_hyb
        if bond_choice == "tetrahedral":
            bond_set = ['tetrahedral', [(0.0, 0.0, 1.0), (0.0, 0.94264149109217832, -0.33380685923377101),
                                        (0.81635168387599166, -0.47165551958557062, -0.33333316528257168),
                                        (-0.81635168387599166, -0.47165551958557062, -0.33333316528257168)]]
        elif bond_choice == "trigonal_planar":
            bond_set = ['trigonal_planar', [(0.0, 0.0, 1.0), (0.0, 0.866025404, -0.5),
                                            (0.0, -0.866025404, -0.5)]]
        elif bond_choice == "linear":
            bond_set = ['linear', [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)]]
        elif bond_choice == "trigonal_pyramidial":
            bond_set = ['trigonal_pyramidial', [(0.0, 0.0, 1.0), (0.0, 0.94264149109217832, -0.33380685923377101),
                                                (0.81635168387599166, -0.47165551958557062, -0.33333316528257168)]]
        elif bond_choice == "bent":
            bond_set = ['bent', [(0.0, 0.0, 1.0), (0.0, 0.94264149109217832, -0.33380685923377101)]]
        elif bond_choice == "s":
            bond_set = ['s', [(0.0, 0.0, 1.0)]]
        else:
            print ('error detected')
        return bond_set


    #Takes the basic geometry obtained in the geometry funtion and rotates the bonds to match the atom it is attached to.
    def __rotation__(self, hybridisation, vector):
        hybridrot = []
        if hybridisation[0] == vector:
            for i in hybridisation:
                vec = (i[0]*(-1),i[1]*(-1),i[2]*(-1))
                hybridrot.append(vec)
        elif hybridisation[0] == (vector[0]*(-1),vector[1]*(-1),vector[2]*(-1)):
            hybridrot = hybridisation
        else:
            initialvect = hybridisation[0]
            u = (initialvect[1]*vector[2]-initialvect[2]*vector[1],
                 initialvect[2]*vector[0]-initialvect[0]*vector[2],
                 initialvect[0]*vector[1]-initialvect[1]*vector[0])
            magu = sqrt(u[0]**2+u[1]**2+u[2]**2)
            normu = (u[0]/magu,u[1]/magu,u[2]/magu)
            angle = arccos(vector[0]*initialvect[0]+vector[1]*initialvect[1]+vector[2]*initialvect[2])
            x = [cos(angle)+(normu[0]**2)*(1-cos(angle)),
                 normu[0]*normu[1]*(1-cos(angle))-normu[2]*sin(angle),
                 normu[0]*normu[2]*(1-cos(angle))+normu[1]*sin(angle)]
            y = [normu[1]*normu[0]*(1-cos(angle))+normu[2]*sin(angle),
                 cos(angle)+(normu[1]**2)*(1-cos(angle)),
                 normu[1]*normu[2]*(1-cos(angle))-normu[0]*sin(angle)]
            z = [normu[2]*normu[0]*(1-cos(angle))-normu[1]*sin(angle),
                 normu[2]*normu[1]*(1-cos(angle))+normu[0]*sin(angle),
                 cos(angle)+(normu[2]**2)*(1-cos(angle))]
            R = matrix((x,y,z))
            for i in hybridisation:
                matr = matrix ([[i[0]],[i[1]],[i[2]]])
                hybx = (R*matr)
                midlist = array(hybx).tolist()
                hyby = ((-1)*midlist[0][0],(-1)*midlist[1][0],(-1)*midlist[2][0])
                magx = sqrt(hyby[0]**2+hyby[1]**2+hyby[2]**2)
                normx = (hyby[0]/magx,hyby[1]/magx,hyby[2]/magx)
                hybridrot.append(normx)
        hybridrot.remove(hybridrot[0])
        return hybridrot


    #Creates molecules to form an initial population for use in genetic algorithm.
    def randomizer(self):
        atom_numbers = self.atom_numbers
        moleculegroup = self.moleculegroup
        hybridlist = []
        molecule = Atoms()
        for i in atom_numbers:
            if len(molecule)!=0:
                unstable = 1
                full_tick = 0
                while unstable != 0:
                    unstable = 0
                    full_tick += 1
                    if full_tick >= 50:
                        print ('unstable attempt, correcting...')
                        return
                    achoose = choice(hybridlist)
                    ticker = 0
                    connected = 0
                    while achoose[1][1] == []:
                        achoose = choice(hybridlist)
                        ticker += 1
                        if ticker >= 2*len(hybridlist):
                            print ('failed')
                            return
                    vectlist = achoose[1][1]
                    vector = choice(vectlist)
                    bondlength = (covalent_radii[achoose[0].number] + covalent_radii[i])*1
                    j = Atom(i, [achoose[0].position[0] + bondlength*vector[0],
                                 achoose[0].position[1] + bondlength*vector[1],
                                 achoose[0].position[2] + bondlength*vector[2]])
                    ticker = 0
                    for a in molecule:
                        close = (covalent_radii[a.number] + covalent_radii[j.number])*0.7
                        dista = sqrt((a.position[0]-j.position[0])**2+(a.position[1]-j.position[1])**2+(a.position[2]-j.position[2])**2)
                        if dista < close:
                            unstable += 1
                            ticker += 1
                            #print 'close atom detecting, correcting...'
                            if ticker >= 2*len(molecule):
                                print ('failed')
                                return
                initialvector = ((-1)*vector[0],(-1)*vector[1],(-1)*vector[2])
                hybridisation = self.__element_dict__(i,achoose[1][0])
                hybridrot = [hybridisation[0],self.__rotation__(hybridisation[1], vector)]
                hybridlist.append([j,hybridrot])
                molecule.append(j)
                for i in hybridlist:
                    if i[1][1]==[]:
                        hybridlist.remove(i)
                if hybridrot[1] != []:
                    hybridlist.remove(achoose)
                    achoose[1].remove(vectlist)
                    vectlist.remove(vector)
                    achoose[1].append(vectlist)
                    hybridlist.append(achoose)
            else:
                first = Atom(i, position = (0,0,0))
                hybridisation = self.__element_dict__(i,'none')
                hybridlist.append([first, hybridisation])
                molecule.append(first)
        moleculegroup.append(molecule)
        print ('number of molecules created:', len(moleculegroup))