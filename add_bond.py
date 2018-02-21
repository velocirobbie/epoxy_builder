import numpy as np
import time
from write_coords import Writer

class MakeBond(object):
    def __init__(self, sim, atom_type1, atom_type2, change_data,
            search_radius = 7, outputfile=0, networkfile=0, Nbonds=20):
        self.sim = sim
        self.a1  = atom_type1
        self.a2  = atom_type2
        self.change_data = change_data
        self.r2 = search_radius**2 
        self.networkfile = networkfile

        searching = True
        count = 0
        while searching:
            self.ids2index = np.zeros((max(sim.ids+1)),dtype=int)
            for i in range(len(sim.ids)):
                self.ids2index[sim.ids[i]] = i

            #self.ids2index = {}
            #for i, id_ in enumerate(sim.ids):
            #    self.ids2index[id_] = i
            
            print '===new search===='
            
            found = False
            all1, xyz1 = self.find_all_atoms_with_type(self.a1)
            all2, xyz2 = self.find_all_atoms_with_type(self.a2)
            if len(all1) == 0 or len(all2) == 0: break
            distances = self.calc_distance_matrix(xyz1,xyz2)
            min_dist = distances.min()
            if min_dist < self.r2:
                i,j = np.where( distances == min_dist )
                a,b = all1[int(i)], all2[int(j)]
                #if (self.sim.atom_labels[a] != self.a1 and
                #        self.sim.atom_labels[b] != self.b1): raise Exception
                #print sim.atom_labels[a], sim.atom_labels[b], a,b 
                #print sim.molecules[a], sim.molecules[b]
                #print self.find_neighbours(a), self.find_neighbours(b)
                self.make_bond(self.sim,a,b)
                count += 1
                found = True
                
            if count >= Nbonds: searching = False
            if not found: searching = False
            """
            self.ids2index = np.zeros((max(sim.ids+1)),dtype=int)
            for i in range(len(sim.ids)):
                self.ids2index[sim.ids[i]] = i


            def distance2(xyz1,xyz2):
                bx = self.sim.xhi - self.sim.xlo
                by = self.sim.yhi - self.sim.ylo
                bz = self.sim.zhi - self.sim.zlo
                def check_periodic2(vec, length):
                    half_length = length / 2
                    over = vec > half_length
                    under= vec < -half_length
                    return vec + (under * length) - (over * length)
        
                x = xyz2[0] - xyz1[0]
                x = check_periodic2(x, bx)
                y = xyz2[1] - xyz1[1]
                y = check_periodic2(y, by)
                z = xyz2[2] - xyz1[2]
                z = check_periodic2(z, bz)
                return np.sqrt( x*x + y*y + z*z )

            for i, bond in enumerate(self.sim.bonds):
                a1 = self.ids2index[bond[0]]
                a2 = self.ids2index[bond[1]]
                d = distance2(self.sim.coords[a1],self.sim.coords[a2])
                #if bond[0] == 54939 and bond[1] == 54941:
                #    print bond, self.sim.coords[54939], self.sim.coords[54941]
                if d > 10:
                    print i,a1,a2,bond,d
                    print self.sim.coords[a1]
                    print self.sim.coords[a2]        
            """
        print atom_type1, atom_type2, 'Bonds made: ', count
        if outputfile:
            with open(outputfile,'a') as f:
                f.write(str(atom_type1)+' '+ str(atom_type2)+'\t Bonds made: '
                        + str(count) + '\n')

    def make_bond(self,sim,a,b):
        change_data = self.change_data
        
        typea = sim.atom_labels[a]
        typeb = sim.atom_labels[b]
        sim.bonds = np.vstack((sim.bonds,[sim.ids[a],sim.ids[b]]))
        print "new bond: index ",a,b,', labels',self.sim.ids[a],self.sim.ids[b]
        if self.networkfile:
            with open(self.networkfile,'a') as f:
                f.write(str(a)+'\t'+str(b)+'\t'+
                        str(sim.molecules[a])+'\t'+
                        str(sim.molecules[b])+'\t'+
                        str(typea)+'\t'+str(typeb)+'\n')

        sim.bond_labels = np.append(sim.bond_labels,change_data['new_bond'])
        
        sim.atom_labels[a] = change_data['atoms'][str(typea)]['label']
        sim.charges[a] = change_data['atoms'][str(typea)]['charge']
        sim.atom_labels[b] = change_data['atoms'][str(typeb)]['label']
        sim.charges[b] = change_data['atoms'][str(typeb)]['charge']
         
        aneighbours = self.find_neighbours(a)
        bneighbours = self.find_neighbours(b)        
        neighbours = aneighbours + bneighbours
        #if len(aneighbours) != 5: raise Exception
        #if len(bneighbours) != 4: raise Exception
        
        if change_data['neighbours']:
          for n in neighbours:
            typen = sim.atom_labels[n]
            if str(typen) in change_data['neighbours']:
                sim.atom_labels[n] =change_data['neighbours'][str(typen)]['label']
                sim.charges[n] = change_data['neighbours'][str(typen)]['charge']
        
        self.add_new_angles(a,b,aneighbours,bneighbours)
        self.add_new_dihedrals(a,b,aneighbours,bneighbours)

        self.update_connection_labels(sim, a,b)
        #if sim.molecules[a] != sim.molecules[b]:
        #    mol_to_change = sim.molecules[b]
        #    for i in range(len(sim.molecules)):
        #        if sim.molecules[i] == mol_to_change:
        #            sim.molecules[i] = sim.molecules[a]
        
        ha = self.find_hydrogen_on(sim,a)
        hb = self.find_hydrogen_on(sim,b)
        #print 'removing',ha,hb,self.sim.ids[ha],self.sim.ids[hb]
        self.remove_atoms(sim, hb,ha)
        """ 	
        for i in range(len(sim.dihedrals)):
            if sim.dihedral_labels[i] == 0:
                print sim.dihedrals[i]
                for atom in sim.dihedrals[i]:
                    print atom, sim.atom_labels[atom-1]
                raise Exception
        for i in range(len(sim.angles)):
            if sim.angle_labels[i] == 0:
                print sim.angles[i]
                for atom in sim.angles[i]:
                    print atom, sim.atom_labels[atom-1]
                raise Exception
	    """
    def update_connection_labels(self, sim, a, b):
        for thing in ['angle','dihedral']:
            qlist = thing+'s'
            qlabels = thing+'_labels'

            connections  = set(np.where((
                           getattr(self.sim,qlist)== sim.ids[a]) )[0])
            connections |= set(np.where((
                           getattr(self.sim,qlist)== sim.ids[b]) )[0])
            for connection in connections:
                types = []
                for atom in getattr(self.sim,qlist)[connection]:
                    types += [self.sim.atom_labels[self.ids2index[atom]]]
                found = 0
                for datum in self.change_data[qlist]:
                    data_connection = self.change_data[qlist][datum]
                    if ( types == data_connection 
                         or list(reversed(types)) == data_connection ):
                        found +=1
                        label = getattr(self.sim,qlabels)#[connection]
                        label[connection] = datum
                if found != 1 and getattr(self.sim,qlabels)[connection] == 0:
                    pass #print 'NOPE',qlist,types

    def add_new_angles(self,a,b,aneighbours,bneighbours):
        for i in set(aneighbours) - {b}:
            self.sim.angles = np.vstack((self.sim.angles, [self.sim.ids[i],
                                                           self.sim.ids[a],
                                                           self.sim.ids[b]]))
            self.sim.angle_labels = np.append(self.sim.angle_labels,0)
        for i in set(bneighbours) - {a}:
            self.sim.angles = np.vstack((self.sim.angles, [self.sim.ids[a],
                                                           self.sim.ids[b],
                                                           self.sim.ids[i]]))
            self.sim.angle_labels = np.append(self.sim.angle_labels,0)

    def add_new_dihedrals(self,a,b,aneighbours,bneighbours):
        def add(dihedral):
            self.sim.dihedrals = np.vstack((self.sim.dihedrals,
                                 [self.sim.ids[index] for index in dihedral] ))
            self.sim.dihedral_labels = np.append(self.sim.dihedral_labels, 0)

        aneighbours = set(aneighbours) - {b}
        bneighbours = set(bneighbours) - {a}
        #add torsions with a,b in the centre
        for i in aneighbours:
            for j in bneighbours:
                add([i,a,b,j])
        #add torsions with b at one end
        for i in aneighbours:
            ineighbours = set( self.find_neighbours(i) ) - {a}
            for j in ineighbours:
                add([j,i,a,b])
        #add torsions with a at one end
        for i in bneighbours:
            ineighbours = set( self.find_neighbours(i) ) - {b}
            for j in ineighbours:
                add([j,i,b,a])

    def find_hydrogen_on(self, sim, centre):
        neighbours = self.find_neighbours(centre)
        found = False
        for neighbour in neighbours:
            #print 'H search',centre,':',neighbour,sim.masses[sim.atom_labels[neighbour]],sim.atom_labels[neighbour]
            if sim.masses[sim.atom_labels[neighbour]] == 1.008:
                #self.remove_atom(sim, neighbour)
                found = True
                break
        if not found: raise Exception(centre, neighbours)
        return neighbour

    def remove_atoms(self, sim, a, b):
        sim.coords = np.delete(sim.coords,[a,b],0)
        
        #Remove connectivity
        connectivity_quantities = ['bond','angle','dihedral','improper']
        for q in connectivity_quantities:
          for i in [a,b]:
            qlist = q+'s'
            qlabels = q+'_labels'
            connections = np.where( getattr(sim, qlist)==sim.ids[i])
            new_array = np.delete( getattr(sim, qlist),connections[0], 0)
            setattr(sim, qlist, new_array)
            new_array = np.delete( getattr(sim, qlabels), connections[0])
            setattr(sim, qlabels, new_array)

        atom_quantities = ['molecules','atom_labels','charges','ids']
        for q in atom_quantities:
            new_array = np.delete( getattr(sim, q),[a,b])
            setattr( sim, q, new_array)
        #Squish connectivity labels
        
        #def squish(x):
        #    if x > atom: x -= 1
        #    return x
        #squish = lambda x: x-1 if x > atom else x
        #vsquish = np.vectorize(squish)
        """
        connection_types = ['bonds','angles','dihedrals','impropers']
        for connection_type in connection_types:
            connections = getattr(sim, connection_type)
            connections -= np.array(connections > a)
            connections -= np.array(connections > b)
        """ 
    def find_neighbours(self, centre):
        #print '----- find_neighbours', centre,self.sim.atom_labels[centre],'-----'
        #print 'id',self.sim.ids[centre]
        bonds = np.where(self.sim.bonds==self.sim.ids[centre])
        neighbours = []
        for bond in np.transpose(bonds):
            neighbours += [self.ids2index[self.sim.bonds[bond[0],bond[1]-1]] ]
        #print 'centre:',centre,', neighbours:',neighbours
        for a in neighbours:
            bonds = np.where(self.sim.bonds==self.sim.ids[a])
            worked = False
            for bond in np.transpose(bonds):
                b = self.ids2index[self.sim.bonds[bond[0],bond[1]-1]]
                bond1 = self.sim.bonds[bond[0],bond[1]]
                bond2 = self.sim.bonds[bond[0],bond[1]-1] 
                #print a, b, [bond1,bond2], [self.ids2index[bond1], self.ids2index[bond2]]
                if b == centre: worked += 1
            if worked != 1: raise Exception(centre,neighbours)
        return neighbours
    
    def distance_between_atoms(self,i,j):
        posi = self.sim.coords[i]
        posj = self.sim.coords[j]
        r = 0
        for k in range(3):
            r += ( posi[k] - posj[k] ) **2
        return r

    def calc_distance_matrix(self,xyz1,xyz2):
        bx = self.sim.xhi - self.sim.xlo
        by = self.sim.yhi - self.sim.ylo
        bz = self.sim.zhi - self.sim.zlo
        def check_periodic(vec, length):
            half_length = length / 2
            over = np.array(vec > half_length)
            under= np.array(vec < -half_length)
            return vec + (under * length) - (over * length)
            
            #for i in range(len(vec)):
            #    for j in range(len(vec[i])):
            #        if vec[i,j] > half_length:
            #            vec[i,j] = vec[i,j] - length
            #        if vec[i,j] < -half_length:
            #            vec[i,j] = vec[i,j] + length
            #return vec

        x = xyz2[:,0] - xyz1[:,0,np.newaxis]
        x = check_periodic(x, bx)
        y = xyz2[:,1] - xyz1[:,1,np.newaxis]
        y = check_periodic(y, by)
        z = xyz2[:,2] - xyz1[:,2,np.newaxis]
        z = check_periodic(z, bz)
        return x*x + y*y + z*z

    #def calc_distance_matrix(self,all1,all2):
    #    x = len(all1); y = len(all2)
    #    distances = np.zeros((x,y))
    #    for i in range(x):
    #        for j in range(y):
    #            distances[i,j] = self.distance_between_atoms(all1[i],all2[j])
    #    return distances

    def find_all_atoms_with_type(self, label):
        atoms = []
        xyz   = []
        for i in range(len(self.sim.atom_labels)):
            if self.sim.atom_labels[i] == label:
                atoms += [i]
                xyz   += [self.sim.coords[i]]
        return np.array(atoms), np.array(xyz)

