from write_coords import Writer
from read_lammpsdata import ReadLammpsData
from add_bond import MakeBond
import yaml
import sys
with open('polyout','a') as f:
    f.write('New epoxy resin formation\n')

a = ReadLammpsData(sys.argv[1])#'reactants.data')
epoxy_amine2 = yaml.load(open('epoxy_close.json'))
b = MakeBond(a,1,4, epoxy_amine2,outputfile = 'polyout')

"""
c = ['bonds','angles','dihedrals']
for attr in c:
    print 'testing', thing
    thing = getattr(a,attr)
    for i in thing:
        sum_ = 0
        for j in range(3):
            sum_ += (a.coords[c,j] - a.coords[d,j])**2
        if sum_ > 100:
            print attr,i,sum_
"""

output = Writer(a)
output.write_xyz()
output.write_lammps()
