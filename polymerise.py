from write_coords import Writer
from read_lammpsdata import ReadLammpsData
from add_bond import MakeBond
from params import Parameterise
import json
import sys
import os.path
import numpy as np


lammps_sim_file       = sys.argv[1]
vdw_defs              = sys.argv[2]
polymerisation_config = sys.argv[3]

lammps_sim = ReadLammpsData( lammps_sim_file )

lammps_sim.vdw_defs = json.load(open( vdw_defs )) 
lammps_sim.vdw_defs = {int(k):v for k,v in lammps_sim.vdw_defs.items()}

config = json.load(open( polymerisation_config ))

# If this is part of an itterative polymerisation, check if new connecions
# have been made so they don't need to be parameterised every time
if os.path.isfile('new_connections.json'):
    new_connections = 'new_connections.json'
else:
    new_connections = {'angles':{},'dihedrals':{}}

b = MakeBond(lammps_sim, 
                config, Nbonds = 20,
                outputfile = 'polyout',
                new_connections=new_connections,
                networkfile = 'nx.dat')
Parameterise(lammps_sim, lammps_sim.vdw_defs, b.new_connections)

json.dump(lammps_sim.vdw_defs,open('vdw.json','w'))

output = Writer(lammps_sim,'polymerised')
output.write_xyz('polymerised.data')
output.write_lammps('polymerised.data')


