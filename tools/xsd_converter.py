import re
import sys

"""A tool put together for extracting the coordinates out of an xsd format file
since ASE doesn;t always correctly handle the scaling and units.  The lattice output
doesn't carry over the lattice vectors but puts 0's in as a place holder.
Output is in extended xyz format.
Input is the name of the file without the '.xsd' suffix.
"""

file = sys.argv[1]
with open(file+'.xsd') as xsdfile:
   data = xsdfile.read()
   atom3d = re.findall('Atom3d*.*', data)

coords = []
atoms = []
for i in atom3d:
    coord_single = re.findall('XYZ="[0-9\.\,-]*',i)
    atom_single = re.findall('Components="[a-zA-Z]*"',i)
    if len(atom_single) == 1:
        atom_single1 = atom_single[0].strip('Components=')
        atom_single2 = atom_single1.strip('"')
        atoms.append(atom_single2)
    if (len(coord_single)==1):
        coord_string = coord_single[0].strip('XYZ="')
        coord_string2 = coord_string.split(',')
        coords.append(coord_string2)
    else:
        pass



with open(file+'.xyz','a') as output:
    output.write(str(len(atoms))+'\n')
    output.write('Lattice=" 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"\n')
    for j in range(len(atoms)):
        output.write(atoms[j]+'   '+coords[j][0]+'  '+coords[j][1]+'  '+coords[j][2]+'\n')
