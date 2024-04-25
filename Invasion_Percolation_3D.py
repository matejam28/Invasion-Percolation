from math import inf
import pandas as pd
from ecl.eclfile import EclFile
from ecl.grid import EclGrid
import numpy as np


def is_neighbour(grid, root_index, index):
    root_dx, root_dy, root_dz = grid.get_ijk(root_index)
    row_dx, row_dy, row_dz = grid.get_ijk(index)

    delta_dx = abs(root_dx - row_dx)
    delta_dy = abs(root_dy - row_dy)
    delta_dz = abs(root_dz - row_dz)

    if delta_dx == 1 and delta_dy == 0 and delta_dz == 0:
        return True
    if delta_dx == 0 and delta_dy == 1 and delta_dz == 0:
        return True
    if delta_dx == 0 and delta_dy == 0 and delta_dz == 1:
        return True

    return False

def findNeigh(grid,listCheck,listAlreadyIn):
    findNeighs=[]
    for index in listCheck:
        root_dx, root_dy, root_dz = grid.get_ijk(index)
        for dx in (root_dx+1,root_dx-1):
            for dy in (root_dy+1,root_dy-1):
                for dz in (root_dz+1,root_dz-1):
                    checkIndex=grid.findIndex(dx,dy,dz)
                    if checkIndex not in listCheck and checkIndex not in listAlreadyIn:
                        findNeighs.append(index)
    return findNeighs


init=EclFile('3DMODEL64.INIT')
fields=['PERMX', 'PERMY', 'PERMZ', 'PORO', 'DX', 'DY', 'DZ', 'DEPTH']
df = pd.DataFrame(zip(*[init.iget_named_kw(field, 0) for field in fields]), columns=fields)
convmDtom2=9.869E-16
df['RADIUS'] = np.sqrt(df['PERMX']*convmDtom2/df['PORO'])
df['SIGMA'] = 0.07
df['GRAVITY'] = -9.819
df['DELTARO'] = 200
df['Pt'] = 2*(df['SIGMA']/df['RADIUS'])-(df['DELTARO']*df['GRAVITY']*df['DEPTH'])

grid=EclGrid('3DMODEL64.EGRID')
poroKw=init.iget_named_kw('PORO',0)
grid3d=grid.create_3d(poroKw)
grid3d[:]=0

root_index = 0
root = df.loc[root_index]
min_pressure = inf
min_index = root_index


def find_neighbours(grid, filled_cells):
    neighbours = []

    for cell in filled_cells:
        x, y, z = grid.get_ijk(cell)
        for dx in (x+1, x-1):
            if dx >= grid.nx:
                continue
            if dx < 0:
                continue
            index = grid.get_global_index(ijk=(dx, y, z))
            if index not in filled_cells and index not in neighbours:
                neighbours.append(index) # not invaded cells next to the invaded cells (boundary of the IP CO2 plum migration)
        for dy in (y+1, y-1):
            if dy >= grid.ny:
                continue
            if dy < 0:
                continue
            index = grid.get_global_index(ijk=(x, dy, z))
            if index not in filled_cells and index not in neighbours:
                neighbours.append(index)
        for dz in (z+1, z-1):
            if dz >= grid.nz:
                continue
            if dz < 0:
                continue
            index = grid.get_global_index(ijk=(x, y, dz))
            if index not in filled_cells and index not in neighbours:
                neighbours.append(index)
    return neighbours


root = 223199
print(grid.get_ijk(root))
filled_cells = [root]

for iter in range(0, 200):
    min_pressure = inf
    min_index = root_index

    neighbours = find_neighbours(grid, filled_cells)
    for index in neighbours:
        pressure = df['Pt'][index]
        if pressure < min_pressure:
            min_pressure = pressure
            min_index = index
    filled_cells.append(min_index)
    grid3d[grid.get_ijk(min_index)]=iter

grid3d=grid3d.astype(int)
grid3d=np.swapaxes(grid3d,0,2)
#print(grid3d)
np.savetxt('invasiongrid3d.txt',grid3d[:,31,:],fmt='%i')
#print(filled_cells, len(filled_cells))
print(list(map(lambda x: grid.get_ijk(x), filled_cells)))
#print(grid3d)

import matplotlib.pyplot as plt
#plt.scatter([0,1],[0,2])
plt.title("Invasion Percolation: Simple Model")
plt.legend(neighbours)
plt.imshow(grid3d[:,31,:])
plt.show()

#for index, row in df.iterrows():
#     break
#     print (grid.get_ijk(index))
#     if index > 100000:
#        break
#     continue

#     if is_neighbour(grid, root_index, index):
#         print(f'Node {index} ({grid.get_ijk(index)}) is a neighbour of {root_index} ({grid.get_ijk(root_index)})')
#         pressure = row['Pc']

#         if pressure < min_pressure:
#             min_pressure = pressure
#             min_index = index

# print(f'Result is node {df.loc[min_index]} with pressure {min_pressure}')


df.head()
df.to_csv(r'IP3D.csv')