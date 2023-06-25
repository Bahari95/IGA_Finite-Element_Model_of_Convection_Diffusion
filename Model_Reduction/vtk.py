from evtk.hl import gridToVTK 
import numpy as np 
import random as rnd 

# Coordinates
x = np.load("mesh_X.npy")
y = np.load("mesh_Y.npy")
z = np.load("mesh_Z.npy")
# We add Jacobian function to make the grid more interesting
Solut_un = np.load("un_sol.npy")
Solut_dx = np.load("ux.npy")
# Dimensions 
nx, ny, nz = x.shape[0], y.shape[0], z.shape[0]
ncells = nx * ny * nz
npoints = (nx + 1) * (ny + 1) * (nz + 1)

# Variables <br>
gridToVTK("./domain", x, y, z, pointData = {"un_solution" : Solut_un, "dx_solution" : Solut_dx})

