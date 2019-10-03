import numpy as np
import random, time
from dolfin import *
import math
print("Hello, I am the T matrix generator")


# ------------------
# --- PARAMETERS ---
# ------------------
outdir = './data/modelF/'
kdir = './data/perm/'
dd = 8 #16#8#4
fv = 1.0e3
def getf(i,j):
    cx = (i+0.5)*hh
    cy = (j+0.5)*hh;
    if cy < hh: return fv
    if cy > (1.0-hh): return -fv
    return 0.0

# -----------------
# --- LOAD k(x) ---
# -----------------
NS = 10
NN = 64
Nf = NN*NS
for ii in range(NS):
    locN = ii*NS + 0
    infile0 = kdir+'/k'+str(locN)+'.txt'
    with open(infile0) as f:
        lines_list = f.readlines()
        my_data0 = [float(val) for val in lines_list[1::2]]# with scipping
    arr0 = np.reshape(my_data0, (NN, NN))
    for jj in range(1, NS):
        locN = ii*NS + jj
        infileij = kdir+'/k'+str(locN)+'.txt'
        with open(infileij) as f:
            lines_list = f.readlines()
            my_dataij = [float(val) for val in lines_list[1::2]]# with scipping
        arrij = np.reshape(my_dataij, (NN, NN))
        arr0 = np.concatenate((arr0, arrij), axis=1)# concatinate col-wise
    if ii==0:
        arrK = arr0 # concatinate row-wise
    else:
        arrK = np.concatenate((arrK, arr0))# concatinate col-wise
print(arrK.min(), arrK.max())
Nx, Ny = arrK.shape
print(Nx, Ny)


# --------------------
# --- RESCALE k(x) ---
# --------------------
# arrK2 = np.log(arrK)
arrK2 = arrK[::dd, ::dd]
Nx2, Ny2 = arrK2.shape

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
fig = plt.figure(figsize=(10, 4))

ax1 = fig.add_subplot(121)
im1 = ax1.imshow(arrK, interpolation='None')
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im1, cax=cax1, orientation='vertical')
major_ticks = np.arange(0, Nx, NN)
ax1.set_xticks(major_ticks)
ax1.set_yticks(major_ticks)

ax2 = fig.add_subplot(122)
im2 = ax2.imshow(arrK2, interpolation='None')
divider = make_axes_locatable(ax2)
cax2 = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im2, cax=cax2, orientation='vertical');
major_ticks = np.arange(0, Nx2, NN/dd)
ax2.set_xticks(major_ticks)
ax2.set_yticks(major_ticks)

plt.show()
plt.close()
print(Nx, Ny, NN)
print(Nx2, Ny2, NN/dd)

NN = int(NN/dd)
Nf = NN*NS
arrK = arrK[::dd, ::dd]
Nx, Ny = arrK.shape
print(Nx, Ny, NN)


# ------------------------
# --- GENERATE T and S ---
# ------------------------
import sys, math
import petsc4py
from petsc4py import PETSc
petsc4py.init(sys.argv)

hh = 1.0/Nx; volK = hh*hh
n = Nx*Ny

T = PETSc.Mat().createAIJ([n, n], nnz=5)
vecS = PETSc.Vec().createSeq(n) 
for i in range(Nx):
    for j in range(Ny):
        I = i*Ny+j
        diagval = 0
        if j!=0:
            val = 2.0/(1.0/arrK[i,j-1]+1.0/arrK[i,j]); diagval += val
            T.setValue(I, I-1, -val)
        if j!=(Ny-1):
            val = 2.0/(1.0/arrK[i,j+1]+1.0/arrK[i,j]); diagval += val
            T.setValue(I, I+1, -val)
        if i!=0:
            val = 2.0/(1.0/arrK[i-1,j]+1.0/arrK[i,j]); diagval += val
            T.setValue(I, I-Ny, -val)
        if (i!=(Nx-1)):
            val = 2.0/(1.0/arrK[i+1,j]+1.0/arrK[i,j]); diagval += val
            T.setValue(I, I+Ny, -val)
        T.setValue(I, I, diagval)
        # S for GMsFEM
        sval = arrK[i,j]*volK
        vecS.setValue(I, sval) 
T.assemblyBegin()
T.assemblyEnd()
print('generate T')


# -------------------------------
# --- SAVE T, S, RHS and DOFs ---
# -------------------------------
# save T
filenameT = outdir + 'mat-K.txt'
fileT = open(filenameT, "w")
bufferT = ''
for I in range(n):
    cols,vals = T.getRow(I)
    for cj in range(len(cols)):
        bufferT += str(I) + ' ' + str(cols[cj]) + ' ' + str(vals[cj]) + '\n'
fileT.write(bufferT)
fileT.close()
print('save mat T into ' + filenameT)

# save S as vec
filenameS = outdir + 'mat-Sdiag.txt'
fileS = open(filenameS, "w")
bufferS = ''
for I in range(n):
    sval = vecS.getValue(I)
    bufferS += str(I) + ' ' + str(sval) + '\n'
fileS.write(bufferS)
fileS.close()
print('save mat S into ' + filenameS)

# save DOF
filenameDof = outdir + 'dof100'
fileDof = open(filenameDof, "w")
bufferDof = ''
for I in range(n):
    dof = I; ci = I
    mci = 0; pi = 0
    bufferDof += str(dof) + ' ' + str(ci) + ' ' + str(mci) + ' ' + str(pi) + ' ' + str(volK) + '\n'
fileDof.write(bufferDof)
fileDof.close()
print('save DOF into ' + filenameDof)

q = PETSc.Vec().createSeq(n) 
for i in range(Nx):
    for j in range(Ny):
        I = i*Ny+j
        q.setValue(I, getf(i,j)*volK)
        
# save Rhs
outRhs = outdir + 'rhs.txt'
fileRhs = open(outRhs, "w")
bufferRhs = ''
for I in range(n):
    sval = q.getValue(I)
    bufferRhs += str(I) + ' ' + str(sval) + '\n'
fileRhs.write(bufferRhs)
fileRhs.close()
print('save rhs into ' + outRhs)        

print('size ', n)


# ------------------------
# --- SAVE figure of k ---
# ------------------------
meshc = UnitSquareMesh(NS, NS)
Vc = FunctionSpace(meshc, 'DG', 0)
uc = Function(Vc)
uarrc = uc.vector().array()

mesh = UnitSquareMesh(Nx, Ny)
V = FunctionSpace(mesh, 'DG', 0)
u = Function(V)
uarr = u.vector().array()
print("functions fenics")

for i in range(Nx):
    for j in range(Ny):
        I = i*Ny + j
        val = arrK[i,j]
        uarr[2*I] = val
        uarr[2*I+1] = uarr[2*I]

u.vector().set_local(uarr)
filef = File(outdir+"results/k.pvd")
filef << u 
print('k saved')

