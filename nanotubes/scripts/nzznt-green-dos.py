# Code for obtaining DOS of a N-ZZNT using Green's function recursive method

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def make_green_hamiltonian(N, alpha, beta):
    NN = 4*N
    H = np.identity(NN, dtype=complex) * alpha

    for i in range(NN):
        for j in range(NN):
            if i == j+1 or i == j-1:
                H[i,j] = beta

    for i in range(1, 2*N, 2):
        H[i, NN-i-1] = H[NN-i-1, i] = beta

    H[0, 2*N-1] = H[2*N-1, 0] = beta
    H[2*N, NN-1] = H[NN-1, 2*N] = beta
    return H

def make_green_V(N, beta):
    NN = 4*N
    V = np.zeros((NN, NN), dtype=complex)

    for i in range(0, 2*N, 2):
        V[NN-i-1, i] = beta
    return V

def make_green_dos(H, V, N, E, eta):
    E_points = len(E)
    NN = 4*N
    dos_array = np.zeros(E_points)

    V_A = V.conj().T

    for i in range(E_points):
        if i % 10 == 0:
            print(f'Energy step: {i} / {E_points}')

        G0 = np.linalg.inv((E[i] + 1j*eta)*np.identity(NN) - H)
        SR = G0
        SL = G0
        for n in range(n_max):
            SR = np.linalg.inv(np.identity(NN) - G0 @ V @ SR @ V_A) @ G0
            SL = np.linalg.inv(np.identity(NN) - G0 @ V_A @ SL @ V) @ G0

        g = np.linalg.inv(np.identity(NN) - SR @ V_A @ SL @ V) @ SR

        dos_array[i] = -np.trace(g).imag
    return dos_array


# General parameters

alpha = 0.0
beta = -2.8

eta = 1e-2
n_max = 800

E_points = 500
E = np.linspace(-9.0, 9.0, E_points)

# N = 9 
N = 9
H = make_green_hamiltonian(N, alpha, beta)
V = make_green_V(N, beta)

dos_g_9 = make_green_dos(H, V, N, E, eta)
AG9 = integrate.simps(dos_g_9, E)

# N = 10
N = 10
H = make_green_hamiltonian(N, alpha, beta)
V = make_green_V(N, beta)

dos_g_10 = make_green_dos(H, V, N, E, eta)
AG10 = integrate.simps(dos_g_10, E)

# N = 11
N = 11
H = make_green_hamiltonian(N, alpha, beta)
V = make_green_V(N, beta)

dos_g_11 = make_green_dos(H, V, N, E, eta)
AG11 = integrate.simps(dos_g_11, E)

plt.figure(figsize=(15,5))

plt.subplot(1,3,1)
plt.plot(E, dos_g_9/AG9, c='k')
plt.xlabel('E')
plt.ylabel('DOS')
plt.ylim(0, 0.3)
plt.grid(ls=':')
plt.title(f'N-ZZNT, N = 9')

plt.subplot(1,3,2)
plt.plot(E, dos_g_10/AG10, c='k')
plt.xlabel('E')
plt.ylabel('DOS')
plt.ylim(0, 0.3)
plt.grid(ls=':')
plt.title(f'N-ZZNT, N = 10', c='k')

plt.subplot(1,3,3)
plt.plot(E, dos_g_11/AG11, c='k')
plt.xlabel('E')
plt.ylabel('DOS')
plt.ylim(0, 0.3)
plt.grid(ls=':')
plt.title(f'N-ZZNT, N = 11')

plt.show()
