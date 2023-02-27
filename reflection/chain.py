import numpy as np


class Atom:

    def __init__(self, u=0, v=0, m=1):
        self.u = u
        self.v = v
        self.f = 0
        self.m = m

    def __repr__(self) -> str:
        return f'Atom(u={self.u}, v={self.v}, m={self.m}'


class Spring:

    def __init__(self, k):
        self.k = k


class Boundary:
    periodic = 0
    nonperiodic = 1    
    

class Chain:

    def __init__(self):
        self.atoms = {}
        self.springs = {}
        self.boundary = Boundary.nonperiodic
        self.n = 0

    def create_atoms(self, n=10, cell_constant=1):
        self.n = n
        betta = 0.1
        n0 = round(n/2)
        U0 = 1
        k1 = round(n/20)
        Omega = 10
        k = 1
        m = 1
        a = 1
        g1 = 0.5*a * np.sqrt(4*k/m-1)

        for i in range(n):
            B0 = U0*np.exp(-betta**2/2*(i-n0)**2)
            u = B0 * np.sin(k1*i)
            v = -B0 * (Omega * np.cos(k1*i) - betta**2*g1/a*(i-n0) * np.sin(k1*i))
            self.atoms[i] = Atom(u, v)

    def add_pair_springs(self, k=1):
        for i in range(self.n-1):
            self.springs[i+0.5] = Spring(k)

    def evaluate_step(self, dt=0.01):
        u = [atom.u for atom in self.atoms.values()]
        if self.boundary==Boundary.nonperiodic:
            for atom in self.atoms.values():
                atom.f = 0
            for i in range(self.n-1):
                self.atoms[i].f += self.springs[i+0.5].k * (u[i+1] - u[i])
                self.atoms[i+1].f += -self.springs[i+0.5].k * (u[i+1] - u[i])
            for i in range(1, self.n-1):
                self.atoms[i].v += self.atoms[i].f / self.atoms[i].m * dt
                self.atoms[i].u += self.atoms[i].v * dt
    
    
    
