import matplotlib.pyplot as plt
import numpy as np
from celluloid import Camera
from progress.bar import IncrementalBar

def tanh(x, a):
    return 0.5+0.5*(np.exp(a*x)-np.exp(-a*x))/(np.exp(a*x)+np.exp(-a*x))

n = 1500
dt = 0.01
max_time = 1500
snap_steps = 1000

m1 = 1
m2 = 1
U0 = 1
beta = 0.02
n0 = 50+4/beta
a = 1
g1 = 0
c1 = 1
c2 = 5
d1 = 1
d2 = 1
omega_min1 = np.sqrt(d1/m1)
omega_min2 = np.sqrt(d2/m2)
omega_max1 = np.sqrt((4*c1+d1)/m1)
omega_max2 = np.sqrt((4*c2+d2)/m2)
omega_min = max(omega_min1, omega_min2)
omega_max = min(omega_max1, omega_max2)
Omega = 0.45 * (omega_max+omega_min)
print(Omega)
k1 = np.arccos(0.5*m1/c1) * ((2*c1 + d1)/m1 - Omega**2)
g1 = 0.5/Omega * np.sqrt((Omega**2 - d1/m1) * (4*c1/m1 + d1/m1 - Omega**2))

dt = 0.05/omega_max
max_steps = round(max_time / dt)
 
A_reflected = []
A_transmitted = []
A_initial = []


#for index, theta in enumerate([0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]):
for index, theta in enumerate([0.1]):

    A_reflected.append(0)
    A_transmitted.append(0)
    A_initial.append(0)


    m = np.ones(n)
    u = np.zeros(n)
    v = np.zeros(n)

    for i in range(n):
        Bn = U0 * np.exp(-0.5 * beta**2 * (i-n0)**2)
        u[i] = Bn * np.sin(k1 * i)
        v[i] = -Bn * ( Omega * np.cos(k1 * i) -  
                      beta**2 * g1/a * (i-n0) * np.sin(k1 * i))

    c = np.ones(n-1)
    d = np.ones(n-1)

    for i in range(n-1):
        c[i] = c1+(c2-c1)*tanh(a*(i-0.5*n), theta)
        d[i] = d1+(d2-d1)*tanh(a*(i-0.5*n), theta)

    figure = plt.figure(figsize=(24, 6))
    ax = plt.subplot(1, 1, 1)
    ax2 = ax.twinx()
    ax.set_ylim((-2*U0, 2*U0))
    camera = Camera(figure)
    ax2.plot(list(range(n-1)), c, color='r', linewidth=1)
    ax.plot(list(range(n)), u, marker='o', color='b', markersize=2, linewidth=1)
    camera.snap()
    
    
    bar = IncrementalBar('Calculating', max=max_steps)

    for step in range(max_steps):

        if step==20000:
            u[:200] = 0
            v[:200] = 0

        # if step < 20000:
        #     A_initial[index] = max(A_initial[index], max(u[:3*n0]))

        # if step > 40000:
        #     A_transmitted[index] = max(A_transmitted[index], max(u[3*n0:])) 

        # if step > 40000 and step < 50000:
        #     A_reflected[index] = max(A_reflected[index], max(u[:3*n0])) 

        u_old = np.copy(u)
        for i in range(1, n-1):
            v[i] += ( -c[i-1] * (u_old[i] - u_old[i-1]) 
                     - c[i] * (u_old[i] - u_old[i+1]) 
                     - d[i]*u_old[i]) / m[i] * dt
            u[i] += v[i]*dt
        
        if step % snap_steps == 0:
            ax2.plot(list(range(n-1)), c, color='r', linewidth=1)
            ax.plot(list(range(n)), u, marker='o', color='b', markersize=2, linewidth=1, label=f'{step}')
            ax.text(5, 0.5, f'{step}', fontsize=12)
            camera.snap()
        bar.next()

    bar.finish()
    print()

    animation = camera.animate()
    animation.save(f'test_{theta}.gif')

print(A_initial)
print(A_reflected)
print(A_transmitted)