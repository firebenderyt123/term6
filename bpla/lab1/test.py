import numpy as np
from bpla import start_bpla

X0 = 0.0
Y0 = 50.0
Z0 = 0.0
Vx0 = 0.0
Vy0 = 0.0
Vz0 = 0.0

psi0 = 0.0
teta0 = 0.0
gamma0 = 0.0
Om1 = 0.0  # вектор кутової швидкості МК
Om2 = 0.0  # вектор кутової швидкості МК
Om3 = 0.0  # вектор кутової швидкості МК
eps0 = 0
eps1 = 0
eps2 = 0
eps3 = 0
eps4 = 0
eps5 = 0
omega0 = 21.0 * 2 * np.pi
omega1 = 10.0 * 2 * np.pi
omega2 = 10.0 * 2 * np.pi
omega3 = 10.0 * 2 * np.pi
omega4 = 10.0 * 2 * np.pi
omega5 = 16.0 * 2 * np.pi


def f(state):
    states = start_bpla(state)
    return states[-1][6:9]


direction = 1  # 1 - увеличение, -1 - уменьшение
prev_psi = -5.19423180e-04
min_psi = -5.19423180e-04
best_omega5 = omega5 / (2 * np.pi)
eps = 0.01  # шаг

directionChanged = False

while True:

    state = np.array([
        X0, Y0, Z0,
        Vx0, Vy0, Vz0,
        psi0, teta0, gamma0,
        Om1, Om2, Om3,
        omega0,
        omega1,
        omega2,
        omega3,
        omega4,
        omega5,
        eps0,
        eps1,
        eps2,
        eps3,
        eps4,
        eps5
    ])

    psi, teta, gamma = f(state)

    print(psi, prev_psi)

    if abs(psi) < abs(min_psi):
        min_psi = psi
        best_omega5 = omega5 / (2 * np.pi)
        print(omega5)
    elif abs(psi) > abs(prev_psi):

        if directionChanged:
            break

        direction = -direction
        directionChanged = True

    else:
        omega5 = (omega5 / (2 * np.pi) + eps * direction) * 2 * np.pi

    prev_psi = psi

print(best_omega5, min_psi)
