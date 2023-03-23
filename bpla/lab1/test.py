import numpy as np
from bpla import start_bpla

time = 0.1

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
    states = start_bpla(time, state)
    return states[-1][6:9]


direction = 1  # 1 - увеличение, -1 - уменьшение
min_psi, min_teta, min_gamma = f(np.array([
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
]))
prev_psi, prev_teta, prev_gamma = min_psi, min_teta, min_gamma

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

    print(abs(teta), abs(min_teta))

    if abs(teta) < abs(min_teta):
        min_teta = teta
        best_omega5 = omega5 / (2 * np.pi)
        directionChanged = False
    elif abs(teta) > abs(prev_teta):

        if directionChanged:
            break

        direction = -direction
        directionChanged = True

    else:
        omega5 = (omega5 / (2 * np.pi) + eps * direction) * 2 * np.pi

    prev_teta = teta

print("Best Omega:", best_omega5, "Min teta:", min_teta)
