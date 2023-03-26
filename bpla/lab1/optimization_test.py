from scipy.optimize import minimize
import numpy as np
from bpla import start_bpla

pos = np.array([
    0.0,
    50.0,
    0.0
])
v = np.array([
    0.0,
    0.0,
    0.0
])
psi = 0.0
teta = 0.0
gamma = 0.0
Om = np.array([0.0, 0.0, 0.0])  # вектор кутової швидкості МК
eps = np.array([
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
])
omega = np.array([
    21.0 * 2 * np.pi,
    10.0 * 2 * np.pi,
    10.0 * 2 * np.pi,
    10.0 * 2 * np.pi,
    10.0 * 2 * np.pi,
    17.0 * 2 * np.pi
])


# def f(x):
#     return x[0]**2 + x[1] ** 2


def f(x):
    # if x[0] < 21 * 2 * np.pi:
    #     x[0] = 100000000000
    # if x[1] < 17 * 2 * np.pi:
    #     x[1] = 100000000000
    # omega = np.array([
    #     0.0 * 2 * np.pi,
    #     x[0],
    #     x[0],
    #     x[0],
    #     x[0],
    #     0.0 * 2 * np.pi
    # ])
    omega = np.array([
        x[0],
        297.1363333,
        297.1363333,
        297.1363333,
        297.1363333,
        x[1]
    ])
    state = [
        pos, v,
        psi, teta, gamma,
        Om, omega, eps
    ]
    states = start_bpla(0.04, state)
    x, y, z = states[-1][0:3]
    p, t, g = states[-1][6:9]
    Om1, Om2, Om3 = states[-1][9:12]
    # print(p, t, g)
    # return abs(p) + abs(t) + abs(g)
    # return abs(y - 50)
    return -y


if __name__ == '__main__':
    init_val = f([omega[0], omega[5]])
    print(init_val)

    initial_guess = np.array([omega[0], omega[5]])
    result = minimize(
        f,
        initial_guess,
        method='Nelder-Mead'
    )
    print(result)
    print(result.x / 2 / np.pi)
