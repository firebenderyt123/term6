import matplotlib.pyplot as plt
from math import pi
import numpy as np

# Головні моменти інерції МК
I = np.array([  # noqa
    1.8,
    2.8,
    3.0
])
m = 18.0  # маса БПЛА, кг
dt = 0.01  # часовий такт оновлення значень
J = 0.00016  # момент інерції рухомої частини електричногодвигуна (кг∙м2)
J0 = 0.04  # Момент інерції рухомою частини ДВЗ
Kf = 0.0005  # Коефіцієнт силової характеристики ЕД
Kf0 = 0.01  # Коефіцієнт силової характеристики ДВЗ
Km = 0.0001  # Коефіцієнт моментної характеристики ЕД
Km0 = 0.0004  # Коефіцієнт моментної характеристики ДВз
p = np.array([
    np.array([0, 0, 0]),
    np.array([0.5, 0, -0.5]),
    np.array([0.5, 0, 0.5]),
    np.array([-0.5, 0, 0.5]),
    np.array([-0.5, 0, -0.5]),
    np.array([-1.5, 0, 0]),
])  # радіус-вектори розташування ДВЗ та ЕД в ЗСК (м)

# Силя тяжіння
G = np.array([
    0.0,
    -9.81,
    0.0
])

# Экспериментально отримані коефіцієнти опору повітря
Ks = np.array([
    0.1,
    0.1,
    0.1
])

# Початкові умови

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
omega0 = 0.0
omega1 = 0.0
omega2 = 0.0
omega3 = 0.0
omega4 = 0.0
omega5 = 0.0
eps0 = 21.0 * 2 * pi
eps1 = 10.0 * 2 * pi
eps2 = 10.0 * 2 * pi
eps3 = 10.0 * 2 * pi
eps4 = 10.0 * 2 * pi
eps5 = 17.0 * 2 * pi


def start_bpla():
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

    bpla = BPLA(
        I, m, dt, J, J0, Kf, Kf0, Km, Km0, p, G, Ks,
        state
    )
    bpla.launch()

    states = bpla.states

    v = [state[3:6] for state in states]
    print(v, len(v))

    plt.plot(v)
    plt.title('Найкраще значення функції на кожній ітерації')
    plt.xlabel('Час, с')
    plt.ylabel('Проекції вектору швидкості, м/с')
    plt.show()
