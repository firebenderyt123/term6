from classes.BPLA import BPLA
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

# Потрібна висота та  потрібні  координати
H_ = 170
pos = np.array([
    1500,
    0,
    800
])

Vx_ = 10  # м/с значення  поперечної  швидкості
Vymax = 2  # м/с Максимально допустима  швидкість  вертикального  руху
Vz_ = 10  # м/с значення  поперечної  швидкості

# Жорсткість пружини ????
c = np.array([1, 1, 1, 1, 1, 1])

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
eps0 = 21.0 * 2 * pi
eps1 = 10.0 * 2 * pi
eps2 = 10.0 * 2 * pi
eps3 = 10.0 * 2 * pi
eps4 = 10.0 * 2 * pi
eps5 = 17.0 * 2 * pi
omega0 = eps0
omega1 = eps1
omega2 = eps2
omega3 = eps3
omega4 = eps4
omega5 = eps5


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
        I, m, dt, J, J0, Kf, Kf0,
        Km, Km0, p, G, Ks, H_, pos,
        Vx_, Vymax, Vz_, c,
        state
    )
    bpla.launch()

    states = bpla.states

    position = [state[0:3] for state in states]
    pos_labels = ['x', 'y', 'z']

    plt.plot(position, label=pos_labels)
    plt.title('Поведінка координат відносно ІСК під час польоту')
    plt.xlabel('Час, с')
    plt.ylabel('Координата, м')

    plt.legend()
    plt.show()

    # v = [state[3:6] for state in states]
    # print(position)

    # plt.plot(v)
    # plt.title('Поведінка проекцій V1, V2,V3 швидкості відповідно на осі X, Y, Z ІСК під час польоту')  # noqa
    # plt.xlabel('Час, с')
    # plt.ylabel('Проекції вектору швидкості, м/с')
    # plt.show()
