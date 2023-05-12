import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import quaternion  # noqa

time = 100

# Головні моменти інерції МК
Ixx = 1.8
Iyy = 2.8
Izz = 3.0

m = 18.0  # маса БПЛА, кг
dt = 0.01  # часовий такт оновлення значень
J = 0.00016  # момент інерції рухомої частини електричногодвигуна (кг∙м2)
J0 = 0.04  # Момент інерції рухомою частини ДВЗ
Kf = 0.0005  # Коефіцієнт силової характеристики ЕД
Kf0 = 0.01  # Коефіцієнт силової характеристики ДВЗ
Km = 0.0001  # Коефіцієнт моментної характеристики ЕД
Km0 = 0.0004  # Коефіцієнт моментної характеристики ДВз
p = np.array([
    np.array([0.0, 0.0, 0.0]),
    np.array([0.5, 0.0, -0.5]),
    np.array([0.5, 0.0, 0.5]),
    np.array([-0.5, 0.0, 0.5]),
    np.array([-0.5, 0.0, -0.5]),
    np.array([-1.5, 0.0, 0.0]),
])  # радіус-вектори розташування ДВЗ та ЕД в ЗСК (м)

# Силя тяжіння
Gx = 0.0
Gy = -9.81
Gz = 0.0

# Экспериментально отримані коефіцієнти опору повітря
Ks = np.array([
    0.1,
    0.1,
    0.1
])

c = np.array([Kf0, Kf, Kf, Kf, Kf, Kf])

# Початкові умови

x, y, z = 0.0, 50.0, 0.0
vx, vy, vz = 0.0, 0.0, 0.0

psi = 0.0
teta = 0.0
gamma = 0.0
Om1, Om2, Om3 = 0.0, 0.0, 0.0
eps0, eps1, eps2, eps3, eps4, eps5 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

omega0 = 20.67 * 2 * np.pi
omega1 = 10 * 2 * np.pi
omega2 = 10 * 2 * np.pi
omega3 = 10 * 2 * np.pi
omega4 = 10 * 2 * np.pi
omega5 = 16.5 * 2 * np.pi

e1 = np.array([1, 0, 0])
e2 = np.array([0, 1, 0])
e3 = np.array([0, 0, 1])


# program
states = []

num_steps = int(time / dt)
for i in range(num_steps):

    currStep = i

    # вивід поточного стану кожну секунду
    if (i % (num_steps / time) == 0):
        t = i / (num_steps / time)
    print(f'{t} [{x} {y} {z}] [{psi}, {teta}, {gamma}]')  # noqa
    states.append([x, y, z, vx, vy, vz, psi, teta, gamma, Om1, Om2, Om3, omega0, omega1, omega2, omega3, omega4, omega5, eps0, eps1, eps2, eps3, eps4, eps5])  # noqa

    if np.isnan(x):
        break

    # Алгоритм використання моделі (1)-(12)
    # 1
    # Вектори сили тяги двигунів
    F = np.array([
        Kf0 * omega0 ** 2 * e2,
        Kf * omega1 ** 2 * e2,
        Kf * omega2 ** 2 * e2,
        Kf * omega3 ** 2 * e2,
        Kf * omega4 ** 2 * e2,
        Kf * omega5 ** 2 * e3
    ])

    Fe_zsk = np.sum(F, axis=0)  # Сумарний вектор

    # кватерніон орієнтації ЗСК відносно ІСК
    Lam_psi = np.quaternion(
        np.cos(psi/2),
        0,
        np.sin(psi/2),
        0
    )
    Lam_teta = np.quaternion(
        np.cos(teta/2),
        0,
        0,
        np.sin(teta/2)
    )
    Lam_gamma = np.quaternion(
        np.cos(gamma/2),
        np.sin(gamma/2),
        0,
        0
    )
    Lam = Lam_psi * Lam_teta * Lam_gamma

    # 4
    Phiez = np.quaternion(
        0,
        Fe_zsk[0],
        Fe_zsk[1],
        Fe_zsk[2]
    )
    Fe = (Lam * Phiez * Lam.inverse()).vec

    # 5
    v = np.quaternion(0, vx, vy, vz)
    v_zsk = (Lam.inverse() * v * Lam).vec

    Fs_zsk = np.array([
        -Ks[0] * v_zsk[0] ** 2 * np.sign(v_zsk[0]),
        -Ks[1] * v_zsk[1] ** 2 * np.sign(v_zsk[1]),
        -Ks[2] * v_zsk[2] ** 2 * np.sign(v_zsk[2])
    ])
    Phis_zsk = np.quaternion(
        0,
        Fs_zsk[0],
        Fs_zsk[1],
        Fs_zsk[2]
    )

    Fs = (Lam * Phis_zsk * Lam.inverse()).vec

    # 6
    # Вектор сумарного кінетичного моменту «МК + ротори з гвинтами»
    H = np.array([
        Ixx * Om1,
        Iyy * Om2 + J * (
            omega1 - omega2 + omega3 - omega4
        ) + J0 * omega0,
        Izz * Om3 + J * omega5
    ])

    Mg_zsk = np.cross(np.array([Om1, Om2, Om3]), H)  # векторний доб

    # 7
    dH = np.array([
        0,
        J * (
            eps1 - eps2 + eps3 - eps4
        ) + J0 * eps0,
        J * eps5
    ])

    # 8
    MFe = np.sum(np.cross(p, F), axis=0)

    # 9
    # Аеродинамічні моменти від кожного гвинта
    Ma = np.array([
        Km0 * omega0 ** 2 * e2,
        Km * omega1 ** 2 * e2,
        -Km * omega2 ** 2 * e2,
        Km * omega3 ** 2 * e2,
        -Km * omega4 ** 2 * e2,
        Km * omega5 ** 2 * e3
    ])

    # 10
    MAe = np.sum(Ma, axis=0)

    # На кожному кроці інтегрування обчислювати праві частини (управління)
    # диференціальних рівнянь (13)-(18) послідовно за формулами (32), (33)
    b = 0.1
    k0 = -b ** 2
    k1 = -2 * b

    dR_isk = np.quaternion(0, x - 0, 0, z - 0)

    dR_zsk = (Lam.inverse() * dR_isk * Lam).vec

    gamma_star = (m * (k1 * v_zsk[2] + k0 * dR_zsk[2]) - Fs_zsk[2] + c[5] * omega5**2) / (-Gy * m)  # noqa

    teta_star = -(m * (k1 * v_zsk[0] + k0 * dR_zsk[0]) - Fs_zsk[0]) / (-Gy * m)

    # print(gamma_star, teta_star)
    # exit()

    # (23), (25), (27), (28) з роз'ясненнями під формулами;  (19) - (22).
    # 23
    vx_ = Fe[0] / m + Gx + Fs[0] / m
    vy_ = Fe[1] / m + Gx + Fs[1] / m
    vz_ = Fe[2] / m + Gx + Fs[2] / m

    b = 0.5
    k0 = -b ** 3
    k1 = -3 * b ** 2
    k2 = -3 * b

    ay_star = k2 * vy_ + k1 * vy + k0 * (y - 50)

    # 25
    b = 5
    k0 = -b ** 3
    k1 = -3 * b ** 2
    k2 = -3 * b

    psi_ = (Om3 * np.sin(gamma) - Om2 * np.cos(gamma)) / np.cos(teta)
    teta_ = Om2 * np.sin(gamma) + Om3 * np.cos(gamma)
    gamma_ = Om1 + np.tan(teta) * (Om3 * np.sin(gamma) - Om2 * np.cos(gamma))

    Om1_ = (-Mg_zsk[0] - dH[0] + MFe[0] + MAe[0]) / Ixx
    Om2_ = (-Mg_zsk[1] - dH[1] + MFe[1] + MAe[1]) / Iyy
    Om3_ = (-Mg_zsk[2] - dH[2] + MFe[2] + MAe[2]) / Izz

    psi__ = (Om3_ * np.sin(gamma) + Om3 * np.cos(gamma) * gamma_ - Om2_ * np.cos(gamma) + Om2 * np.sin(gamma) * gamma_) / np.cos(teta) + (Om3 * np.sin(gamma) - Om2 * np.cos(gamma)) * np.sin(teta) * teta_ / np.cos(teta)**2  # noqa

    sigma_y_star = k2 * psi__ + k1 * psi_ + k0 * (psi - 0)

    # 27, 28
    b = 15
    k0 = -b ** 3
    k1 = -3 * b ** 2
    k2 = -3 * b

    teta__ = Om2_ * np.sin(gamma) + Om3_ * np.cos(gamma) + gamma_ * (Om2 * np.cos(gamma) - Om3 * np.sin(gamma)) # noqa

    gamma__ = Om1_ + teta_ * (
        Om3 * np.sin(gamma) - Om2 * np.cos(gamma)
    ) / np.cos(teta) ** 2 + np.tan(teta) * (
        Om3_ * np.sin(gamma) - Om2_ * np.cos(gamma)
    ) + np.tan(teta) * gamma_ * (
        Om3 * np.cos(gamma) + Om2 * np.sin(gamma)
    )

    sigma_z_star = k2 * teta__ + k1 * teta_ + k0 * (teta - teta_star)

    sigma_x_star = k2 * gamma__ + k1 * gamma_ + k0 * (gamma - gamma_star)

    # 19 - 22
    eps0 = (ay_star * m) / (2 * c[0] * omega0)
    eps5 = (sigma_y_star * Iyy) / (3 * c[5] * omega5)

    if (currStep % 2 == 0):
        eps1 = (sigma_z_star * Izz) / (c[1] * omega1 + c[2] * omega2 + c[3] * omega3 + c[4] * omega4) # noqa
        eps2 = eps1
        eps3 = -eps1
        eps4 = -eps1
    else:
        eps1 = (sigma_x_star * Ixx) / (c[1] * omega1 + c[2] * omega2 + c[3] * omega3 + c[4] * omega4) # noqa
        eps2 = -eps1
        eps3 = -eps1
        eps4 = eps1

    # print(currStep % 2 == 0, sigma_z_star, eps1, omega1, omega2)
    # exit()

    # Введемо змінні стану МК та  керуючі змінні.
    # 1 - 12
    x_dot = vx
    y_dot = vy
    z_dot = vz

    vx_dot = Fe[0] / m + Gx + Fs[0] / m
    vy_dot = Fe[1] / m + Gy + Fs[1] / m
    vz_dot = Fe[2] / m + Gz + Fs[2] / m

    psi_dot = (Om3 * np.sin(gamma) - Om2 * np.cos(gamma)) / np.cos(teta)
    teta_dot = Om2 * np.sin(gamma) + Om3 * np.cos(gamma)
    gamma_dot = Om1 + np.tan(teta) * (Om3 * np.sin(gamma)-Om2 * np.cos(gamma))

    Om1_dot = (-Mg_zsk[0] - dH[0] + MFe[0] + MAe[0]) / Ixx
    Om2_dot = (-Mg_zsk[1] - dH[1] + MFe[1] + MAe[1]) / Iyy
    Om3_dot = (-Mg_zsk[2] - dH[2] + MFe[2] + MAe[2]) / Izz

    x = x + x_dot * dt
    y = y + y_dot * dt
    z = z + z_dot * dt

    vx = vx + vx_dot * dt
    vy = vy + vy_dot * dt
    vz = vz + vz_dot * dt

    psi = psi + psi_dot * dt
    teta = teta + teta_dot * dt
    gamma = gamma + gamma_dot * dt

    Om1 = Om1 + Om1_dot * dt
    Om2 = Om2 + Om2_dot * dt
    Om3 = Om3 + Om3_dot * dt

    # 13 - 18
    omega0_dot = eps0
    omega1_dot = eps1
    omega2_dot = eps2
    omega3_dot = eps3
    omega4_dot = eps4
    omega5_dot = eps5

    omega0 = omega0 + omega0_dot * dt
    omega1 = omega1 + omega1_dot * dt
    omega2 = omega2 + omega2_dot * dt
    omega3 = omega3 + omega3_dot * dt
    omega4 = omega4 + omega4_dot * dt
    omega5 = omega5 + omega5_dot * dt


# save to excel
titles = [
    'x',
    'y',
    'z',
    'Vx',
    'Vy',
    'Vz',
    'psi',
    'teta',
    'gamma',
    'Om1',
    'Om2',
    'Om3',
    'omega0',
    'omega1',
    'omega2',
    'omega3',
    'omega4',
    'omega5',
    'eps0',
    'eps1',
    'eps2',
    'eps3',
    'eps4',
    'eps5'
]

df = pd.DataFrame(states, columns=titles)

df.to_excel('out/states.xlsx', index=False)

# draw plot
position = [state[0:3] for state in states]
pos_labels = ['x', 'y', 'z']

v = [state[3:6] for state in states]
v_labels = ['Vx', 'Vy', 'Vz']

angle = [np.rad2deg(state[6:9]) for state in states]
angle_labels = ['psi', 'teta', 'gamma']

Om = [np.rad2deg(state[9:12]) for state in states]
Om_labels = ['Om1', 'Om2', 'Om3']

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

axs[0, 0].plot(position, label=pos_labels)
axs[0, 0].set_title('Поведінка координат відносно ІСК під час польоту')
axs[0, 0].set_xlabel('Час, мс')
axs[0, 0].set_ylabel('Координата, м')
axs[0, 0].legend()

axs[0, 1].plot(v, label=v_labels)
axs[0, 1].set_title('Поведінка проекцій швидкості відповідно ІСК')
axs[0, 1].set_xlabel('Час, мс')
axs[0, 1].set_ylabel('Проекції вектору швидкості, м/с')
axs[0, 1].legend()

axs[1, 0].plot(angle, label=angle_labels)
axs[1, 0].set_title('Поведінка кутів відповідно ІСК')
axs[1, 0].set_xlabel('Час, мс')
axs[1, 0].set_ylabel('Кут тангажу та крену, град')
axs[1, 0].legend()

axs[1, 1].plot(Om, label=Om_labels)
axs[1, 1].set_title('Поведінка кутових швидкостей відповідно ІСК')
axs[1, 1].set_xlabel('Час, мс')
axs[1, 1].set_ylabel('Кутова швидкість, град/с')
axs[1, 1].legend()

plt.show()
