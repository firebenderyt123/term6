from scipy.integrate import odeint
import numpy as np


class BPLA:

    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])
    e3 = np.array([0, 0, 1])

    def __init__(
            self,
            I, m, dt, J, J0, Kf, Kf0,  # noqa
            Km, Km0, p, G, dH,
            state
            ):
        self.t = 0.0
        self.I = I  # noqa
        self.m = m
        self.dt = dt
        self.J = J
        self.J0 = J0
        self.Kf = Kf
        self.Kf0 = Kf0
        self.Km = Km
        self.Km0 = Km0
        self.p = p
        self.G = G
        self.dH = dH
        self.state = state  # стан беспілотника

        # Рівнодіючs силb тяги двигунів у проекціях на осі ІСК
        self.Fe = np.array([0.0, 0.0, 0.0])

        # Проекції сили опору повітря
        self.Fs = np.array([0.0, 0.0, 0.0])

        # Вектор гіроскопічного моменту в проекціях на осі ЗСК
        self.Mg = np.array([0.0, 0.0, 0.0])

        # Cумарний момент сил тяги двигунів
        self.MFe = np.array([0.0, 0.0, 0.0])

        # Cумарний аеродинамічний момент від шести двигунів
        self.MAe = np.array([0.0, 0.0, 0.0])

        self.next_step()

    # Функція, що повертає похідні від вектора стану
    def state_dot(self, state, t):
        x, y, z, vx, vy, vz, psi, teta, gamma, Om1, Om2, Om3, eps0, eps1, eps2, eps3, eps4, eps5 = state  # noqa

        # розрахунок похідних
        x_dot = vx
        y_dot = vy
        z_dot = vz
        v_dot = self.Fe / self.m + self.G + self.Fs / self.m

        psi_dot = (Om3 * np.sin(gamma) - Om2 * np.cos(gamma)) / np.cos(teta) # noqa
        teta_dot = Om2 * np.sin(gamma) + Om3 * np.cos(gamma)
        gamma_dot = Om1 + np.tan(teta) * (Om3 * np.sin(gamma) - Om3 * np.cos(gamma)) # noqa

        Om_dot = (-self.Mg - self.dH + self.MFe + self.MAe) / self.I

        omega0_dot = eps0
        omega1_dot = eps1
        omega2_dot = eps2
        omega3_dot = eps3
        omega4_dot = eps4
        omega5_dot = eps5

        return np.array([
            x_dot,
            y_dot,
            z_dot,
            v_dot[0],
            v_dot[1],
            v_dot[2],
            psi_dot,
            teta_dot,
            gamma_dot,
            Om_dot[0],
            Om_dot[1],
            Om_dot[2],
            omega0_dot,
            omega1_dot,
            omega2_dot,
            omega3_dot,
            omega4_dot,
            omega5_dot
        ])

    def next_step(self):
        self.t += self.dt

        # чисельне інтегрування
        new_state_dot = BPLA.euler_integration(self.state_dot, self.state, self.dt, self.t + self.dt)[0]  # noqa

        test_state_dot = odeint(self.state_dot, self.state, [self.t + self.dt])[0]  # noqa

        if test_state_dot.all() == new_state_dot.all():
            print('equal')
        else:
            print('not equal', test_state_dot)

        print(new_state_dot)

    def calc_new_state(self):
        x, y, z, vx, vy, vz, psi, teta, gamma, Om1, Om2, Om3, omega0, omega1, omega2, omega3, omega4, omega5 = self.state # noqa

        # 1
        # Вектори сили тяги двигунів
        F = np.array([
            self.Kf0 * self.omega0 ** 2 * self.e2,
            self.Kf * self.omega1 ** 2 * self.e2,
            self.Kf * self.omega2 ** 2 * self.e2,
            self.Kf * self.omega3 ** 2 * self.e2,
            self.Kf * self.omega4 ** 2 * self.e2,
            self.Kf * self.omega5 ** 2 * self.e3
        ])
        Fez = np.sum(F)  # Сумарний вектор

        # кватерніон орієнтації ЗСК відносно ІСК
        q_psi = np.array([np.cos(psi/2), 0, -np.sin(psi/2), 0])
        q_teta = np.array([np.cos(teta/2), 0, 0, np.sin(teta/2)])
        q_gamma = np.array([np.cos(gamma/2), np.sin(gamma/2), 0, 0])
        q = q_psi @ q_teta @ q_gamma

        # 4
        Phiez = np.hstack([0, Fez])
        self.Fe = q * Phiez * np.conj(q)

        # 5
        # vx =
        # Fsick =

        self.Fs = [
            -self.Ks[0] * vx ** 2 * np.sign(vx),
            -self.Ks[1] * vy ** 2 * np.sign(vy),
            -self.Ks[2] * vz ** 2 * np.sign(vz)
        ]

        # 6

        # Вектор сумарного кінетичного моменту «МК + ротори з гвинтами»
        H = [
            self.I[0] * Om1,
            self.I[1] * Om2 + self.J(omega1 - omega2 + omega3 - omega4) + self.J0 * omega0,  # noqa
            self.I[2] * Om3 + self.J * omega5
        ]

        self.Mg = np.cross(np.array([Om1, Om2, Om3]), H)  # векторний доб

        # 7
        self.dH = [
            0,
            self.J * (self.eps[1] - self.eps[2] + self.eps[3] - self.eps[4]) + self.J[0] * self.eps[0],  # noqa
            self.J * self.eps[5]
        ]

        # 8
        self.MFe = np.sum(np.cross(self.p, F), axis=0)

        # 9

        # Аеродинамічні моменти від кожного гвинта
        Ma = np.array([
            self.Km0 * omega0 ** 2 * self.e2,
            self.Km * omega1 ** 2 * self.e2,
            -self.Km * omega2 ** 2 * self.e2,
            self.Km * omega3 ** 2 * self.e2,
            -self.Km * omega4 ** 2 * self.e2,
            self.Km * omega5 ** 2 * self.e3
        ])

        # 10
        self.MAe = np.sum(Ma)

    @staticmethod
    def euler_integration(f, x0, h, xn):
        """
        f: функція правої частини диференціального рівняння
        x0: початкове значення змінної x
        h: крок інтегрування
        xn: кінцеве значення змінної x
        return: масив значень y на кожному кроці інтегрування
        """
        n = int((xn - x0[0]) / h) + 1
        y = np.zeros((n, len(x0)))
        y[0] = x0
        x = x0[0]

        for i in range(1, n):
            y[i] = y[i-1] + h * f(y[i-1], x)
            x += h

        return y
