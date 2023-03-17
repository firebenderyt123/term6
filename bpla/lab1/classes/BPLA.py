import numpy as np
from math import sqrt


class BPLA:

    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])
    e3 = np.array([0, 0, 1])

    def __init__(
            self,
            I, m, dt, J, J0, Kf, Kf0,  # noqa
            Km, Km0, p, G, Ks,
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
        self.Ks = Ks
        self.state = state  # стан беспілотника

        # Вектор локальної похідної кінетичного моменту двигунів
        self.dH = np.array([0.0, 0.0, 0.0])

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

    def launch(self):
        np.seterr(all='raise')
        for i in range(int(50 / self.dt)):  # 50 cек
            self.next_step()
            self.calc_new_state()
            print(self.t)

    # Функція, що повертає похідні від вектора стану
    def state_dot(self, state):
        x, y, z, vx, vy, vz, psi, teta, gamma, Om1, Om2, Om3, omega0, omega1, omega2, omega3, omega4, omega5, eps0, eps1, eps2, eps3, eps4, eps5 = state  # noqa

        # розрахунок похідних
        x_dot = vx
        y_dot = vy
        z_dot = vz
        v_dot = self.Fe / self.m + self.G + self.Fs / self.m

        psi_dot = (Om3 * np.sin(gamma) - Om2 * np.cos(gamma)) / np.cos(teta) # noqa
        teta_dot = Om2 * np.sin(gamma) + Om3 * np.cos(gamma)
        gamma_dot = Om1 + np.tan(teta) * (Om3 * np.sin(gamma) - Om2 * np.cos(gamma)) # noqa

        Om_dot = (-self.Mg - self.dH + self.MFe + self.MAe) / self.I

        omega0_dot = eps0
        omega1_dot = eps1
        omega2_dot = eps2
        omega3_dot = eps3
        omega4_dot = eps4
        omega5_dot = eps5

        print(Om_dot, self.Mg, self.dH, self.MFe, self.MAe, self.I)

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
            omega5_dot,
            eps0,
            eps1,
            eps2,
            eps3,
            eps4,
            eps5
        ], dtype=np.float64)

    def next_step(self):
        self.t += self.dt

        # чисельне інтегрування
        new_state_dot = BPLA.euler_integration(self.state_dot, self.state, [self.t, self.t + self.dt], self.dt)  # noqa

        self.prev_state = self.state
        self.state = new_state_dot

    def calc_new_state(self):
        x, y, z, vx, vy, vz, psi, teta, gamma, Om1, Om2, Om3, omega0, omega1, omega2, omega3, omega4, omega5, eps0, eps1, eps2, eps3, eps4, eps5 = self.state # noqa

        # 1
        # Вектори сили тяги двигунів
        F = np.array([
            self.Kf0 * omega0 ** 2 * self.e2,
            self.Kf * omega1 ** 2 * self.e2,
            self.Kf * omega2 ** 2 * self.e2,
            self.Kf * omega3 ** 2 * self.e2,
            self.Kf * omega4 ** 2 * self.e2,
            self.Kf * omega5 ** 2 * self.e3
        ])

        Fe_zsk = np.sum(F, axis=0)  # Сумарний вектор

        # кватерніон орієнтації ЗСК відносно ІСК
        Lam_psi = np.array([
            np.cos(psi/2),
            0,
            -np.sin(psi/2),
            0
        ])
        Lam_teta = np.array([
            np.cos(teta/2),
            0,
            0,
            np.sin(teta/2)
        ])
        Lam_gamma = np.array([
            np.cos(gamma/2),
            np.sin(gamma/2),
            0,
            0
        ])

        # Lam1 = np.multiply(Lam_psi, Lam_teta, Lam_gamma)

        Lam = BPLA.quat_mult(
            BPLA.quat_mult(Lam_psi, Lam_teta),
            Lam_gamma
        )

        # 4
        Phiez = np.hstack([0, Fe_zsk])
        Lam_inv = BPLA.quat_inverse(Lam)
        self.Fe = BPLA.vect(  # ???
            BPLA.quat_mult(
                BPLA.quat_mult(Lam, Phiez),
                Lam_inv
             )
        )

        # 5
        v_zsk = BPLA.quat_mult(
            BPLA.quat_mult(
                Lam_inv,
                np.array([0, vx, vy, vz])
            ),
            Lam
         )

        Fs_zsk = np.array([
            0,  # ???
            -self.Ks[0] * v_zsk[0] ** 2 * np.sign(v_zsk[0]),
            -self.Ks[1] * v_zsk[1] ** 2 * np.sign(v_zsk[1]),
            -self.Ks[2] * v_zsk[2] ** 2 * np.sign(v_zsk[2])
        ])
        self.Fs = BPLA.vect(  # ???
            BPLA.quat_mult(
                BPLA.quat_mult(Lam, Fs_zsk),
                Lam_inv
            )
        )

        # 6

        # Вектор сумарного кінетичного моменту «МК + ротори з гвинтами»
        H = np.array([
            self.I[0] * Om1,
            self.I[1] * Om2 + self.J * (omega1 - omega2 + omega3 - omega4) + self.J0 * omega0,  # noqa
            self.I[2] * Om3 + self.J * omega5
        ])

        Om = np.array([Om1, Om2, Om3])

        self.Mg = np.cross(Om, H)  # векторний доб

        print(Om, H)

        # 7
        self.dH = np.array([
            0,
            self.J * (eps1 - eps2 + eps3 - eps4) + self.J0 * eps0,  # noqa
            self.J * eps5
        ])

        # 8 ??? двигунів 5, а сил 6
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
        self.MAe = np.sum(Ma, axis=0)

    @staticmethod
    def vect(q):  # ???
        """
        q: кватерніон
        return: 3-х мірний вектор
        """
        return q[1:]

    @staticmethod
    def quat_mult(a, b):
        """
        Множення кватерніонів
        a: кватерніон
        b: кватерніон
        return: кватерніон
        """
        q0 = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
        q1 = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2]
        q2 = a[0] * b[2] + a[2] * b[0] + a[3] * b[1] - a[1] * b[3]
        q3 = a[0] * b[3] + a[3] * b[0] + a[1] * b[2] - a[2] * b[1]

        return np.array([q0, q1, q2, q3])

    @staticmethod
    def quat_inverse(q):
        """
        Інверсія кватерніону
        q: кватерніон
        return: кватерніон
        """
        return np.array([q[0], -q[1], -q[2], -q[3]])

    @staticmethod
    def euler_integration(f, state, time_span, time_step):
        """
        f: функція правої частини диференціального рівняння
        state: початкове значення стану
        time_span: проміжок часу
        time_step: крок інтегрування
        return: масив значень y на кожному кроці інтегрування
        """
        num_steps = int((time_span[1] - time_span[0]) / time_step) + 1
        for i in range(num_steps - 1):
            state_derivative = f(state)
            state += state_derivative * time_step

        return state
