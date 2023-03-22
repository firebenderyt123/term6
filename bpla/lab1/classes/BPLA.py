import numpy as np
from scipy.integrate import quad
import sympy


class BPLA:

    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])
    e3 = np.array([0, 0, 1])

    def __init__(
            self,
            I, m, dt, J, J0, Kf, Kf0,  # noqa
            Km, Km0, p, G, Ks, H_, pos,
            Vx_, Vymax, Vz_, c,
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
        self.H_ = H_
        self.pos = pos
        self.Vx_ = Vx_
        self.Vymax = Vymax
        self.Vz_ = Vz_
        self.c = c

        [
            self.x,
            self.y,
            self.z,
            self.vx,
            self.vy,
            self.vz,
            self.psi,
            self.teta,
            self.gamma,
            self.Om1,
            self.Om2,
            self.Om3,
            self.omega0,
            self.omega1,
            self.omega2,
            self.omega3,
            self.omega4,
            self.omega5,
            self.eps0,
            self.eps1,
            self.eps2,
            self.eps3,
            self.eps4,
            self.eps5
        ] = state
        self.state = state

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

        ts = 50 - 0
        num_steps = int(ts / self.dt)
        state = self.state
        states = [state]
        for i in range(num_steps):
            # вивід поточного стану кожну секунду
            if (i % (num_steps / ts) == 0):
                self.t = i / (num_steps / ts)
                states.append(state)
                print(f'{self.t}/{ts}', state)

                # self.main_controller()
                self.calc_forces()

            prev_state = state
            state = self.derivative(state)
            state = prev_state + state * self.dt
            # [
            #     self.x,
            #     self.y,
            #     self.z,
            #     self.vx,
            #     self.vy,
            #     self.vz,
            #     self.psi,
            #     self.teta,
            #     self.gamma,
            #     self.Om1,
            #     self.Om2,
            #     self.Om3,
            #     self.omega0,
            #     self.omega1,
            #     self.omega2,
            #     self.omega3,
            #     self.omega4,
            #     self.omega5,
            #     self.eps0,
            #     self.eps1,
            #     self.eps2,
            #     self.eps3,
            #     self.eps4,
            #     self.eps5
            # ] = state

        self.states = states

    # Функція, що повертає похідні від вектора стану
    def derivative(self, state):
        x, y, z, vx, vy, vz, psi, teta, gamma, Om1, Om2, Om3, omega0, omega1, omega2, omega3, omega4, omega5, eps0, eps1, eps2, eps3, eps4, eps5 = state  # noqa

        # розрахунок похідних
        x_dot = vx
        y_dot = vy
        z_dot = vz
        v_dot = self.Fe / self.m + self.G + self.Fs / self.m

        psi_dot = (
            Om3 * np.sin(gamma) - Om2 * np.cos(gamma)
        ) / np.cos(teta)

        teta_dot = Om2 * np.sin(
            gamma
        ) + Om3*np.cos(gamma)

        gamma_dot = Om1 + np.tan(teta) * (
            Om3 * np.sin(gamma) - Om2 * np.cos(gamma)
        )

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
            omega5_dot,
            eps0,
            eps1,
            eps2,
            eps3,
            eps4,
            eps5
        ])

    # Управління  МК
    def main_controller(self):
        ay_ = self.vertical_controller()
        # self.values_formation()
        # self.angle_controller()

        self.eps0 = ay_ * self.m / (2 * self.c[0] * self.omega0)
        # self.eps5 =

    # Регулятор  вертикального  каналу  системи управління
    def vertical_controller(self):
        vy_dot = self.Fe[1] / self.m + self.G[1] + self.Fs[1] / self.m

        if abs(self.H_ - self.y) > 20:
            b = 1
            k = np.array([
                -b ** 2,
                -2 * b
            ])
            ay_ = k[1] * self.vy + k[0] * (
                vy_dot - self.Vymax*np.sign(self.H_ - self.y)
            )
        else:
            b = 0.5
            k = np.array([
                -b ** 3,
                -3 * b ** 2,
                -3 * b
            ])
            ay_ = k[2] * vy_dot + k[1] * self.vy + k[0] * (self.y - self.H_)

        return ay_

    # Регулятор  каналу  курса  системи управління
    def angle_controller(self):
        b = 5
        k = np.array([
            -b ** 3,
            -3 * b ** 2,
            -3 * b
        ])

        f = (
            self.Om3 * np.sin(self.gamma) - self.Om2 * np.cos(self.gamma)
        ) / np.cos(self.teta)
        df = sympy.diff(f, x)

        # psi__  # ??????

        # self.sigma_y_ = k[2] * psi__

    # Формування  заданих  значень керованих  змінних,  виходячи  з  запланованих режимів  # noqa
    def values_formation(self):
        '''
        Формування  заданого  значення  кута  курсу
        для націлювання МК на задану точку призначення
        здійснюється за формулою
        '''
        dx, _, dz = self.pos - np.array([self.x, self.y, self.z])
        dx_plus_dz = (dx ** 2 + dz ** 2) ** 0.5

        if dx_plus_dz < 10:
            self.psi_ = self.psi
        else:
            self.psi_ = np.arccos(dx / dx_plus_dz) * np.sign(dz)

        '''
        Формування потрібного  значення  кута крену
        для  реалізації  заданої  поперечної  швидкості  МК
        здійснюється за формулою
        '''
        b = 0.5
        k = np.array([
            -b ** 2,
            -2 * b
        ])
        int_a = self.t  # ???
        int_b = self.t + self.dt  # ???
        self.gamma_ = 1 / (self.G[1] * self.m) * (
            self.m * (
                k[1] * (self.v_zsk - self.Vz_) + k[0] * quad(
                    self.v_zsk - self.Vz_, int_a, int_b  # ???
                )
            ) - self.Fs_zsk[2] + abs(self.c[5]) * self.omega5 ** 2
        )

        '''
        Формування  потрібного  значення  кута  тангажу
        для  реалізації  заданої  поздовжньої  швидкості
        МК здійснюється за формулою
        '''
        k1 = -0.1
        self.teta_ = -(
            self.m * k1 * (self.vx - self.Vx_) - self.Fs[0]
        ) / (self.G[1] * self.m)

        '''
         Формування змінних значень кутів  крену та
        тангажу для стабілізації МК у заданій точці
        здійснюється за формулами
        '''
        b = 0.1
        k = np.array([
            -b ** 2,
            -2 * b
        ])

        d_R = np.array([0, self.x - self.x_, 0, self.z - self.z_])  # ???
        d_R_zsk = BPLA.quat_mult(
            BPLA.quat_mult(self.Lam_inv, d_R),
            self.Lam
        )

        d_vx, _, d_vz = d_R_zsk

        self.gamma_ = (
            self.m * (
                self.k[1] * self.vz + k[0] * d_vz
            ) - self.Fs[2] + abs(self.c[5]) * self.omega5 ** 2
        ) / (self.G[1] * self.m)

        self.teta_ = -(
            self.m * (
                self.k[1] * self.vx + self.k[0] * d_vx
            ) - self.Fs[0]
        ) / (self.G[1] * self.m)

    def calc_forces(self):

        # 1
        # Вектори сили тяги двигунів
        self.F = np.array([
            self.Kf0 * self.omega0 ** 2 * self.e2,
            self.Kf * self.omega1 ** 2 * self.e2,
            self.Kf * self.omega2 ** 2 * self.e2,
            self.Kf * self.omega3 ** 2 * self.e2,
            self.Kf * self.omega4 ** 2 * self.e2,
            self.Kf * self.omega5 ** 2 * self.e3
        ])

        self.Fe_zsk = np.sum(self.F, axis=0)  # Сумарний вектор

        # кватерніон орієнтації ЗСК відносно ІСК
        self.Lam_psi = np.array([
            np.cos(self.psi/2),
            0,
            -np.sin(self.psi/2),
            0
        ])
        self.Lam_teta = np.array([
            np.cos(self.teta/2),
            0,
            0,
            np.sin(self.teta/2)
        ])
        self.Lam_gamma = np.array([
            np.cos(self.gamma/2),
            np.sin(self.gamma/2),
            0,
            0
        ])

        # Lam1 = np.multiply(Lam_psi, Lam_teta, Lam_gamma)

        self.Lam = BPLA.quat_mult(
            BPLA.quat_mult(self.Lam_psi, self.Lam_teta),
            self.Lam_gamma
        )

        # 4
        self.Phiez = np.hstack([0, self.Fe_zsk])
        self.Lam_inv = BPLA.quat_inverse(self.Lam)
        self.Fe = BPLA.vect(  # ???
            BPLA.quat_mult(
                BPLA.quat_mult(self.Lam, self.Phiez),
                self.Lam_inv
             )
        )

        # 5
        self.v_zsk = BPLA.quat_mult(
            BPLA.quat_mult(
                self.Lam_inv,
                np.array([0, self.vx, self.vy, self.vz])
            ),
            self.Lam
         )

        self.Fs_zsk = np.array([
            0,  # ???
            -self.Ks[0] * self.v_zsk[0] ** 2 * np.sign(self.v_zsk[0]),
            -self.Ks[1] * self.v_zsk[1] ** 2 * np.sign(self.v_zsk[1]),
            -self.Ks[2] * self.v_zsk[2] ** 2 * np.sign(self.v_zsk[2])
        ])
        self.Fs = BPLA.vect(  # ???
            BPLA.quat_mult(
                BPLA.quat_mult(self.Lam, self.Fs_zsk),
                self.Lam_inv
            )
        )

        # 6

        # Вектор сумарного кінетичного моменту «МК + ротори з гвинтами»
        self.H = np.array([
            self.I[0] * self.Om1,
            self.I[1] * self.Om2 + self.J * (
                self.omega1 - self.omega2 + self.omega3 - self.omega4
            ) + self.J0 * self.omega0,
            self.I[2] * self.Om3 + self.J * self.omega5
        ])

        self.Om = np.array([self.Om1, self.Om2, self.Om3])

        self.Mg = np.cross(self.Om, self.H)  # векторний доб

        # print(Om, H)

        # 7
        self.dH = np.array([
            0,
            self.J * (
                self.eps1 - self.eps2 + self.eps3 - self.eps4
            ) + self.J0 * self.eps0,
            self.J * self.eps5
        ])

        # 8 ???
        self.MFe = np.sum(np.cross(self.p, self.F), axis=0)

        # 9

        # Аеродинамічні моменти від кожного гвинта
        self.Ma = np.array([
            self.Km0 * self.omega0 ** 2 * self.e2,
            self.Km * self.omega1 ** 2 * self.e2,
            -self.Km * self.omega2 ** 2 * self.e2,
            self.Km * self.omega3 ** 2 * self.e2,
            -self.Km * self.omega4 ** 2 * self.e2,
            self.Km * self.omega5 ** 2 * self.e3
        ])

        # 10
        self.MAe = np.sum(self.Ma, axis=0)

    @staticmethod
    def vect(q):
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
