import numpy as np
import quaternion  # noqa
from scipy.integrate import quad
import sympy


class BPLA:

    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])
    e3 = np.array([0, 0, 1])

    def __init__(
        self, time,
        I, m, dt, J, J0, Kf, Kf0,  # noqa
        Km, Km0, p, G, Ks, H_, pos,
        Vx_, Vymax, Vz_, c,
        state
    ):
        self.time = time
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

        # Вектор локальної похідної кінетичного моменту двигунів
        self.dH = np.array([0.0, 0.0, 0.0])

        # Рівнодіючs силb тяги двигунів у проекціях на осі ІСК
        self.Fe = np.array([0.0, 0.0, 0.0])

        # Проекції сили опору повітря
        self.Fs = np.array([0.0, 0.0, 0.0])

        # Вектор гіроскопічного моменту в проекціях на осі ЗСК
        self.Mg_zsk = np.array([0.0, 0.0, 0.0])

        # Cумарний момент сил тяги двигунів
        self.MFe = np.array([0.0, 0.0, 0.0])

        # Cумарний аеродинамічний момент від шести двигунів
        self.MAe = np.array([0.0, 0.0, 0.0])

    def launch(self):
        np.seterr(all='raise')

        ts = self.time
        num_steps = int(ts / self.dt)
        states = []
        for i in range(num_steps):

            state = [
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
            ]

            # вивід поточного стану кожну секунду
            if (i % (num_steps / ts) == 0):
                self.t = i / (num_steps / ts)
                print(f'{self.t}/{ts}', state)

            states.append(state)

            self.calc_forces()
            # self.main_controller()

            self.derivative()

            # print(state[6:9])

        self.states = states

    # Функція, що повертає похідні від вектора стану
    def derivative(self):
        # розрахунок похідних
        self.x = self.x + self.vx * self.dt
        self.y = self.y + self.vy * self.dt
        self.z = self.z + self.vz * self.dt

        self.vx = self.vx + (
            self.Fe[0] / self.m + self.G[0] + self.Fs[0] / self.m
        ) * self.dt
        self.vy = self.vy + (
            self.Fe[1] / self.m + self.G[1] + self.Fs[1] / self.m
        ) * self.dt
        self.vz = self.vz + (
            self.Fe[2] / self.m + self.G[2] + self.Fs[2] / self.m
        ) * self.dt

        psi = (
            self.Om3 * np.sin(self.gamma) - self.Om2 * np.cos(self.gamma)
        ) / np.cos(self.teta)

        teta = self.Om2 * np.sin(
            self.gamma
        ) + self.Om3*np.cos(self.gamma)

        gamma = self.Om1 + np.tan(self.teta) * (
            self.Om3 * np.sin(self.gamma) - self.Om2 * np.cos(self.gamma)
        )

        self.psi = self.psi + psi * self.dt
        self.teta = self.teta + teta * self.dt
        self.gamma = self.gamma + gamma * self.dt

        # print(psi_dot, teta_dot, gamma_dot)

        # print(-self.Mg_zsk, self.dH, self.MFe, self.MAe)

        self.Om1 = self.Om1 + (
            (
                -self.Mg_zsk[0] - self.dH[0] + self.MFe[0] + self.MAe[0]
            ) / self.I
        ) * self.dt
        self.Om2 = self.Om2 + (
            (
                -self.Mg_zsk[1] - self.dH[1] + self.MFe[1] + self.MAe[1]
            ) / self.I
        ) * self.dt
        self.Om3 = self.Om3 + (
            (
                -self.Mg_zsk[2] - self.dH[2] + self.MFe[2] + self.MAe[2]
            ) / self.I
        ) * self.dt

        # print(Om_dot)

        self.omega0 = self.omega0 + self.eps0 * self.dt
        self.omega1 = self.omega1 + self.eps1 * self.dt
        self.omega2 = self.omega2 + self.eps2 * self.dt
        self.omega3 = self.omega3 + self.eps3 * self.dt
        self.omega4 = self.omega4 + self.eps4 * self.dt
        self.omega5 = self.omega5 + self.eps5 * self.dt

    # Управління  МК
    def main_controller(self):
        ay_ = self.vertical_controller()
        self.values_formation()
        # self.angle_controller()

        self.eps0 = ay_ * self.m / (2 * self.c[0] * self.omega0)
        # self.eps5 = * self.I[1] / (3 * self.c[5] * self.omega5)

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
        f = lambda t: self.v_zsk - self.Vz_  # noqa
        int_a = self.t  # ???
        int_b = self.t + self.dt  # ???

        self.gamma_ = 1 / (self.G[1] * self.m) * (
            self.m * (
                k[1] * (self.v_zsk - self.Vz_) + k[0] * np.array([
                    quad(lambda t: f(t)[i], int_a, int_b)[0] for i in range(3)
                ])
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

        x_, z_ = self.vx, self.vz

        d_R = np.quaternion(0, self.x - x_, 0, self.z - z_)  # ???
        d_R_zsk = (self.Lam.inverse() * d_R * self.Lam).vec

        d_vx, _, d_vz = d_R_zsk

        self.gamma_ = (
            self.m * (
                k[1] * self.vz + k[0] * d_vz
            ) - self.Fs[2] + abs(self.c[5]) * self.omega5 ** 2
        ) / (self.G[1] * self.m)

        self.teta_ = -(
            self.m * (
                k[1] * self.vx + k[0] * d_vx
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

        Fe_zsk = np.sum(self.F, axis=0)  # Сумарний вектор

        # кватерніон орієнтації ЗСК відносно ІСК
        print(self.psi)
        Lam_psi = np.quaternion(
            np.cos(self.psi/2),
            0,
            -np.sin(self.psi/2),
            0
        )
        Lam_teta = np.quaternion(
            np.cos(self.teta/2),
            0,
            0,
            np.sin(self.teta/2)
        )
        Lam_gamma = np.quaternion(
            np.cos(self.gamma/2),
            np.sin(self.gamma/2),
            0,
            0
        )
        self.Lam = Lam_psi * Lam_teta * Lam_gamma

        # 4

        Phiez = np.quaternion(
            0,
            Fe_zsk[0],
            Fe_zsk[1],
            Fe_zsk[2]
        )
        self.Fe = (self.Lam * Phiez * self.Lam.inverse()).vec

        # 5
        v = np.quaternion(0, self.vx, self.vy, self.vz)
        self.v_zsk = (self.Lam.inverse() * v * self.Lam).vec

        self.Fs_zsk = np.array([
            -self.Ks[0] * self.v_zsk[0] ** 2 * np.sign(self.v_zsk[0]),
            -self.Ks[1] * self.v_zsk[1] ** 2 * np.sign(self.v_zsk[1]),
            -self.Ks[2] * self.v_zsk[2] ** 2 * np.sign(self.v_zsk[2])
        ])
        Fs_zsk = np.quaternion(
            0,
            -self.Ks[0] * self.v_zsk[0] ** 2 * np.sign(self.v_zsk[0]),
            -self.Ks[1] * self.v_zsk[1] ** 2 * np.sign(self.v_zsk[1]),
            -self.Ks[2] * self.v_zsk[2] ** 2 * np.sign(self.v_zsk[2])
        )

        self.Fs = (self.Lam * Fs_zsk * self.Lam.inverse()).vec

        # 6

        # Вектор сумарного кінетичного моменту «МК + ротори з гвинтами»
        Om = np.array([self.Om1, self.Om2, self.Om3])
        H = np.array([
            self.I[0] * Om[0],
            self.I[1] * Om[1] + self.J * (
                self.omega1 - self.omega2 + self.omega3 - self.omega4
            ) + self.J0 * self.omega0,
            self.I[2] * Om[2] + self.J * self.omega5
        ])

        self.Mg_zsk = np.cross(Om, H)  # векторний доб

        # if self.Mg_zsk[0] != 0:
        #     self.Mg_zsk = self.Mg_zsk / np.linalg.norm(self.Mg_zsk)

        # print(self.Mg_zsk)
        # print(Om, H, self.Mg_zsk)

        # 7
        self.dH = np.array([
            0,
            self.J * (
                self.eps1 - self.eps2 + self.eps3 - self.eps4
            ) + self.J0 * self.eps0,
            self.J * self.eps5
        ])

        # 8
        self.MFe = np.sum(np.cross(self.p, self.F), axis=0)

        # 9

        # Аеродинамічні моменти від кожного гвинта
        Ma = np.array([
            self.Km0 * self.omega0 ** 2 * self.e2,
            self.Km * self.omega1 ** 2 * self.e2,
            -self.Km * self.omega2 ** 2 * self.e2,
            self.Km * self.omega3 ** 2 * self.e2,
            -self.Km * self.omega4 ** 2 * self.e2,
            self.Km * self.omega5 ** 2 * self.e3
        ])

        # 10
        self.MAe = np.sum(Ma, axis=0)
