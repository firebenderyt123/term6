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
        Km, Km0, p, G, Ks, H_, end_pos,
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
        self.end_pos = end_pos
        self.Vx_ = Vx_
        self.Vymax = Vymax
        self.Vz_ = Vz_
        self.c = c

        [
            self.pos,
            self.v,
            self.psi,
            self.teta,
            self.gamma,
            self.Om,
            self.omega,
            self.eps
        ] = state
        self.current_state = []

        self.F = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0]
        ])
        self.Lam = np.quaternion(0.0, 0.0, 0.0, 0.0)
        self.v_zsk = np.array([0.0, 0.0, 0.0])
        self.Fs_zsk = np.array([0.0, 0.0, 0.0])

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

    # повертає поточний стан
    def getCurrentState(self):
        return [
            self.pos[0],
            self.pos[1],
            self.pos[2],
            self.v[0],
            self.v[1],
            self.v[2],
            self.psi,
            self.teta,
            self.gamma,
            self.Om[0],
            self.Om[1],
            self.Om[2],
            self.omega[0],
            self.omega[1],
            self.omega[2],
            self.omega[3],
            self.omega[4],
            self.omega[5],
            self.eps[0],
            self.eps[1],
            self.eps[2],
            self.eps[3],
            self.eps[4],
            self.eps[5],
            self.F[0][0],
            self.F[0][1],
            self.F[0][2],
            self.F[1][0],
            self.F[1][1],
            self.F[1][2],
            self.F[2][0],
            self.F[2][1],
            self.F[2][2],
            self.F[3][0],
            self.F[3][1],
            self.F[3][2],
            self.F[4][0],
            self.F[4][1],
            self.F[4][2],
            self.F[5][0],
            self.F[5][1],
            self.F[5][2],
            self.Lam.w,
            self.Lam.x,
            self.Lam.y,
            self.Lam.z,
            self.Fe[0],
            self.Fe[1],
            self.Fe[2],
            self.v_zsk[0],
            self.v_zsk[1],
            self.v_zsk[2],
            self.Fs_zsk[0],
            self.Fs_zsk[1],
            self.Fs_zsk[2],
            self.Fs[0],
            self.Fs[1],
            self.Fs[2],
            self.Mg_zsk[0],
            self.Mg_zsk[1],
            self.Mg_zsk[2],
            self.dH[0],
            self.dH[1],
            self.dH[2],
            self.MFe[0],
            self.MFe[1],
            self.MFe[2],
            self.MAe[0],
            self.MAe[1],
            self.MAe[2]
        ]

    def launch(self):
        np.seterr(all='raise')

        states = []

        ts = self.time
        num_steps = int(ts / self.dt)
        for i in range(num_steps):

            state = self.getCurrentState()

            # вивід поточного стану кожну секунду
            if (i % (num_steps / ts) == 0):
                self.t = i / (num_steps / ts)
                # print(f'{self.t}/{ts}', state)

            self.calc_forces()
            # self.main_controller()

            self.derivative()

            states.append(state)

        self.states = states

    # Функція, що повертає похідні від вектора стану
    def derivative(self):
        # розрахунок похідних
        pos = self.v

        v = self.Fe / self.m + self.G + self.Fs / self.m

        psi = (
            self.Om[2] * np.sin(self.gamma) - self.Om[1] * np.cos(self.gamma)
        ) / np.cos(self.teta)

        teta = self.Om[1] * np.sin(
            self.gamma
        ) + self.Om[2] * np.cos(self.gamma)

        gamma = self.Om[0] + np.tan(self.teta) * (
            self.Om[2] * np.sin(self.gamma) - self.Om[1] * np.cos(self.gamma)
        )

        Om = (-self.Mg_zsk - self.dH + self.MFe + self.MAe) / self.I

        omega = self.eps

        self.pos = self.pos + pos * self.dt
        self.v = self.v + v * self.dt
        self.psi = self.psi + psi * self.dt
        self.teta = self.teta + teta * self.dt
        self.gamma = self.gamma + gamma * self.dt
        self.Om = self.Om + Om * self.dt
        self.omega = self.omega + omega * self.dt

    # Управління  МК
    def main_controller(self):
        self.vertical_controller()
        self.values_formation()
        # self.angle_controller()

        self.eps0 = self.ay * self.m / (2 * self.c[0] * self.omega0)
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
                vy_dot - self.Vymax * np.sign(self.H_ - self.y)
            )
        else:
            b = 0.5
            k = np.array([
                -b ** 3,
                -3 * b ** 2,
                -3 * b
            ])
            ay_ = k[2] * vy_dot + k[1] * self.vy + k[0] * (self.y - self.H_)

        self.ay = ay_

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
        dx, _, dz = self.end_pos - self.pos
        dx_plus_dz = (dx ** 2 + dz ** 2) ** 0.5

        if dx_plus_dz < 10:
            psi_ = self.psi
        else:
            psi_ = np.arccos(dx / dx_plus_dz) * np.sign(dz)

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
        F = np.array([
            self.Kf0 * self.omega[0] ** 2 * self.e2,
            self.Kf * self.omega[1] ** 2 * self.e2,
            self.Kf * self.omega[2] ** 2 * self.e2,
            self.Kf * self.omega[3] ** 2 * self.e2,
            self.Kf * self.omega[4] ** 2 * self.e2,
            self.Kf * self.omega[5] ** 2 * self.e3
        ])

        Fe_zsk = np.sum(F, axis=0)  # Сумарний вектор

        # кватерніон орієнтації ЗСК відносно ІСК
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
        v = np.quaternion(0, self.v[0], self.v[1], self.v[2])
        v_zsk = (Lam.inverse() * v * Lam).vec

        Fs_zsk = np.array([
            -self.Ks[0] * v_zsk[0] ** 2 * np.sign(v_zsk[0]),
            -self.Ks[1] * v_zsk[1] ** 2 * np.sign(v_zsk[1]),
            -self.Ks[2] * v_zsk[2] ** 2 * np.sign(v_zsk[2])
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
            self.I[0] * self.Om[0],
            self.I[1] * self.Om[1] + self.J * (
                self.omega[1] - self.omega[2] + self.omega[3] - self.omega[4]
            ) + self.J0 * self.omega[0],
            self.I[2] * self.Om[2] + self.J * self.omega[5]
        ])

        Mg_zsk = np.cross(self.Om, H)  # векторний доб

        # if self.Mg_zsk[0] != 0:
        #     self.Mg_zsk = self.Mg_zsk / np.linalg.norm(self.Mg_zsk)

        # print(self.Mg_zsk)
        # print(Om, H, self.Mg_zsk)

        # 7
        dH = np.array([
            0,
            self.J * (
                self.eps[1] - self.eps[2] + self.eps[3] - self.eps[4]
            ) + self.J0 * self.eps[0],
            self.J * self.eps[5]
        ])

        # 8
        MFe = np.sum(np.cross(self.p, F), axis=0)

        # 9

        # Аеродинамічні моменти від кожного гвинта
        Ma = np.array([
            self.Km0 * self.omega[0] ** 2 * self.e2,
            self.Km * self.omega[1] ** 2 * self.e2,
            -self.Km * self.omega[2] ** 2 * self.e2,
            self.Km * self.omega[3] ** 2 * self.e2,
            -self.Km * self.omega[4] ** 2 * self.e2,
            self.Km * self.omega[5] ** 2 * self.e3
        ])

        # 10
        MAe = np.sum(Ma, axis=0)

        self.F = F
        self.Lam = Lam
        self.Fe = Fe
        self.v_zsk = v_zsk
        self.Fs_zsk = Fs_zsk
        self.Fs = Fs
        self.Mg_zsk = Mg_zsk
        self.dH = dH
        self.MFe = MFe
        self.MAe = MAe
