import numpy as np
import matplotlib.pyplot as plt

file_in = 'data4app.txt'
file_out = 'out.txt'

E = np.eye(4)
print("E =", E)


def calc():
    # Ініціалізація
    F_ = E
    delta_ = np.array([
        [1],
        [1],
        [1],
        [1]
    ])
    # print(delta_, delta_.T)

    S = np.array([1, 2, 3, 4])
    eps = 1e-06
    N = 4

    data = np.loadtxt(file_in)
    L = len(data)
    P_list = np.zeros((L, 4))
    taus = np.zeros((L, 1))
    df_list = np.zeros(L)
    f_list = np.zeros((L, 4))

    for i, (df, T) in enumerate(data):
        tau = T / 50
        taus[i] = tau
        df_list[i] = df  # похибка з файлу
        f = np.array([
            [1],
            [tau],
            [tau**2],
            [tau**3]
        ])
        f_list[i] = f.T

        # print(F_ * f, f.T * F_)
        # print(F_ * f * f.T * F_ == (F_ * f) * (f.T * F_))
        # exit()

        gamma = 1 + np.dot(np.dot(f.T, F_), f)
        F_ = F_ - 1 / gamma * np.dot(F_, f) * np.dot(f.T, F_)
        delta_ = delta_ + np.dot(f, df)

        if N > 0:
            for s in S:
                e_s = E[:, int(s-1)]
                e_s = np.reshape(e_s, (len(e_s), 1))
                g_prev = gamma  # для вивіда попередньої gamma
                gamma = 1 - np.dot(np.dot(e_s.T, F_), e_s)
                if abs(gamma) > eps:
                    print('Итерация, gamma', i, g_prev)
                    print('Итерация, F matrix', i, F_)
                    F_ = F_ + 1 / gamma * np.dot(F_, e_s) * np.dot(e_s.T, F_)
                    print('Итерация, gamma', i, gamma)
                    print('Итерация, F matrix', i, F_)
                    delta_ = delta_ - e_s
                    N = N - 1
                    S = [z for z in S if z != s]
                    print("Итерация, S:", i, S)
                    # exit()
                    break

        P = np.dot(F_, delta_)
        P_list[i] = P.T

        # exit()

    dT = P[0] + P[1] * taus + P[2] * taus**2 + P[3] * taus**3

    da = 1 / L * np.sum(df_list - dT)
    ds = np.sqrt(1 / L * np.sum((df_list - dT) - da)**2)

    print(P, da, ds)

    np.savetxt(file_out, P_list)

    return L, df_list, dT


def draw_plot(time, df_list, dT):
    plt.plot(time, df_list)
    plt.xlabel('Час, с')
    plt.ylabel('Похибка, град/год')
    plt.grid()

    plt.plot(time, dT)
    plt.xlabel('Час, с')
    plt.ylabel('Модельне значення')
    plt.grid()
    plt.show()


def start():
    N, df_list, dT = calc()
    draw_plot(np.arange(N), df_list, dT)


if __name__ == "__main__":
    start()
