import math
import json
import pandas as pd

file = "data.json"
out_file = "out.xlsx"

dt = 0.1  # крок інтегрування рівнянь
T = 800  # тривалість моделювання, cек
dmu_ = 1  # безрозмірний крок інтегрування чисельника та знаменника
mumin = -math.pi / 36  # рад./с2
mumax = math.pi / 36  # рад./с2

# Я задаю
fimin = -math.pi  # мінімальниц кут відхилення ШСЗ
fimax = math.pi  # максимальний кут відхилення ШСЗ
omegamin = -10 / 180 * math.pi  # мінімальна кутова швидкість
omegamax = 10 / 180 * math.pi  # максимальна кутова швидкість


df = pd.DataFrame(columns=["t", "fi", "omega", "mu"])  # save results here

with open(file, 'r', encoding="utf-8") as f:
    data = json.load(f)

# clear out file
with open(out_file, "w", encoding="utf-8") as f:
    f.write("")


def Ft(x, data):
    a, b, c, d = data
    if x <= a:
        return 0
    elif x >= a and x <= b:
        return (x - a) / (b - a)
    elif x >= b and x <= c:
        return 1
    elif x >= c and x <= d:
        return (d - x) / (d - c)
    elif d <= x:
        return 0
    else:
        print("Шось пішло не по плану.. x =", x)


# формула 2.1
def F(n, x):
    return Ft(x, data["F"][n-1])


def zfp(fi_, omega_, mu_):
    chistar = 0

    for n1, n2, m in data["data"]:
        alpha = F(n1, fi_)
        beta = F(n2, omega_)
        gamma = min(alpha, beta)
        delta = F(m, mu_)
        chi = min(gamma, delta)

        chistar = chi if (chistar < chi) else chistar

    return chistar


def main(fi, omega, t=0):
    iteration = 0
    while (t < T - 0.5 * dt):
        fi_ = 200 / (fimax - fimin) * (fi - fimin) - 100
        omega_ = 200 / (omegamax - omegamin) * (omega - omegamin) - 100

        mu_, s1, s2 = -100, 0, 0

        while mu_ < 100 - 0.5 * dmu_:
            chistar = zfp(fi_, omega_, mu_)
            s1 = s1 + mu_ * chistar * dmu_
            s2 = s2 + chistar * dmu_
            mu_ = mu_ + dmu_

        mu_star = s1 / s2
        mu = (mu_star + 100) / 200 * (mumax - mumin) + mumin

        fi += omega * dt
        omega += mu * dt
        t += dt

        # save to results array
        df.loc[iteration] = [t, fi, omega, mu]
        iteration += 1

        print(t, T)

    df.to_excel(out_file, index=False)


if __name__ == "__main__":
    fi = 60 / 180 * math.pi  # Початкове значення кут, рад
    omega = 5 / 180 * math.pi  # Початкове значення кутової швидкості, рад/c
    main(fi, omega)
