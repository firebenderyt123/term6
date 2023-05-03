import pandas as pd
import math


def main(omega, fi, omega__, fi__, dt, eps, m1, m2):
    t = 0
    data = []

    s1 = -0.2 + 0.1j
    s2 = -0.2 - 0.1j

    k1 = (s1 * s2).real
    k2 = -(s1 + s2).real

    while omega__ > eps:
        u = -k1 * omega__ - k2 * fi__

        omega_new = fi
        fi_new = u
        omega__new = fi__ + m1 * (omega - omega__)
        fi__new = u + m2 * (omega - omega__)

        omega = omega + omega_new * dt
        fi = fi + fi_new * dt
        omega__ = omega__ + omega__new * dt
        fi__ = fi__ + fi__new * dt
        print(omega, omega__)

        t += dt

        data.append([t, omega, fi, omega__, fi__])

    df = pd.DataFrame(data, columns=['t', 'omega', 'fi', 'omega^', 'fi^'])
    df.to_excel('output.xlsx', index=False)


if __name__ == "__main__":
    dt = 0.1
    fi = 60 / 180 * math.pi  # Початкове значення кут, рад
    omega = 5 / 180 * math.pi  # Початкове значення кутової швидкості, рад/c
    omega__ = 5 / 180 * math.pi
    fi__ = 60 / 180 * math.pi
    eps = 0.01

    # коефіцієнти спостерігаючого пристрою
    m1 = 0.6
    m2 = 0.36

    main(omega, fi, omega__, fi__, dt, eps, m1, m2)
