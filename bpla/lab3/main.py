import pandas as pd


def main(h, v, h__, v__, dt, eps, m1, m2):
    t = 0
    data = []

    s1 = -0.2 + 0.1j
    s2 = -0.2 - 0.1j

    k0 = (s1 * s2).real
    k1 = -(s1 + s2).real

    while h__ > eps:
        u = -k0 * h__ - k1 * v__

        h_new = v
        v_new = u
        h__new = v__ + m1 * (h - h__)
        v__new = u + m2 * (h - h__)

        h = h + h_new * dt
        v = v + v_new * dt
        h__ = h__ + h__new * dt
        v__ = v__ + v__new * dt
        print(h, h__)

        t += dt

        data.append([t, h, v, h__, v__])

    df = pd.DataFrame(data, columns=['t', 'h', 'v', 'h^', 'v^'])
    df.to_excel('output.xlsx', index=False)


if __name__ == "__main__":
    dt = 0.001
    h = 100
    v = 0
    h__ = 110
    v__ = 0
    eps = 0.01

    m1 = 0.6
    m2 = 0.36

    main(h, v, h__, v__, dt, eps, m1, m2)
