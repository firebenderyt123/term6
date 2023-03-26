from geneticalgorithm import geneticalgorithm as ga
from bpla import start_bpla
import numpy as np

pos = np.array([
    0.0,
    50.0,
    0.0
])
v = np.array([
    0.0,
    0.0,
    0.0
])
psi = 0.0
teta = 0.0
gamma = 0.0
Om = np.array([0.0, 0.0, 0.0])  # вектор кутової швидкості МК
eps = np.array([
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
])
omega = np.array([
    21.0 * 2 * np.pi,
    10.0 * 2 * np.pi,
    10.0 * 2 * np.pi,
    10.0 * 2 * np.pi,
    10.0 * 2 * np.pi,
    17.0 * 2 * np.pi
])

def func(x):
    omeg = np.array([
        x[0],
        10.0 * 2 * np.pi,
        10.0 * 2 * np.pi,
        10.0 * 2 * np.pi,
        10.0 * 2 * np.pi,
        x[1]
    ])
    state = [
        pos, v,
        psi, teta, gamma,
        Om, omeg, eps
    ]
    states = start_bpla(0.04, state)
    x, y, z = states[-1][0:3]
    p, t, g = states[-1][6:9]
    Om1, Om2, Om3 = states[-1][9:12]
    # print(p, t, g)
    return abs(p) + abs(t) + abs(g)
    # return -y


def ai_work():
    init_val = func([omega[0], omega[5]])
    print(init_val)
    # exit()

    varbound = np.array([
        (20.0, 22.0),
        (16.0, 18.0)
    ]) * 2 * np.pi
    algorithm_param = {'max_num_iteration': 500,
                       'population_size': 50,
                       'mutation_probability': 0.01,
                       'elit_ratio': 0.01,
                       'crossover_probability': 0.95,
                       'parents_portion': 0.3,
                       'crossover_type': 'uniform',
                       'max_iteration_without_improv': None}

    # Створення об'єкту алгоритму
    model = ga(
        function=func,
        dimension=2,
        variable_type='real',
        variable_boundaries=varbound,
        algorithm_parameters=algorithm_param
    )

    # Запуск алгоритму
    model.run()

    # Виведення результатів
    print("Best solution:", model.output_dict['variable'])


if __name__ == "__main__":
    ai_work()
