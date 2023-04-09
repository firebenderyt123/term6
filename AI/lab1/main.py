from genetic import GeneticAlgorithm
import matplotlib.pyplot as plt


def start():
    func = lambda x, y: x + y # noqa
    x = (-2.0, 4.0)
    y = (0.0, 4.0)
    q = 0.1  # точність

    N = 200  # кількість особин
    iters = 500  # кількість ітерацій
    pc = 0.95  # ймовірність схрещування [0.9, 1]
    pm = 0.01  # ймовірність мутації [0, 0.01]
    be = 0.01  # полезная хрень [0, 0.01]

    ga = GeneticAlgorithm(func, x, y, q, N, iters, pc, pm, be)
    ga.start()

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 8))

    axs[0].plot(ga.best_fitness)
    axs[0].set_title('Найкраще значення функції на кожній ітерації')
    axs[0].set_xlabel('Ітерація')
    axs[0].set_ylabel('Значення')

    axs[1].plot(ga.F_list)
    axs[1].set_title('Пристосованість функції на кожній ітерації')
    axs[1].set_xlabel('Ітерація')
    axs[1].set_ylabel('Пристосованість')

    print(ga.getAnswer())

    plt.show()


if __name__ == "__main__":
    start()
