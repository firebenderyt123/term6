import math
import random
from typing import Union, List, Tuple, Callable


class GeneticAlgorithm:

    def __init__(
        self,
        func: Callable,
        x_int: Tuple[float, float],
        y_int: Tuple[float, float],
        q: float,  # точність
        N: int,  # кількість  особин
        iters: int,  # кількість ітерацій
        pc: float = 0.95,  # ймовірність схрещування [0.9, 1]
        pm: float = 0.005,  # ймовірність мутації [0, 0.01]
        be: float = 0.005  # коефф для замены долбоебов на рандомных- [0, 0.01]
    ) -> None:
        self.func = func
        self.x_int = x_int
        self.y_int = y_int
        self.N = N
        self.iters = iters
        self.pc = pc
        self.pm = pm
        self.be = be

        # Кращий чел
        self.best_individ: Union[int, None] = None

        # Кращий результат
        self.max_fitness: Union[int, None] = None

        # список пристосованості на ітераціях
        self.F_list: List[float] = []

        # тут буду хранить лучшие значения функции на каждой итерации
        self.best_fitness: List[float] = []

        # тута будуть зберігатися x, y
        self.roots_list: List[List[float]] = []

        a, b = self.x_int
        c, d = self.y_int

        # шукається рішення з точністю q знаків після коми
        q = len(str(q).split('.')[1])

        # Найменше натуральне число L
        self.Lx = math.ceil(math.log2((b - a) * 10 ** q + 1))
        self.Ly = math.ceil(math.log2((d - c) * 10 ** q + 1))

        # N - число  підінтервалів розбиття
        Nx = 2 ** self.Lx - 1
        Ny = 2 ** self.Ly - 1

        # крок квантування інтервалу a, b
        self.hx = (b - a) / Nx
        self.hy = (d - c) / Ny

        # Початкові  значення  хромосом  особин

        # Find = []
        # ind = []
        # x_num = 4
        # y_num = 3
        # x_step = (b - a) / x_num
        # y_step = (d - c) / y_num
        # for i in range(1, x_num):
        #     for j in range(y_num-1, 0, -1):
        #         x = a + i * x_step
        #         y = c + j * y_step
        #         ind.append([x, y])
        #         Find.append(self.F(x, y))

        Find = []
        ind = []
        for i in range(1, Nx):
            for j in range(Ny-1, 0, -1):
                x = a + i * self.hx
                y = c + j * self.hy
                ind.append([x, y])
                Find.append(self.F(x, y))

        self.ind = ind
        self.Find = Find

    # fitness функція
    def F(self, x: float, y: float):
        return self.func(x, y)

    def start(self):

        self.presence()

        for i in range(self.iters):
            # print(i, self.Find)
            self.selection()

            ind = []
            for pair in self.parents:
                children = self.crossover(pair[0], pair[1])
                for child in children:
                    chromosome1 = self.mutate(child[0])
                    chromosome2 = self.mutate(child[1])
                    ind.append([chromosome1, chromosome2])

            self.ch_ind = ind

            self.newPresence()
            self.addBestFitness()

            # print(i, self.ch_ind)
            print(i, self.Y, self.F_, self.best_fitness[-1])

    # Обчислимо пристосованість усієї популяції
    def presence(self) -> None:
        self.Y = sum(y for y in self.Find)
        self.F_ = self.Y / self.N

        a, b = self.x_int
        c, d = self.y_int

        # Обчислюємо значення хромосом особин
        ind = []

        for x, y in self.ind:
            x_ = math.ceil((x - a) / self.hx)
            y_ = math.ceil((y - c) / self.hy)
            ch_x = bin(x_)[2:].rjust(self.Lx, '0')
            ch_y = bin(y_)[2:].rjust(self.Ly, '0')

            ind.append([ch_y, ch_x])

        self.ch_ind = ind

        # print(self.ch_ind)

    # Виконаємо селекцію за допомогою рулетки
    def selection(self) -> None:
        Fmin = min(self.Find)

        F_ind = [val + 2 * abs(Fmin) for val in self.Find]

        self.Y_ = sum(F_ind)
        self.F__ = self.Y_ / self.N

        # print(Fmin, Y_, F__)

        Sind = [val / self.Y_ * 100 for val in F_ind]
        S = []
        x = 0
        for y in Sind:
            S.append([x, x + y])
            x += y

        # Генеруємо випадковим чином 6 дійсних чисел з діапазону [0;100]
        nums = [random.uniform(0, 100) for i in range(self.N)]

        parents = []
        for num in nums:
            i = 0
            for a, b in S:
                if num >= a and num <= b:
                    parents.append(self.ch_ind[i])
                    break
                i += 1

        # розіб'ємо особин по парах
        random.shuffle(parents)
        self.parents = [
            [parents[i], parents[i+1]]
            for i in range(0, len(parents), 2)
        ]

    # Схрещування
    def crossover(
        self,
        parent1: str,
        parent2: str
    ) -> Tuple[str, str]:
        L = len(parent1)

        if random.random() <= self.pc:
            K = random.randint(1, L - 1)
            child1 = parent1[:K] + parent2[K:]
            child2 = parent2[:K] + parent1[K:]
            return child1, child2
        else:
            return parent1, parent2

    # Мутація
    def mutate(
        self,
        chromosome: str
    ) -> str:

        mutated_chromosome = ''
        for gen in chromosome:
            if random.random() < self.pm:
                mutated_chromosome += '0' if gen == '1' else '1'
            else:
                mutated_chromosome = chromosome
                break
        return mutated_chromosome

    # Обчислимо пристосованість нових особин, які формують нову популяцію
    def newPresence(self):
        a, b = self.x_int
        c, d = self.y_int

        Find = []
        ind = []
        best_individ = self.best_individ
        max_fitness = self.max_fitness
        for i in self.ch_ind:
            y_ = int(i[0], 2)
            x_ = int(i[1], 2)

            y = c + y_ * self.hy
            x = a + x_ * self.hx

            individ = [x, y]
            fitness = self.F(x, y)

            ind.append(individ)
            Find.append(fitness)

            if max_fitness is None or max_fitness < fitness:
                best_individ = individ
                max_fitness = fitness

        self.best_individ = best_individ
        self.max_fitness = max_fitness

        # Модификация
        zipped = list(zip(ind, Find, self.ch_ind))
        zipped_sorted = sorted(zipped, key=lambda x: x[1], reverse=True)

        ind_s, Find_s, ch_ind_s = zip(*zipped_sorted)
        elite = list(ind_s[:self.N - math.floor(self.N * self.be)])
        elite_F = list(Find_s[:self.N - math.floor(self.N * self.be)])
        elite_chind = list(ch_ind_s[:self.N - math.floor(self.N * self.be)])

        new_ind = elite
        new_chind = elite_chind
        new_Find = elite_F
        for i in range(math.floor(self.N * self.be)):
            x = random.uniform(0, b)
            y = random.uniform(0, d)

            x_ = math.ceil((x - a) / self.hx)
            y_ = math.ceil((y - c) / self.hy)

            ch_x = bin(x_)[2:].rjust(self.Lx, '0')
            ch_y = bin(y_)[2:].rjust(self.Ly, '0')

            new_ind.append([x, y])
            new_chind.append([ch_y, ch_x])
            new_Find.append(self.F(x, y))

        self.ind = new_ind
        self.ch_ind = new_chind
        self.Find = new_Find

        # self.ind = ind
        # self.Find = Find

        self.Y = sum(Find)
        self.F_ = self.Y / self.N

    # Добавляем лучшее значение функции в массив
    def addBestFitness(self):
        self.F_list.append(self.F_)
        self.best_fitness.append(self.max_fitness)
        self.roots_list.append(self.best_individ)

    # отримати найкращє значення функції
    def getAnswer(self):
        return {
            'best_fitness': max(self.best_fitness),
            'best_roots': max(self.roots_list)
        }


if __name__ == "__main__":

    func = lambda x, y: x ** 2 - y ** 2  # noqa
    x = (0, 4)
    y = (1, 2)
    q = 0.1  # точність

    N = 6  # кількість особин
    iters = 2  # кількість ітерацій

    ga = GeneticAlgorithm(func, x, y, q, N, iters)
    ga.start()
