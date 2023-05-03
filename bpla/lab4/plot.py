import matplotlib.pyplot as plt
import numpy as np


def draw_plots(states):
    position = [state[0:3] for state in states]
    pos_labels = ['x', 'y', 'z']

    v = [state[3:6] for state in states]
    v_labels = ['Vx', 'Vy', 'Vz']

    angle = [np.rad2deg(state[6:9]) for state in states]
    angle_labels = ['psi', 'teta', 'gamma']

    Om = [np.rad2deg(state[9:12]) for state in states]
    Om_labels = ['Om1', 'Om2', 'Om3']

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

    axs[0, 0].plot(position, label=pos_labels)
    axs[0, 0].set_title('Поведінка координат відносно ІСК під час польоту')
    axs[0, 0].set_xlabel('Час, мс')
    axs[0, 0].set_ylabel('Координата, м')
    axs[0, 0].legend()

    axs[0, 1].plot(v, label=v_labels)
    axs[0, 1].set_title('Поведінка проекцій швидкості відповідно ІСК')
    axs[0, 1].set_xlabel('Час, мс')
    axs[0, 1].set_ylabel('Проекції вектору швидкості, м/с')
    axs[0, 1].legend()

    axs[1, 0].plot(angle, label=angle_labels)
    axs[1, 0].set_title('Поведінка кутів відповідно ІСК')
    axs[1, 0].set_xlabel('Час, мс')
    axs[1, 0].set_ylabel('Кут тангажу та крену, град')
    axs[1, 0].legend()

    axs[1, 1].plot(Om, label=Om_labels)
    axs[1, 1].set_title('Поведінка кутових швидкостей відповідно ІСК')
    axs[1, 1].set_xlabel('Час, мс')
    axs[1, 1].set_ylabel('Кутова швидкість, град/с')
    axs[1, 1].legend()

    plt.show()

    # omega = [state[12:18] for state in states]
    # omega_labels = ['omega0', 'omega1', 'omega2', 'omega3', 'omega4', 'omega5']  # noqa

    # plt.plot(omega, label=omega_labels)
    # plt.title('Угли')
    # plt.xlabel('Час, с')
    # plt.ylabel('Кутова швидкість, рад/с')

    # plt.legend()
    # plt.show()
