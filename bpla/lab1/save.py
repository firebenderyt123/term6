import os
import pandas as pd


def save_array(data, titles=[], path="out/states.xlsx"):

    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))

    if len(titles) != len(data):
        titles = [
            'X', 'Y', 'Z',
            'Vx', 'Vy', 'Vz',
            'psi', 'teta', 'gamma',
            'Om1', 'Om2', 'Om3',
            'omega0',
            'omega1',
            'omega2',
            'omega3',
            'omega4',
            'omega5',
            'eps0',
            'eps1',
            'eps2',
            'eps3',
            'eps4',
            'eps5'
        ]

    df = pd.DataFrame(data, columns=titles)

    df.to_excel(path, index=False)
