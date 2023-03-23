import pandas as pd
import numpy as np

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

path = 'out/states.xlsx'
df = pd.read_excel(path)


def getRow(lineNum):
    return df.iloc[lineNum]


def printDiff(diff):
    for title, series in diff.iteritems():
        if abs(series) > 1:
            print(title, series)


if __name__ == '__main__':
    row0 = getRow(826)
    row1 = getRow(827)
    row2 = getRow(828)
    row3 = getRow(829)

    printDiff(row1 - row0)
    print('\n')
    printDiff(row2 - row1)
    print('\n', 'Шось не то')
    printDiff(row3 - row2)

    psi = row0.psi + (
        row1.Om3 * np.sin(row1.gamma) - row1.Om2 * np.cos(row1.gamma)
    ) / np.cos(row1.teta) * 0.01

    print(psi)

    psi = row1.psi + (
        row2.Om3 * np.sin(row2.gamma) - row2.Om2 * np.cos(row2.gamma)
    ) / np.cos(row2.teta) * 0.01

    print(psi)
