import matplotlib.pyplot as plt

# Задаем точки изотермы
C = [0.4, 1.2, 2.3, 2.7, 5.0]  # концентрация метилоранжу, мг/дм3
T = [36.6, 61.6, 85.0, 115, 116.6]  # величина адсорбции метилоранжу, мг/г

# Строим график
plt.plot(C, T, 'bo-')
plt.xlabel('Концентрація метилоранжу, мг/дм3')
plt.ylabel('Адсорбція метилоранжу, мг/г')
plt.title('Ізотерма адсорбції метилоранжу активованим вугіллям')
plt.show()
