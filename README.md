# Optimized High Fidelity Sequence Search (OHFSS)
Реализация генетического алгоритма для управления кубитом с помощью биполярной/униполярной последовательности импульсов

# Запуск main.cpp
Для корректной работы программы необходимо указать следующие аргументы. Если аргумент не будет указан, то будет использоваться значение по умолчанию.

| Аргумент | Значение по умолчанию | Описание                                |
| ---------|:---------------------:| ---------------------------------------:|
| len      | 120  | Длина последовательности для генерации популяции |
| module   | 120  | Максимальная длина последовательности, должна быть кратна аргументу len |
| w01      | 3    | Коэффициент частоты кубита |
| w12      | 0.25 | Коэффициент w12 |
| wt      | 25 | Коэффициент частоты генератора |
| angle      | 0.024 | Желаемый угол последовательности |
| module      | 0.0001 | Максимальный модуль разности нужного угла и вычисленного |
| cp      | 0.8 | Вероятность скрещивания |
| mp      | 0.8 | Вероятность мутации |
| iter      | 500 | Количество итераций |
| type      | unipolar | Тип вычисляемой последовательности |

# Результат работы main.cpp
Результат работы программы дописывается в конец файла с названием "result_<аргументы запуска>". Если файла не существует, то он будет создан.
Формат вывода программы:

<Длина последовательности> \t <Последовательность> \t <Количество циклов> \t <Угол> \t <Утечка> \t <Время выполнения> \n

# processing.py - Обработка результатов main.cpp
Для корректной работы программы необходимо указать следующие аргументы:

-h, --help           Показать сообщение с подсказкой.

--file FILENAME      Путь до обрабатываемого файла.

--angle ANGLE        Необходимый угол. Если не указывать, то фильтрации по этому параметру не будет.

--module MODULE      Максимальный модуль разности между необходимым углом и вычисленным.

--fidelity FIDELITY  Верхняя граница утечки. Если не указывать, то фильтрации по этому параметру не будет.

Вывод осуществляется с помощью стандартного потока вывода.
