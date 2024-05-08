import matplotlib.pyplot as plt
import numpy as np

colors = {
    'wigner_6j': 'green',
    'gsl_6j': 'red',
}

data = {key: [] for key in colors.keys()}

with open('data/bench_6j.txt', encoding='utf-8') as fp:
    key = None
    for line in fp:
        line = line.strip()
        if not line:
            continue
        if line[:-1] in colors.keys():
            key = line[:-1]
        else:
            data[key].append(line)

for key in colors.keys():
    data[key] = np.fromstring('\n'.join(data[key]), sep=' ').reshape(-1, 5)
    start = 20
    dj = data[key][start:, 0]/2
    mean_err = data[key][start:, 2]
    std_err = data[key][start:, 3]
    max_err = data[key][start:, 4]
    plt.yscale('log')
    plt.plot(dj, max_err, marker='x', color=colors[key], label=key+" max error")
    plt.plot(dj, std_err, marker = '.', color=colors[key], label=key+" std error")

plt.xlabel('Jmax')
plt.ylabel('Error')
plt.legend()
plt.savefig('data/bench_6j.svg')
