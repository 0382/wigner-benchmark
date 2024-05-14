import matplotlib.pyplot as plt
import numpy as np

def load_data(fname, keys):
    data = {key: [] for key in keys}
    with open(fname, encoding='utf-8') as fp:
        key = None
        for line in fp:
            line = line.strip()
            if not line:
                continue
            if line[:-1] in keys:
                key = line[:-1]
            else:
                data[key].append(line)
    for key in keys:
        data[key] = np.fromstring('\n'.join(data[key]), sep=' ').reshape(-1, 6)
    return data

def plot_data(fname, data, colors, start=0):
    fig1 = plt.figure()
    fig2 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    for key in data.keys():
        j = data[key][start:, 0]/2
        max_err = data[key][start:, 2]
        max_rel_err = data[key][start:, 3]
        ax1.plot(j, max_err, marker='.', color=colors[key], label=key+" max error")
        ax2.plot(j, max_rel_err, marker = '.', color=colors[key], label=key+" max relative error")
    
    ax1.set_xlabel('Jmax')
    ax1.set_ylabel('Error')
    ax1.legend()
    fig1.savefig(fname+'_err.svg')
    ax2.set_xlabel('Jmax')
    ax2.set_ylabel('Relative Error')
    ax2.legend()
    fig2.savefig(fname+'_rel_err.svg')


if __name__ == '__main__':
    colors = {'wigner_3j': 'green', 'gsl_3j': 'red'}
    data = load_data('data/bench_3j.txt', colors.keys())
    plot_data('data/bench_3j', data, colors)
    colors = {'wigner_6j': 'green', 'gsl_6j': 'red'}
    data = load_data('data/bench_6j.txt', colors.keys())
    plot_data('data/bench_6j', data, colors)