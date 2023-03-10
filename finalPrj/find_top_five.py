import os
import sys


def read_pageranks(fname):
    results = {}
    with open(fname) as f:
        content = f.readlines()
        for line in content:
            data = line.replace('\n', '')
            data = data.split(',')
            results[int(data[0])] = float(data[1])
    return results


def add_dicts(D1, D2):
    for key in D1.keys():
        if key in D2:
            D2[key] += D1[key]
        else:
            D2[key] = D1[key]
    return D2


def compile_top_five(path):
    D = {}
    for file in os.listdir(path):
        if len(file.split('_')) == 1:
            os.remove(os.path.join(path, file))
            continue
        graph_name = file.split('_')[0]
        results = read_pageranks(path + '/' + file)
        if graph_name not in D:
            D[graph_name] = (results, 1)
        else:
            results = add_dicts(results, D[graph_name][0])
            D[graph_name] = (results, D[graph_name][1] + 1)
    for graph_name in D.keys():
        for node in D[graph_name][0].keys():
            D[graph_name][0][node] /= D[graph_name][1]
    for graph_name in D.keys():
        D[graph_name] = [(node, rank) for node, rank in D[graph_name][0].items()]
        D[graph_name] = sorted(
            D[graph_name], key=lambda pair: pair[1], reverse=True)
        D[graph_name] = D[graph_name][:5]
    return D


def main():
    # path = sys.argv[1]
    path = '/home/ian/Repos/cpts411-fork/final-proj-ian/out2/runtimes/top_results'
    top = compile_top_five(path)
    for graph_name, five in top.items():
        print(graph_name)
        for node, rank in five:
            print('\t%d : %lf' % (node, rank))
    for graph_name, five in top.items():
        outfile = os.path.join(path, (graph_name+'.csv'))
        with open(outfile, 'w') as f:
            f.write('vertex,rank\n')
            for node, rank in five:
                f.write('%d,%lf\n' % (node, rank))


if __name__ == "__main__":
    main()
