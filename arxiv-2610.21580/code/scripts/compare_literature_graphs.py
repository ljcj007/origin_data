"""核对五格点模型与文献六格点模型的连接差异。"""
import argparse
from collections import Counter, deque
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GRAPH = ROOT / 'graphs/super_honeycomb_JA_JB.dat'


def current_edges():
    rows = GRAPH.read_text(encoding='utf-8').splitlines()
    sites, count, _ = map(int, rows[0].split())
    edges = []
    for row in rows[1:count+1]:
        a, b, dx, dy, dz, kind = map(int, row.split())
        if dz != 0:
            raise ValueError('本检验只支持二维连接')
        edges.append((a-1, b-1, dx, dy, kind))
    return sites, edges


def adhikary_edges():
    # arXiv:2102.12291v2，式 (1a)-(1c)，转为零起始编号。
    return 6, [(5, 0, 0, 0, 1), (1, 2, 0, 0, 1), (3, 4, 0, 0, 1),
               (0, 1, 0, 0, 2), (2, 3, 0, 0, 2), (4, 5, 0, 0, 2),
               (0, 3, 0, 1, 3), (2, 5, 1, 0, 3), (4, 1, -1, -1, 3)]


def invariants(sites, edges, length=3, active_types=None):
    if length < 3:
        raise ValueError('使用 L>=3 避免小环面重边混淆')
    graph = [set() for _ in range(sites*length**2)]
    def index(x, y, site):
        return ((x % length)*length + y % length)*sites + site
    for x in range(length):
        for y in range(length):
            for a, b, dx, dy, kind in edges:
                if active_types is not None and kind not in active_types:
                    continue
                i, j = index(x, y, a), index(x+dx, y+dy, b)
                if i == j or j in graph[i]:
                    raise ValueError('出现自环或重复无向键')
                graph[i].add(j)
                graph[j].add(i)
    colors, components, partitions = {}, [], []
    for start in range(len(graph)):
        if start in colors:
            continue
        queue, sizes = deque([start]), [1, 0]
        colors[start] = 0
        while queue:
            i = queue.popleft()
            for j in graph[i]:
                if j in colors:
                    if colors[j] == colors[i]:
                        raise ValueError('发现非二分连接')
                    continue
                colors[j] = 1-colors[i]
                sizes[colors[j]] += 1
                queue.append(j)
        components.append(sum(sizes))
        partitions.append(sorted(sizes))
    return dict(L=length, sites=len(graph), bonds=sum(map(len, graph))//2,
                degree_counts=dict(sorted(Counter(map(len, graph)).items())),
                component_size_counts=dict(sorted(Counter(components).items())),
                component_partitions=sorted(partitions), bipartite=True)


def report():
    current, other = current_edges(), adhikary_edges()
    return dict(scope='连接图检验，不是量子基态或热力学计算',
                source='https://arxiv.org/pdf/2102.12291', source_equations='1a-1c',
                provenance={str(GRAPH.relative_to(ROOT)): hashlib.sha256(GRAPH.read_bytes()).hexdigest(),
                            'scripts/compare_literature_graphs.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                full={str(l): dict(super_honeycomb=invariants(*current, length=l),
                                  adhikary_2021=invariants(*other, length=l)) for l in [3, 4]},
                strong_only=dict(super_honeycomb=invariants(*current, active_types={1}),
                                 adhikary_2021=invariants(*other, active_types={1})))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = report()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
    print(json.dumps(result, ensure_ascii=False))
