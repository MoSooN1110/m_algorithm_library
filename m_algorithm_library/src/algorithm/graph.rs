use std::cmp::Ordering;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::ops::Add;

#[derive(Clone, Eq, PartialEq)]
struct Graph<T> {
    g: Vec<Vec<(usize, T)>>,
}

impl<T> From<Vec<Vec<(usize, T)>>> for Graph<T> {
    fn from(v: Vec<Vec<(usize, T)>>) -> Self {
        Self { g: v }
    }
}

impl<T> Graph<T>
where
    T: Clone + Ord + Add<Output = T> + Default + PartialEq + Eq + PartialOrd + Ord + Copy,
{
    fn new(n: usize) -> Self {
        let g = vec![vec![]; n];
        Self { g }
    }
    fn add_edge(&mut self, s: usize, t: usize, cost: T) {
        self.g[s].push((t, cost.clone()));
    }
    fn add_edge_undirected(&mut self, s: usize, t: usize, cost: T) {
        self.add_edge(s, t, cost.clone());
        self.add_edge(t, s, cost);
    }
    fn from_adjacency_matrix(adjacency_matrix: Vec<Vec<T>>) -> Self {
        let n = adjacency_matrix.len();
        let mut g = vec![vec![]; n];
        for i in 0..n {
            for j in 0..n {
                if adjacency_matrix[i][j] != T::default() {
                    g[i].push((j, adjacency_matrix[i][j].clone()));
                }
            }
        }
        Self { g }
    }

    fn out_adjacency_matrix(&self) -> Vec<Vec<T>> {
        let n = self.g.len();
        let mut adjacency_matrix = vec![vec![T::default(); n]; n];
        for i in 0..n {
            for &(j, cost) in &self.g[i] {
                adjacency_matrix[i][j] = cost.clone();
            }
        }
        adjacency_matrix
    }

    fn out_adjacency_list(&self) -> Vec<Vec<(usize, T)>> {
        self.g.clone()
    }

    fn from_adjacency_list(adjacency_list: Vec<Vec<(usize, T)>>) -> Self {
        Self { g: adjacency_list }
    }
}

#[derive(Copy, Clone, Eq, PartialEq)]
struct State<T> {
    cost: T,
    position: usize,
}

impl<T: Ord> Ord for State<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        other.cost.cmp(&self.cost)
    }
}

impl<T: Ord> PartialOrd for State<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn dijkstra<T>(start: usize, graph: &Graph<T>) -> (Vec<Option<T>>, Vec<Option<usize>>)
where
    T: Copy + Ord + Add<Output = T> + Default,
{
    let n = graph.g.len();
    let mut dist: Vec<Option<T>> = vec![None; n];
    let mut predecessor: Vec<Option<usize>> = vec![None; n];
    let mut heap = BinaryHeap::new();

    heap.push(Reverse(State {
        cost: T::default(),
        position: start,
    }));

    while let Some(Reverse(State { cost, position })) = heap.pop() {
        if dist[position] > Some(cost) {
            continue;
        }
        dist[position] = Some(cost);

        for &(next_position, next_cost) in &graph.g[position] {
            if dist[next_position].is_none() || dist[next_position] > Some(cost + next_cost) {
                heap.push(Reverse(State {
                    cost: cost + next_cost,
                    position: next_position,
                }));
                dist[next_position] = Some(cost + next_cost);
                predecessor[next_position] = Some(position);
            }
        }
    }

    (dist, predecessor)
}
fn reconstruct_path(predecessor: &Vec<Option<usize>>, start: usize, end: usize) -> Vec<usize> {
    let mut path = Vec::new();
    let mut current = end;

    while current != start {
        path.push(current);
        current = predecessor[current].unwrap();
    }

    path.push(start);
    path.reverse();

    path
}
