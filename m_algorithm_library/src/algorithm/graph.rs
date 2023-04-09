use std::cmp::Ordering;
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
    //debug
    fn Debug(&self) {}
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

    heap.push(State {
        cost: T::default(),
        position: start,
    });

    while let Some(State { cost, position }) = heap.pop() {
        if dist[position].is_some() {
            continue;
        }

        dist[position] = Some(cost);

        for &(next_position, next_cost) in &graph.g[position] {
            if dist[next_position].is_none() {
                heap.push(State {
                    cost: cost + next_cost,
                    position: next_position,
                });
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
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dijkstra() {
        let mut graph = vec![
            vec![(1, 7), (2, 9), (5, 14)],
            vec![(0, 7), (2, 10), (3, 15)],
            vec![(0, 9), (1, 10), (3, 11), (5, 2)],
            vec![(1, 15), (2, 11), (4, 6)],
            vec![(3, 6), (5, 9)],
            vec![(0, 14), (2, 2), (4, 9)],
        ];
        let start = 0;
        let graph = Graph::from(graph);
        let dist = dijkstra::<i32>(start, &graph).0;
        let expected = vec![Some(0), Some(7), Some(9), Some(20), Some(20), Some(11)];
        assert_eq!(dist, expected);
    }

    #[test]
    fn test_dijkstra_path_restore() {
        let mut graph = vec![
            vec![(1, 7), (2, 9), (5, 14)],
            vec![(0, 7), (2, 10), (3, 15)],
            vec![(0, 9), (1, 10), (3, 11), (5, 2)],
            vec![(1, 15), (2, 11), (4, 6)],
            vec![(3, 6), (5, 9)],
            vec![(0, 14), (2, 2), (4, 9)],
        ];
        let graph = Graph::from(graph);
        let start = 0;
        let end = 4;
        let (dist, predecessor) = dijkstra::<i32>(start, &graph);
        let expected_dist = vec![Some(0), Some(7), Some(9), Some(20), Some(20), Some(11)];
        assert_eq!(dist, expected_dist);
        let path = reconstruct_path(&predecessor, start, end);
        let expected_path = vec![0, 1, 2, 5, 4];
        assert_eq!(path, expected_path);
    }
}
