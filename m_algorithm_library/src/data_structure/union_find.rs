use std::collections::*;
pub struct UnionFind {
    n: usize,
    // root node: -1 * component size
    // otherwise: parent
    parent_or_size: Vec<i32>,
}

impl UnionFind {
    // 0 <= size <= 10^8 is constrained.
    pub fn new(size: usize) -> Self {
        Self {
            n: size,
            parent_or_size: vec![-1; size],
        }
    }
    pub fn merge(&mut self, a: usize, b: usize) -> usize {
        assert!(a < self.n);
        assert!(b < self.n);
        let (mut x, mut y) = (self.leader(a), self.leader(b));
        if x == y {
            return x;
        }
        if -self.parent_or_size[x] < -self.parent_or_size[y] {
            std::mem::swap(&mut x, &mut y);
        }
        self.parent_or_size[x] += self.parent_or_size[y];
        self.parent_or_size[y] = x as i32;
        x
    }

    pub fn same(&mut self, a: usize, b: usize) -> bool {
        assert!(a < self.n);
        assert!(b < self.n);
        self.leader(a) == self.leader(b)
    }
    pub fn leader(&mut self, a: usize) -> usize {
        assert!(a < self.n);
        if self.parent_or_size[a] < 0 {
            return a;
        }
        self.parent_or_size[a] = self.leader(self.parent_or_size[a] as usize) as i32;
        self.parent_or_size[a] as usize
    }
    pub fn size(&mut self, a: usize) -> usize {
        assert!(a < self.n);
        let x = self.leader(a);
        -self.parent_or_size[x] as usize
    }
    pub fn groups(&mut self) -> Vec<Vec<usize>> {
        let mut leader_buf = vec![0; self.n];
        let mut group_size = vec![0; self.n];
        for i in 0..self.n {
            leader_buf[i] = self.leader(i);
            group_size[leader_buf[i]] += 1;
        }
        let mut result = vec![Vec::new(); self.n];
        for i in 0..self.n {
            result[i].reserve(group_size[i]);
        }
        for i in 0..self.n {
            result[leader_buf[i]].push(i);
        }
        result
            .into_iter()
            .filter(|x| !x.is_empty())
            .collect::<Vec<Vec<usize>>>()
    }

    pub fn connected_component_num(&mut self) -> usize {
        let mut set = HashSet::new();
        for i in 0..self.n {
            set.insert(self.leader(i));
        }
        set.len()
    }

    pub fn leader_vec(&mut self) -> Vec<usize> {
        let mut leader_buf = vec![0; self.n];
        for i in 0..self.n {
            leader_buf[i] = self.leader(i);
        }
        leader_buf
    }
}
#[cfg(test)]
mod tests {
    use super::UnionFind;

    #[test]
    fn test_new() {
        let uf = UnionFind::new(5);
        assert_eq!(uf.n, 5);
        assert_eq!(uf.parent_or_size, vec![-1, -1, -1, -1, -1]);
    }

    #[test]
    fn test_merge() {
        let mut uf = UnionFind::new(5);
        uf.merge(1, 2);
        assert!(uf.same(1, 2));
    }

    #[test]
    fn test_same() {
        let mut uf = UnionFind::new(5);
        uf.merge(1, 2);
        assert!(uf.same(1, 2));
        assert!(!uf.same(1, 3));
    }

    #[test]
    fn test_leader() {
        let mut uf = UnionFind::new(5);
        uf.merge(1, 2);
        uf.merge(3, 4);
        assert_eq!(uf.leader(1), 1);
        assert_eq!(uf.leader(2), 1);
        assert_eq!(uf.leader(3), 3);
        assert_eq!(uf.leader(4), 3);
    }

    #[test]
    fn test_size() {
        let mut uf = UnionFind::new(5);
        uf.merge(1, 2);
        uf.merge(3, 4);
        assert_eq!(uf.size(1), 2);
        assert_eq!(uf.size(2), 2);
        assert_eq!(uf.size(3), 2);
        assert_eq!(uf.size(4), 2);
    }

    #[test]
    fn test_groups() {
        let mut uf = UnionFind::new(5);
        uf.merge(1, 2);
        uf.merge(3, 4);
        let groups = uf.groups();
        assert_eq!(groups.len(), 3);
        assert!(groups.contains(&vec![0]));
        assert!(groups.contains(&vec![1, 2]));
        assert!(groups.contains(&vec![3, 4]));
    }

    #[test]
    fn test_connected_component_num() {
        let mut uf = UnionFind::new(5);
        uf.merge(1, 2);
        uf.merge(3, 4);
        assert_eq!(uf.connected_component_num(), 3);
    }

    #[test]
    fn test_leader_vec() {
        let mut uf = UnionFind::new(5);
        uf.merge(1, 2);
        uf.merge(3, 4);
        let leaders = uf.leader_vec();
        assert_eq!(leaders, vec![0, 1, 1, 3, 3]);
    }
}
