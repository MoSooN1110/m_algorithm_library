pub trait Monoid {
    type T: Clone;
    fn id() -> Self::T;
    fn op(a: &Self::T, b: &Self::T) -> Self::T;
}

pub enum MONXOR {}
impl Monoid for MONXOR {
    type T = usize;
    fn id() -> Self::T {
        0
    }
    fn op(a: &Self::T, b: &Self::T) -> Self::T {
        *a ^ *b
    }
}

pub struct SegTree<M: Monoid> {
    n: usize,
    size: usize,
    log: usize,
    buf: Vec<M::T>,
}

impl<M: Monoid> SegTree<M> {
    #[allow(dead_code)]
    fn new(n: usize) -> Self {
        let log = (n as f64).log2().ceil() as usize;
        let size = 2 << log;
        let buf = vec![M::id(); size];
        Self { n, size, log, buf }
    }

    pub fn set(&mut self, mut p: usize, x: M::T) {
        assert!(p < self.n);
        p += self.size / 2;
        self.buf[p] = x;
        while p > 1 {
            p >>= 1;
            self.buf[p] = M::op(&self.buf[p * 2], &self.buf[p * 2 + 1]);
        }
    }

    pub fn from_vec(v: &[M::T]) -> Self {
        let n = v.len();
        let mut seg = Self::new(n);
        seg.buf.resize(2 * seg.size, M::id());
        // for (i, &x) in v.iter().enumerate() {
        //     seg.buf[i + seg.size / 2] = x;
        // }
        for i in 0..n {
            seg.buf[i + seg.size / 2] = v[i].clone();
        }
        for i in (1..(seg.size / 2)).rev() {
            seg.buf[i] = M::op(&seg.buf[i * 2], &seg.buf[i * 2 + 1]);
        }
        seg
    }

    pub fn querry(&self, mut l: usize, mut r: usize) -> M::T {
        assert!(l <= r && r <= self.n);
        let mut sml = M::id();
        let mut smr = M::id();
        l += self.size / 2;
        r += self.size / 2;

        while l < r {
            if l & 1 != 0 {
                sml = M::op(&sml, &self.buf[l]);
                l += 1;
            }
            if r & 1 != 0 {
                r -= 1;
                smr = M::op(&self.buf[r], &smr);
            }
            l >>= 1;
            r >>= 1;
        }
        M::op(&sml, &smr)
    }
}

pub enum MONADD {}
impl Monoid for MONADD {
    type T = usize;
    fn id() -> Self::T {
        0
    }
    fn op(a: &Self::T, b: &Self::T) -> Self::T {
        *a + *b
    }
}

pub enum MONADDI {}
impl Monoid for MONADDI {
    type T = i64;
    fn id() -> Self::T {
        0
    }
    fn op(a: &Self::T, b: &Self::T) -> Self::T {
        *a + *b
    }
}
// existing code

#[cfg(test)]
mod tests {
    use super::{SegTree, MONADD};

    #[test]
    fn test_segtree() {
        //https://onlinejudge.u-aizu.ac.jp/problems/DSL_2_B
        let mut seg = SegTree::<MONADD>::new(10);
        assert_eq!(seg.buf, vec![0; 16 * 2]);

        let data = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let mut seg2 = SegTree::<MONADD>::from_vec(&data);
        assert_eq!(seg2.querry(0, 10), 55);
        seg2.set(0, 5);
        assert_eq!(seg2.querry(0, 10), 59);
        assert_eq!(seg2.querry(1, 5), 14);
        assert_eq!(seg2.querry(3, 8), 30);
    }
}
