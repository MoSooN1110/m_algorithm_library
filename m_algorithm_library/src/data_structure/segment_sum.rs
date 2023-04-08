#![allow(dead_code)]
struct SegmentSum<T> {
    arr: Vec<T>,
}

use std::ops::{Add, Sub};

impl<T: Copy + Eq + Add<Output = T> + Sub<Output = T> + Default> SegmentSum<T> {
    fn new(arr: Vec<T>) -> Self {
        let mut carr = arr.clone();
        carr.insert(0, T::default());
        Self { arr: carr }
    }
    fn build(&mut self) {
        for i in 1..self.arr.len() {
            self.arr[i] = self.arr[i] + self.arr[i - 1];
        }
    }
    fn qerry(&self, l: usize, r: usize) -> T {
        // [l, r)
        assert!(l <= r && r < self.arr.len());
        self.arr[r] - self.arr[l]
    }
}

struct SegmentSum2D<T> {
    arr: Vec<Vec<T>>,
}
impl<T: Copy + Eq + Add<Output = T> + Sub<Output = T> + Default> SegmentSum2D<T> {
    fn new(arr: Vec<Vec<T>>) -> Self {
        let mut carr = arr.clone();
        carr.insert(0, vec![T::default(); carr[0].len()]);
        for i in 0..carr.len() {
            carr[i].insert(0, T::default());
        }
        Self { arr: carr }
    }
    fn build(&mut self) {
        for i in 1..self.arr.len() {
            for j in 1..self.arr[i].len() {
                self.arr[i][j] = self.arr[i][j] + self.arr[i - 1][j] + self.arr[i][j - 1]
                    - self.arr[i - 1][j - 1];
            }
        }
    }
    fn qerry(&self, i1: usize, i2: usize, j1: usize, j2: usize) -> T {
        // [i1, i2) * [j1, j2)
        assert!(i1 <= i2 && i2 < self.arr.len());
        assert!(j1 <= j2 && j2 < self.arr[0].len());
        self.arr[i2][j2] - self.arr[i1][j2] - self.arr[i2][j1] + self.arr[i1][j1]
    }
}

// テストモジュール
#[cfg(test)]
mod tests {
    use super::SegmentSum;

    // エイリアスを作成
    type SegmentSumI64 = SegmentSum<i64>;

    #[test]
    fn test_segment_sum() {
        let input: Vec<i64> = vec![1, 2, 3, 4, 5];
        let mut segment_sum = SegmentSumI64::new(input.clone());

        segment_sum.build();
        // dbg!(segment_sum.arr.clone());
        assert_eq!(segment_sum.qerry(0, 2 + 1), 6); // 1 + 2 + 3 = 6
        assert_eq!(segment_sum.qerry(1, 3 + 1), 9); // 2 + 3 + 4 = 9
        assert_eq!(segment_sum.qerry(2, 4 + 1), 12); // 3 + 4 + 5 = 12
        assert_eq!(segment_sum.qerry(0, 4 + 1), 15); // 1 + 2 + 3 + 4 + 5 = 15
    }

    // segment sum2d test
    use super::SegmentSum2D;
    #[test]
    fn test_segment_sum_2d() {
        let input: Vec<Vec<i64>> = vec![
            vec![1, 2, 3, 4, 5],
            vec![2, 3, 4, 5, 6],
            vec![3, 4, 5, 6, 7],
            vec![4, 5, 6, 7, 8],
            vec![5, 6, 7, 8, 9],
        ];
        let mut segment_sum = SegmentSum2D::new(input.clone());

        segment_sum.build();
        // dbg!(segment_sum.arr.clone());
        assert_eq!(segment_sum.qerry(0, 2 + 1, 0, 2 + 1), 27); // 1 + 2 + 3 + 2 + 3 + 4 + 3 + 4 + 5 = 30
        let sum = input.iter().flatten().sum::<i64>();
        assert_eq!(segment_sum.qerry(1, 4 + 1, 0, 4 + 1), sum - 15); //
        assert_eq!(segment_sum.qerry(0, 4 + 1, 0, 4 + 1), sum); //
    }
}
