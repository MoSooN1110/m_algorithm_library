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

// テストモジュール
#[cfg(test)]
mod tests {
    use super::SegmentSum;

    // エイリアスを作成
    type SegmentSumI64 = SegmentSum<i64>;

    #[test]
    fn test_segment_sum() {
        let input: Vec<i64> = vec![1, 2, 3, 4, 5];
        let mut segment_sum = SegmentSumI64::new(input);

        segment_sum.build();
        // dbg!(segment_sum.arr.clone());
        assert_eq!(segment_sum.qerry(0, 2 + 1), 6); // 1 + 2 + 3 = 6
        assert_eq!(segment_sum.qerry(1, 3 + 1), 9); // 2 + 3 + 4 = 9
        assert_eq!(segment_sum.qerry(2, 4 + 1), 12); // 3 + 4 + 5 = 12
        assert_eq!(segment_sum.qerry(0, 4 + 1), 15); // 1 + 2 + 3 + 4 + 5 = 15
    }
}
