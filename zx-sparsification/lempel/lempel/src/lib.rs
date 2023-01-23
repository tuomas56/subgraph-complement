use std::ops;

pub trait Vector: ops::BitXor<Self, Output=Self> + ops::Not<Output=Self> + PartialEq + Sized + Clone + Copy  {
    const WIDTH: usize;
    const ZERO: Self;

    fn from_bits(bits: impl Iterator<Item=bool>) -> Self;
    fn extract(&self, idx: usize) -> bool;
    fn set(&mut self, idx: usize, val: bool);
}

impl Vector for bool {
    const WIDTH: usize = 1;
    const ZERO: bool = false;
    
    fn from_bits(mut bits: impl Iterator<Item=bool>) -> bool {
        bits.next().unwrap_or(false)
    }

    fn extract(&self, _: usize) -> bool {
        *self
    }

    fn set(&mut self, _: usize, val: bool) {
        *self = val
    }
}

impl Vector for u64 {
    const WIDTH: usize = 64;
    const ZERO: u64 = 0;

    fn from_bits(bits: impl Iterator<Item=bool>) -> u64 {
        let mut out = 0;
        for (i, bit) in bits.take(64).enumerate() {
            out |= (bit as u64) << i;
        }
        out
    }

    fn extract(&self, idx: usize) -> bool {
        (self & (1u64 << idx as u64)) != 0
    }

    fn set(&mut self, idx: usize, val: bool) {
        *self &= !(1u64 << idx);
        *self |= (val as u64) << idx;
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatF2<V: Vector = u64> {
    pitch: usize,
    mpitch: usize,
    pub rows: usize,
    pub cols: usize,
    words: Vec<V>
}

impl<V: Vector> MatF2<V> {
    pub fn new(rows: usize, cols: usize, mut bits: impl Iterator<Item = bool>) -> MatF2<V> {
        let pitch = (cols + V::WIDTH - 1) / V::WIDTH; 
        let mut words = Vec::new();
        for _ in 0..rows {
            let mut row = (&mut bits).take(cols);
            for _ in 0..pitch {
                words.push(V::from_bits(&mut row));
            }
        }

        MatF2 {
            pitch, mpitch: pitch, rows, cols, words
        }
    }

    pub fn build(rows: usize, cols: usize, f: impl Fn(usize, usize) -> bool) -> MatF2<V> {
        Self::new(rows, cols, (0..rows).map(|row| {
            let f = &f;
            (0..cols).map(move |col| {
                f(row, col)
            })
        }).flatten())
    }

    pub fn row_add(&mut self, from: usize, to: usize) {
        let fidx = from * self.pitch;
        let tidx = to * self.pitch;
        
        for word in 0..self.mpitch {
            self.words[tidx + word] = self.words[tidx + word] ^ self.words[fidx + word];
        }
    }

    pub fn row_swap(&mut self, a: usize, b: usize) {
        let aidx = a * self.pitch;
        let bidx = b * self.pitch;
        
        for word in 0..self.mpitch {
            let temp = self.words[aidx + word];
            self.words[aidx + word] = self.words[bidx + word];
            self.words[bidx + word] = temp;
        }
    }

    pub fn col_swap(&mut self, a: usize, b: usize) {
        let aword = a / V::WIDTH;
        let aidx = a % V::WIDTH;
        let bword = b / V::WIDTH;
        let bidx = b % V::WIDTH;

        for row in 0..self.rows {
            let aval = self.words[row * self.pitch + aword].extract(aidx);
            let bval = self.words[row * self.pitch + bword].extract(bidx);
            self.words[row * self.pitch + aword].set(aidx, bval);
            self.words[row * self.pitch + bword].set(bidx, aval);
        }
    }

    fn find_pivot(&mut self, col: usize, row: usize) -> Option<usize> {
        let word = col / V::WIDTH;
        let bit = col % V::WIDTH;
        
        (row..self.rows)
            .find(|i| self.words[i * self.pitch + word].extract(bit))
    }

    fn eliminate_col(&mut self, col: usize, row: usize, pivot: usize) {
        if row != pivot {
            self.row_swap(row, pivot);
        }

        let word = col / V::WIDTH;
        let bit = col % V::WIDTH;

        for i in (row + 1)..self.rows {
            let val = self.words[i * self.pitch + word].extract(bit);
            if val {
                self.row_add(row, i);
            }
        }

        for i in 0..row {
            let val = self.words[i * self.pitch + word].extract(bit);
            if val {
                self.row_add(row, i);
            }
        }
    }

    pub fn gauss_elimination(&mut self) -> usize {
        let mut col = 0;
        let mut row = 0;

        while col < self.cols && row < self.rows {
            if let Some(pivot) = self.find_pivot(col, row) {
                self.eliminate_col(col, row, pivot);
                row += 1;
            }
            col += 1;
        }

        row
    }

    pub fn get(&self, row: usize, col: usize) -> bool {
        self.words[self.pitch * row + col / V::WIDTH].extract(col % V::WIDTH)
    }

    pub fn transpose(&self) -> Self {
        MatF2::build(self.cols, self.rows, |i, j| self.get(j, i))
    }

    pub fn swap_delete_col(&mut self, idx: usize) {
        self.col_swap(idx, self.cols - 1);

        if (self.cols - 1) % V::WIDTH == 0 {
            self.mpitch -= 1;
        }

        self.cols -= 1;
    }

    pub fn add_col(&mut self, val: bool) {
        self.cols += 1;

        if (self.cols - 1) % V::WIDTH == 0  {
            self.mpitch += 1;
        }

        if self.cols - 1 == self.pitch * V::WIDTH  {
            for i in (0..self.rows).rev() {
                let mut element = V::ZERO;
                element.set(0, val);
                self.words.insert((i + 1) * self.pitch, element);
            }

            self.pitch += 1;
        } else {
            let idx = (self.cols - 1) % V::WIDTH;
            let word = (self.cols - 1) / V::WIDTH;
            for i in 0..self.rows {
                self.words[i * self.pitch + word].set(idx, val);
            }
        }
    }

    pub fn row_add_from(&mut self, row: usize, other: &MatF2<V>) {
        for word in 0..self.mpitch.min(other.mpitch) {
            self.words[row * self.pitch + word] = self.words[row * self.pitch + word] ^ other.words[word];
        }
    }

    pub fn add_from(&mut self, other: &MatF2<V>) {
        for row in 0..self.rows.min(other.rows) {
            for word in 0..self.mpitch.min(other.mpitch) {
                self.words[row * self.pitch + word] = self.words[row * self.pitch + word] ^ other.words[row * other.pitch + word];
            }
        }
    }

    pub fn nullspace(&self, k: usize) -> MatF2<V> {
        let mut b = self.clone();
        b.gauss_elimination();
    
        let mut free = Vec::new();
        let mut pivot = Vec::new();
        let mut row = 0;
        let mut col = 0;
        while row < b.rows && col < b.cols {
            if !b.get(row, col) {
                free.push(col);
                col += 1;
            } else {
                pivot.push(col);
                col += 1;
                row += 1;
            }
        }
    
        if row == b.rows {
            free.extend(col..b.cols);
        }
    
        MatF2::build(k.min(free.len()), b.cols, |j, i| {
            if free.binary_search(&i).is_ok() {
                i == free[j]
            } else {
                let k = pivot.binary_search(&i).unwrap();
                b.get(k, free[j])
            }
        })
    }
    
    pub fn rank_decomposition(&self, quiet: bool) -> MatF2<V> {
        fn initial_factor<V: Vector>(a: &MatF2<V>) -> MatF2<V> {
            let n = a.rows;
        
            let n1 = (0..n)
                .filter(|&k| {
                    let sum = (0..n) 
                        .filter(|&j| j != k)
                        .map(|j| a.get(k, j))
                        .fold(false, |a, b| a ^ b);
                    a.get(k, k) != sum
                })
                .collect::<Vec<_>>();
        
            let n2 = (0..n)
                .map(|i| (0..n).map(move |j| (i, j)))
                .flatten()
                .filter(|&(i, j)| i < j && a.get(i, j))
                .collect::<Vec<_>>();
        
            MatF2::build(n, n1.len() + n2.len(), |row, col| {
                if col < n1.len() {
                    row == n1[col]
                } else {
                    let col = col - n1.len();
                    row == n2[col].0 || row == n2[col].1
                }
            })
        }        

        fn improve_factor<V: Vector>(b: &MatF2<V>) -> Option<MatF2<V>> {
            let mut b = b.clone();
            let mut y = b.nullspace(1);
            if y.rows == 0 {
                return None
            }
        
            let sum = (0..b.cols)
                .map(|i| y.get(0, i))
                .fold(false, |a, b| a ^ b);
        
            if sum {
                b.add_col(false);
                y.add_col(true);
            };
        
            let (ca, cb) = (0..y.cols)
                .map(|b| (0..b).map(move |a| (a, b)))
                .flatten()
                .find(|&(a, b)| y.get(0, a) ^ y.get(0, b))?;
        
            for i in 0..b.rows {
                let val = b.get(i, ca) ^ b.get(i, cb);
                if val {
                    b.row_add_from(i, &y);
                }
            }
        
            b.swap_delete_col(cb);
            b.swap_delete_col(ca);
        
            Some(b)
        }

        let a = self;
        let rank = a.clone().gauss_elimination();
        let delta = (0..a.rows)
            .all(|i| !a.get(i, i)) as usize;
        let target = rank + delta;
    
        let mut b = initial_factor(&a);
    
        if !quiet {
            let pb = indicatif::ProgressBar::new((b.cols - target) as u64);
            pb.set_style(indicatif::ProgressStyle::default_bar()
                .template(&format!("[n = {}, {{elapsed:.bold}} elapsed] {{bar:40.cyan/blue}} {{pos:>5}}/{{len}} [{{eta:>3.bold}} remaining]", a.rows))
                .unwrap()
                .progress_chars("##-"));
            let mut prev = b.cols as u64;
    
            while b.cols > target {
                let delta = prev - b.cols as u64;
                prev -= delta;
                pb.inc(delta);
                if let Some(b1) = improve_factor(&b) {
                    b = b1
                } else {
                    break
                }
            }
    
            pb.finish_and_clear();
        } else {
            while b.cols > target {
                if let Some(b1) = improve_factor(&b) {
                    b = b1
                } else {
                    break
                }
            }
        }
    
        b
    }
}

impl<V: Vector> std::fmt::Display for MatF2<V> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.rows {
            write!(f, "[ ")?;
            for j in 0..self.cols {
                let word = j / V::WIDTH;
                let bit = j % V::WIDTH;
                write!(f, "{} ", self.words[i * self.pitch + word].extract(bit) as u8)?;
            }
            write!(f, "]\n")?;
        }

        Ok(())
    }
}

impl<V: Vector> ops::Mul<MatF2<V>> for MatF2<V> {
    type Output = MatF2<V>;

    fn mul(self, rhs: MatF2<V>) -> Self::Output {
        assert_eq!(self.cols, rhs.rows);
        
        MatF2::build(self.rows, rhs.cols, |i, j| {
            (0..self.cols)
                .map(|k| self.get(i, k) & rhs.get(k, j))
                .fold(false, |a, b| a ^ b)
        })
    }
}
