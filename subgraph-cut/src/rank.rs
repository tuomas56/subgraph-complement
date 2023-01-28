use quizx::linalg::Mat2;

pub trait RankDecomposition {
    fn rank_decomposition(&self) -> (Mat2, Mat2);
}

impl RankDecomposition for Mat2 {
    fn rank_decomposition(&self) -> (Mat2, Mat2) {
        let mut r = self.clone();
        r.gauss_x(true, self.num_cols(), &mut ());

        let mut pivots = Vec::new();
        let mut row = 0;
        let mut col = 0;
        while row < r.num_rows() && col < r.num_cols() {
            if r[(row, col)] != 0 {
                pivots.push(col);
                row += 1;
            }
            col += 1;
        }

        let mut c = Mat2::new(vec![vec![0; pivots.len()]; self.num_rows()]);
        for (i, p) in pivots.into_iter().enumerate() {
            for j in 0..self.num_rows() {
                c[(j, i)] = self[(j, p)];
            }
        }

        let nonzero = (0..r.num_rows())
            .filter(|&i| (0..r.num_cols()).any(|j| r[(i, j)] != 0))
            .collect::<Vec<_>>();

        let mut f = Mat2::new(vec![vec![0; self.num_cols()]; nonzero.len()]);
        for (i, p) in nonzero.into_iter().enumerate() {
            for j in 0..self.num_cols() {
                f[(i, j)] = r[(p, j)];
            }
        }

        (c, f)
    }
}

#[test]
fn rank_decomposition_test() {
    for _ in 0..1000 {
        let mut m = Mat2::new(vec![vec![0; 10]; 10]);
        for i in 0..10 {
            for j in 0..10 {
                m[(i, j)] = (rand::random::<f32>() < 0.2) as u8;
            }
        }

        let (a, b) = m.rank_decomposition();
        let k = m.rank();
        assert_eq!(a.num_cols(), k);
        assert_eq!(b.num_rows(), k);
        assert_eq!(a * b, m);
    }
}
