use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ff::{FftField, Field};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, Evaluations, Polynomial, Radix2EvaluationDomain};
use std::collections::BTreeSet;

pub struct BivariatePolynomial<F: Field>(pub Vec<Vec<F>>);

impl<F> BivariatePolynomial<F>
where
    F: Field,
{
    pub fn fold_by_col(&self, alpha: F) -> DensePolynomial<F> {
        let transposed = transpose(self.0.clone());

        let mut res = DensePolynomial::from_coefficients_vec(vec![]);

        let mut pow = F::ONE;
        for c in transposed {
            res += &DensePolynomial::from_coefficients_vec(c.iter().map(|f| pow * f).collect());
            pow *= alpha;
        }

        res
    }
}

pub fn transpose<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>> {
    assert!(!v.is_empty());
    let len = v[0].len();
    let mut iters: Vec<_> = v.into_iter().map(|n| n.into_iter()).collect();
    (0..len)
        .map(|_| {
            iters
                .iter_mut()
                .map(|n| n.next().unwrap())
                .collect::<Vec<T>>()
        })
        .collect()
}

// Compute the quotient
pub fn poly_quotient<F: FftField>(poly: &DensePolynomial<F>, points: &[F]) -> DensePolynomial<F> {
    let evaluations: Vec<_> = points.iter().map(|x| (*x, poly.evaluate(x))).collect();
    let ans_polynomial = naive_interpolation(evaluations.iter());
    let vanishing_poly = vanishing_poly(points);
    let numerator = poly + &ans_polynomial;

    // TODO: Is this efficient or should FFT?
    &numerator / &vanishing_poly
}

// Computes a polynomial that vanishes on points
pub fn vanishing_poly<'a, F: Field>(points: impl IntoIterator<Item = &'a F>) -> DensePolynomial<F> {
    // Compute the denominator (which is \prod_a(x - a))
    let mut vanishing_poly: DensePolynomial<_> =
        DensePolynomial::from_coefficients_slice(&[F::ONE]);
    for a in points {
        vanishing_poly =
            vanishing_poly.naive_mul(&DensePolynomial::from_coefficients_slice(&[-*a, F::ONE]));
    }
    vanishing_poly
}

// Computes a polynomial that interpolates the given points with the given answers
pub fn naive_interpolation<'a, F: Field>(
    points: impl IntoIterator<Item = &'a (F, F)>,
) -> DensePolynomial<F> {
    let points: Vec<_> = points.into_iter().collect();
    let vanishing_poly = vanishing_poly(points.iter().map(|(a, _)| a));

    // Compute the ans polynomial (this is just a naive interpolation)
    let mut ans_polynomial = DensePolynomial::from_coefficients_slice(&[]);
    for (a, eval) in points.iter() {
        // Computes the vanishing (apart from x - a)
        let vanishing_adjusted =
            &vanishing_poly / &DensePolynomial::from_coefficients_slice(&[-*a, F::ONE]);

        // Now, we can scale to get the right weigh
        let scale_factor = *eval / vanishing_adjusted.evaluate(a);
        ans_polynomial = ans_polynomial
            + DensePolynomial::from_coefficients_vec(
                vanishing_adjusted
                    .iter()
                    .map(|x| *x * scale_factor)
                    .collect(),
            );
    }
    ans_polynomial
}

// Allows to amortize the evaluation of the quotient polynomial
pub fn quotient_with_hint<'a, F: Field>(
    claimed_eval: F,
    evaluation_point: F,
    quotient_set: impl IntoIterator<Item = &'a F>,
    //ans_polynomial: &DensePolynomial<F>,
    denom_hint: F,
    ans_eval: F,
) -> F {
    let quotient_set: Vec<_> = quotient_set.into_iter().copied().collect();

    // Check if the evaluation point is in the domain
    for dom in quotient_set.iter() {
        if evaluation_point == *dom {
            panic!("Evaluation point is in the domain");
        }
    }

    let num = claimed_eval - ans_eval;

    num * denom_hint
}

pub fn poly_fold<F: Field>(
    f: &DensePolynomial<F>,
    folding_factor: usize,
    folding_randomness: F,
) -> DensePolynomial<F> {
    let degree = f.degree() + 1;
    let q_poly = to_coefficient_matrix(f, degree.div_ceil(folding_factor), folding_factor);
    q_poly.fold_by_col(folding_randomness)
}

// Given a generator and a coset offset, computes the interpolating offset
// Requires to be given the inversion of the generator and coset offset (and thus can be more
// efficient)
pub fn fft_interpolate<'a, F: FftField>(
    generator: F,
    coset_offset: F,
    generator_inv: F,
    coset_offset_inv: F,
    size_inv: F,
    points: impl IntoIterator<Item = &'a F>,
) -> DensePolynomial<F> {
    let points: Vec<_> = points.into_iter().cloned().collect();
    let folding_factor = points.len();
    assert!(is_power_of_two(folding_factor));

    let size_as_field_element = F::from(folding_factor as u64);

    let domain = Radix2EvaluationDomain {
        size: folding_factor as u64,
        log_size_of_group: folding_factor.ilog2(),
        size_as_field_element,
        size_inv,
        group_gen: generator,
        group_gen_inv: generator_inv,
        offset: coset_offset,
        offset_inv: coset_offset_inv,
        offset_pow_size: coset_offset.pow([folding_factor as u64]),
    };

    let evaluations = Evaluations::from_vec_and_domain(points, domain);

    evaluations.interpolate()
}

pub fn is_power_of_two(n: usize) -> bool {
    n & (n - 1) == 0
}

// Takes the vector of evaluations (assume that evals[i] = f(omega^i))
// and folds them into a vector of such that folded_evals[i] = [f(omega^(i + k * j)) for j in 0..folding_factor]
pub fn stack_evaluations<F: Copy>(evals: Vec<F>, folding_factor: usize) -> Vec<Vec<F>> {
    assert!(evals.len() % folding_factor == 0);
    let size_of_new_domain = evals.len() / folding_factor;

    let mut stacked_evaluations = vec![];
    for i in 0..size_of_new_domain {
        let mut new_evals = vec![];
        for j in 0..folding_factor {
            new_evals.push(evals[i + j * size_of_new_domain]);
        }
        stacked_evaluations.push(new_evals);
    }

    stacked_evaluations
}

// Deduplicates AND orders a vector
pub fn dedup<T: Ord>(v: impl IntoIterator<Item = T>) -> Vec<T> {
    Vec::from_iter(BTreeSet::from_iter(v))
}

pub fn squeeze_integer(sponge: &mut impl CryptographicSponge, range: usize) -> usize {
    assert!(is_power_of_two(range));
    let mut bytes_array = [0; 8];
    let bytes = sponge.squeeze_bytes(8);
    bytes_array.copy_from_slice(&bytes);
    let candidate = usize::from_le_bytes(bytes_array);
    // This is uniform as long as the range is a power of two
    candidate % range
}

pub fn proof_of_work_verify(
    sponge: &mut impl CryptographicSponge,
    proof_of_work_bits: usize,
    pow_nonce: Option<usize>,
) -> bool {
    assert!(proof_of_work_bits <= 32);
    if proof_of_work_bits == 0 {
        return true;
    }

    if pow_nonce.is_none() {
        return false;
    }
    let nonce = pow_nonce.unwrap();
    sponge.absorb(&nonce.to_le_bytes().as_slice());
    let pow_bytes = sponge.squeeze_bytes(4);
    let mut buf = [0; 4];
    buf.copy_from_slice(&pow_bytes[..]);
    let pow = u32::from_le_bytes(buf);
    pow.trailing_zeros() as usize >= proof_of_work_bits
}

pub fn to_coefficient_matrix<F: Field>(
    f: &DensePolynomial<F>,
    rows: usize,
    cols: usize,
) -> BivariatePolynomial<F> {
    if f.degree() + 1 > rows * cols {
        panic!("Degree of polynomial is too large for matrix");
    }

    let mut matrix = vec![vec![F::ZERO; cols]; rows];

    for (i, coeff) in f.coeffs.iter().enumerate() {
        matrix[i / cols][i % cols] = *coeff;
    }

    BivariatePolynomial(matrix)
}

pub fn proof_of_work(
    sponge: &mut impl CryptographicSponge,
    proof_of_work_bits: usize,
) -> Option<usize> {
    assert!(proof_of_work_bits <= 32);
    if proof_of_work_bits == 0 {
        return None;
    }

    let mut buf = [0; 4];
    let mut nonce: usize = 0;
    loop {
        let mut new_sponge = sponge.clone();
        let nonce_bytes = nonce.to_le_bytes();
        new_sponge.absorb(&nonce_bytes.as_slice());
        let pow_bytes = new_sponge.squeeze_bytes(4);
        buf.copy_from_slice(&pow_bytes[..]);
        let pow = u32::from_le_bytes(buf);
        if pow.trailing_zeros() as usize >= proof_of_work_bits {
            sponge.absorb(&nonce_bytes.as_slice());
            sponge.squeeze_bytes(4);
            return Some(nonce);
        }
        nonce += 1;
    }
}
