use alloc::vec;
use alloc::vec::Vec;
use core::fmt::Debug;
use core::iter;
use core::marker::PhantomData;

use itertools::{Itertools, zip_eq};

use p3_challenger::{CanObserve, CanSample, CanSampleBits};
use p3_commit::Mmcs;

use p3_field::{ExtensionField, Field};
use p3_matrix::{Dimensions, Matrix};
use p3_matrix::dense::RowMajorMatrix;
use p3_maybe_rayon::prelude::*;
use p3_util::log2_strict_usize;

struct Committed<Val, InputMmcs: Mmcs<Val>>
where Val: Send + Sync {
    // Merkle root
    commitment: InputMmcs::Commitment,
    // Merkle tree and leaves
    data: InputMmcs::ProverData<RowMajorMatrix<Val>>,
}


struct Proof<Val, Ext, InputMmcs, AccMmcs>
where
    Val: Field,
    Ext: ExtensionField<Val>,
    InputMmcs: Mmcs<Val>,
    AccMmcs: Mmcs<Ext>
{
    acc_commit: AccMmcs::Commitment,
    openings: Vec<(Vec<Vec<Val>>, InputMmcs::Proof, AccMmcs::Proof)>,
}

#[derive(Debug)]
enum Error<Val, Ext, InputMmcs, AccMmcs>
where
    Val: Field,
    Ext: ExtensionField<Val>,
    InputMmcs: Mmcs<Val>,
    AccMmcs: Mmcs<Ext> {
    SizeError(),
    InputError(InputMmcs::Error),
    AccError(AccMmcs::Error),
}


pub trait AccumulationScheme<Challenger>
{
    type Input;
    type InputCommitment;

    type Acc;
    type AccCommitment;


    type Proof;
    // type Proof: Serialize + DeserializeOwned;

    type Error;
    // type Error: Debug;

    fn prove(
        &self,
        input: Self::Input,
        challenger: &mut Challenger,
    ) -> Result<(Self::Acc, Self::Proof), Self::Error>;

    fn verify(
        &self,
        input: Self::InputCommitment,
        dimensions: &[Dimensions],
        proof: &Self::Proof,
        challenger: &mut Challenger,
    ) -> Result<Self::AccCommitment, Self::Error>;
}

struct TestAccumulationScheme<Val, Ext, InputMmcs, AccMmcs, Challenger> {
    num_queries: usize,
    input_mmcs: InputMmcs,
    acc_mmcs: AccMmcs,
    _marker: PhantomData<(Val, Ext, Challenger)>,
}

impl<Val, Ext, InputMmcs, AccMmcs, Challenger> AccumulationScheme<Challenger> for TestAccumulationScheme<Val, Ext, InputMmcs, AccMmcs, Challenger>
where
    Val: Field,
    Ext: ExtensionField<Val>,
    InputMmcs: Mmcs<Val>,
    AccMmcs: Mmcs<Ext>,
    Challenger: CanObserve<AccMmcs::Commitment> + CanSample<Ext> + CanSampleBits<usize>,
{
    type Input = Committed<Val, InputMmcs>;
    type InputCommitment = InputMmcs::Commitment;
    type Acc = Committed<Ext, AccMmcs>;
    type AccCommitment = AccMmcs::Commitment;
    type Proof = Proof<Val, Ext, InputMmcs, AccMmcs>;
    type Error = Error<Val, Ext, InputMmcs, AccMmcs>;

    /// Given a commitment to a list of matrices of the same height,
    fn prove(&self, input: Self::Input, challenger: &mut Challenger) -> Result<(Self::Acc, Self::Proof), Self::Error> {
        let matrices = self.input_mmcs.get_matrices(&input.data);

        let alpha: Ext = challenger.sample();

        let height = matrices.iter().map(|mat| mat.height())
            .all_equal_value().map_err(|_| Error::SizeError().into())?;
        let log_height = log2_strict_usize(height);

        // We compute the RLC of the columns of each matrix using powers of alpha.
        // We then combine them into a single column by multiplying by a shifted power of alpha
        // which is given by the alpha raised to the power equal to the first column index
        let alpha_offsets: Vec<Ext> = matrices.iter().scan(0, |running_width, mat| {
            let current_width = *running_width;
            *running_width += mat.width();
            let alpha_offset = alpha.exp_u64(current_width as u64);
            Some(alpha_offset)
        }).collect();

        let mut acc_col = vec![Ext::zero(); height];

        // acc_col = ∑i rlc_col[i] * alpha_offset[i]
        // where rlc_col[i] = ∑_j matrix[i].column(j) * alpha^j
        for (matrix, alpha_offset) in zip_eq(matrices, alpha_offsets) {
            let rlc_col = matrix.dot_ext_powers(alpha);
            acc_col.par_iter_mut().zip_eq(rlc_col).for_each(|(acc, rlc)| {
                *acc += rlc * alpha_offset
            })
        }

        // Commit to acc, the RLC of all columns
        let (acc_commit, acc_data) = self.acc_mmcs.commit(vec![RowMajorMatrix::new_col(acc_col)]);

        // Observe commitment to acc
        challenger.observe(acc_commit.clone());

        // Sample spot-checking indices
        let query_indices = iter::repeat_with(|| challenger.sample_bits(log_height)).take(self.num_queries);

        // Open both input and acc at each query index
        let openings: Vec<_> = query_indices.map(|index| {
            let (input_values, input_proof) = self.input_mmcs.open_batch(index, &input.data);
            // We ignore the opening of the acc column, since this will be derived by the verifier
            let (_, acc_proof) = self.acc_mmcs.open_batch(index, &acc_data);
            (input_values, input_proof, acc_proof)
        }).collect();

        let acc = Committed { commitment: acc_commit.clone(), data: acc_data };
        let proof = Proof {
            acc_commit,
            openings,
        };
        Ok((acc, proof))
    }

    fn verify(
        &self,
        input: Self::InputCommitment,
        dimensions: &[Dimensions],
        proof: &Self::Proof,
        challenger: &mut Challenger,
    ) -> Result<Self::AccCommitment, Self::Error> {
        let Proof { acc_commit, openings } = proof;

        let height = dimensions.iter().map(|dim| dim.height)
            .all_equal_value()
            .map_err(|_| Error::SizeError())?;
        let log_height = log2_strict_usize(height);

        let alpha: Ext = challenger.sample();

        // Observe commitment to acc
        challenger.observe(acc_commit.clone());

        // Sample spot-checking indices
        let query_indices = iter::repeat_with(|| challenger.sample_bits(log_height)).take(self.num_queries);

        let acc_dimensions = [Dimensions { width: 1, height }];

        for (query_index, opening) in query_indices.zip_eq(openings) {
            let (input_values, input_proof, acc_proof) = opening;

            // Compute the expected the value of the opening of acc at query index
            let acc_value = input_values
                .iter()
                .flatten()
                .rev()
                .fold(Ext::zero(), |rlc, input_value| {
                    (rlc * alpha) + Ext::from_base(*input_value)
                });

            // Verify the opened input values
            self.input_mmcs.verify_batch(&input, dimensions, query_index, input_values, input_proof).map_err(Error::InputError)?;

            // Verify the correct computation of the value of acc
            self.acc_mmcs.verify_batch(acc_commit, &acc_dimensions, query_index, &vec![vec![acc_value]], acc_proof).map_err(Error::AccError)?;
        }
        Ok(acc_commit.clone())
    }
}

#[cfg(test)]
mod tests {
    use p3_baby_bear::{BabyBear, DiffusionMatrixBabyBear};
    use p3_challenger::DuplexChallenger;
    use p3_commit::Mmcs;
    use p3_commit::ExtensionMmcs;
    use p3_field::extension::BinomialExtensionField;
    use p3_poseidon2::{Poseidon2, Poseidon2ExternalMatrixGeneral};
    use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
    use p3_merkle_tree::FieldMerkleTreeMmcs;

    use rand::{SeedableRng};
    use rand_chacha::ChaCha20Rng;

    use super::*;

    type Val = BabyBear;
    type Ext = BinomialExtensionField<Val, 4>;

    type Perm = Poseidon2<Val, Poseidon2ExternalMatrixGeneral, DiffusionMatrixBabyBear, 16, 7>;
    type MyHash = PaddingFreeSponge<Perm, 16, 8, 8>;
    type MyCompress = TruncatedPermutation<Perm, 2, 8, 16>;

    type ValMmcs = FieldMerkleTreeMmcs<
        <Val as Field>::Packing,
        <Val as Field>::Packing,
        MyHash,
        MyCompress,
        8,
    >;
    type ExtMmcs = ExtensionMmcs<Val, Ext, ValMmcs>;
    type Challenger = DuplexChallenger<Val, Perm, 16>;

    type MyScheme = TestAccumulationScheme<Val, Ext, ValMmcs, ExtMmcs, Challenger>;

    #[test]
    fn simple() {
        let mut rng = ChaCha20Rng::seed_from_u64(0);

        let perm = Perm::new_from_rng_128(
            Poseidon2ExternalMatrixGeneral,
            DiffusionMatrixBabyBear,
            &mut rng,
        );
        let hash = MyHash::new(perm.clone());
        let compress = MyCompress::new(perm.clone());

        let input_mmcs = ValMmcs::new(hash, compress);
        let acc_mmcs = ExtMmcs::new(input_mmcs.clone());

        let mut challenger_p = Challenger::new(perm.clone());
        let mut challenger_v = Challenger::new(perm.clone());

        let config = MyScheme {
            num_queries: 256,
            input_mmcs,
            acc_mmcs,
            _marker: Default::default(),
        };


        let height: usize = 1 << 16;

        let widths = [1, 4, 9];

        let matrices = widths.into_iter().map(|width| RowMajorMatrix::<Val>::rand(&mut rng, height, width)).collect_vec();
        let dimensions = matrices.iter().map(|mat| mat.dimensions()).collect_vec();

        let (input_commit, input_data) = config.input_mmcs.commit(matrices);
        let input = Committed { commitment: input_commit.clone(), data: input_data };

        let (acc, proof) = config.prove(input, &mut challenger_p).unwrap();

        let acc_commit = config.verify(input_commit, &dimensions, &proof, &mut challenger_v).unwrap();

        assert_eq!(acc.commitment, acc_commit);
    }
}