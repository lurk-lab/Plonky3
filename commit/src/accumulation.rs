use alloc::vec::Vec;
use crate::Mmcs;
use core::fmt::Debug;
use core::marker::PhantomData;
use p3_field::{ExtensionField, Field};
use p3_matrix::dense::RowMajorMatrix;
use serde::de::DeserializeOwned;
use serde::Serialize;
use p3_challenger::{CanObserve, CanSample};
use p3_matrix::Matrix;

struct Committed<Val, InputMmcs: Mmcs<Val>> {
    // Merkle root
    commitment: InputMmcs::Commitment,
    // Merkle tree and leaves
    data: InputMmcs::ProverData<RowMajorMatrix<Val>>,
}

const NUM_QUERIES: usize = 32;


struct Proof<Val, Ext, InputMmcs, AccMmcs>
where
    Val: Field,
    Ext: ExtensionField<Val>,
    InputMmcs: Mmcs<Val>,
    AccMmcs: Mmcs<Ext>
{
    input_openings: [(Vec<Val>, InputMmcs::Proof); NUM_QUERIES],
    acc_openings: [(Vec<Ext>, AccMmcs::Proof); NUM_QUERIES],

    out_commitment: AccMmcs::Commitment,
    out_openings: [(Vec<Ext>, AccMmcs::Proof); NUM_QUERIES],
}


pub trait AccumulationScheme<Challenger>
{
    type Input;
    type InputCommitment;

    type Acc;
    type AccCommitment;


    type Proof: Serialize + DeserializeOwned;

    type Error: Debug;

    fn prove(
        &self,
        acc: Self::Acc,
        input: Self::Input,
        challenger: &mut Challenger,
    ) -> Result<(Self::Acc, Self::Proof), Self::Error>;

    fn verify(
        &self,
        acc: Self::AccCommitment,
        input: Self::InputCommitment,
        proof: &Self::Proof,
        challenger: &mut Challenger,
    ) -> Result<Self::AccCommitment, Self::Error>;
}

struct TestAccumulationScheme<Val, Ext, InputMmcs, AccMmcs, Challenger> {
    _marker: PhantomData<(Val, Ext, InputMmcs, AccMmcs, Challenger)>,
}

impl<Val, Ext, InputMmcs, AccMmcs, Challenger> AccumulationScheme<Challenger> for TestAccumulationScheme<Val, Ext, InputMmcs, AccMmcs, Challenger>
where
    Val: Field,
    Ext: ExtensionField<Val>,
    InputMmcs: Mmcs<Val>,
    AccMmcs: Mmcs<Ext>,
    Challenger: CanObserve<AccMmcs::Commitment> + CanSample<Ext>,
{
    type Input = Committed<Val, InputMmcs>;
    type InputCommitment = InputMmcs::Commitment;
    type Acc = Committed<Ext, AccMmcs>;
    type AccCommitment = AccMmcs::Commitment;
    type Proof = Proof<Val, Ext, InputMmcs, AccMmcs>;
    type Error = ();

    fn prove(&self, acc: Self::Acc, input: Self::Input, challenger: &mut Challenger) -> Result<(Self::Acc, Self::Proof), Self::Error> {
        let acc_matrix: &RowMajorMatrix<Ext> = {
            let matrices = AccMmcs::get_matrices(&acc.data);
            assert_eq!(matrices.len(), 1);
            matrices[0]
        };
        let input_matrix: &RowMajorMatrix<Val> = {
            let matrices = InputMmcs::get_matrices(&input.data);
            assert_eq!(matrices.len(), 1);
            matrices[0]
        };

        assert_eq!(acc_matrix.dimensions(), input_matrix.dimensions());

        let alpha = challenger.sample();

        // Compute the random linear combination of the acc and input matrices

        // Commit to the resulting `out` matrix and observe

        // Sample row indices for spot-checking

        // open acc, input, out at the sampled rows

        // return { out.commitment, [acc.proof], [input.proof], [out.proof].

        todo!()
    }

    fn verify(&self, acc: Self::AccCommitment, input: Self::InputCommitment, proof: &Self::Proof, challenger: &mut Challenger) -> Result<Self::AccCommitment, Self::Error> {
        // sample challenge alpha
        // observe proof.out_commitment
        // sample row indices
        // verify openings of acc, input, out
        // verify RLC of each opening
        // return proof.out_commitment
        todo!()
    }
}

