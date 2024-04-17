use crate::Mmcs;
use core::fmt::Debug;
use p3_field::{ExtensionField, Field};
use p3_matrix::dense::RowMajorMatrix;
use serde::de::DeserializeOwned;
use serde::Serialize;

struct Committed<Val, InputMmcs: Mmcs<Val>> {
    commitment: InputMmcs::Commitment,
    data: InputMmcs::ProverData<RowMajorMatrix<Val>>,
}

pub trait AccumulationScheme<Val, Ext, InputMmcs, AccMmcs, Challenger>
where
    Val: Field,
    Ext: ExtensionField<Val>,
    InputMmcs: Mmcs<Val>,
    AccMmcs: Mmcs<Ext>,
{
    type Proof: Serialize + DeserializeOwned;

    type Error: Debug;
    fn prove(
        &self,
        acc: Committed<Ext, AccMmcs>,
        witness: Committed<Val, InputMmcs>,
        challenger: &mut Challenger,
    ) -> Result<(Committed<Ext, AccMmcs>, Self::Proof), Self::Error>;

    fn verify(
        &self,
        acc: AccMmcs::Commitment,
        witness: InputMmcs::Commitment,
        proof: &Self::Proof,
        challenger: &mut Challenger,
    ) -> Result<AccMmcs::Commitment, Self::Error>;
}
