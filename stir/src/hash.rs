use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::merkle_tree::{Config, IdentityDigestConverter};
use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use ark_crypto_primitives::Error;
use ark_ff::MontConfig;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use blake3::Hasher;
use lazy_static::lazy_static;
use rand::{Rng, RngCore};
use sha3::Digest;
use std::borrow::Borrow;
use std::marker::PhantomData;
use std::sync::atomic::AtomicUsize;

#[derive(Default, Clone, Copy)]
pub struct Blake3Config;

pub fn default_fs_config() -> Blake3Config {
    Blake3Config
}

pub type LeafH<F> = SHA3LeafHash<F>;
pub type CompressH = SHA3TwoToOneCRHScheme;

#[derive(Debug, Default, Clone)]
pub struct MerkleTreeParams<F>(PhantomData<F>);

impl<F: CanonicalSerialize + Send> Config for MerkleTreeParams<F> {
    type Leaf = Vec<F>;
    type LeafDigest = <LeafH<F> as CRHScheme>::Output;
    type LeafInnerDigestConverter = IdentityDigestConverter<SHA3Digest>;
    type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;
    type LeafHash = LeafH<F>;
    type TwoToOneHash = CompressH;
}

#[derive(
    Debug, Default, Clone, Copy, Eq, PartialEq, Hash, CanonicalSerialize, CanonicalDeserialize,
)]
pub struct SHA3Digest(pub [u8; 32]);
pub struct SHA3LeafHash<F>(PhantomData<F>);
pub struct SHA3TwoToOneCRHScheme;

pub fn merkle_tree_default_config<F: CanonicalSerialize + Send>(
    rng: &mut impl RngCore,
    _leaf_arity: usize,
) -> (
    <LeafH<F> as CRHScheme>::Parameters,
    <CompressH as TwoToOneCRHScheme>::Parameters,
) {
    let leaf_hash_params = <LeafH<F> as CRHScheme>::setup(rng).unwrap();
    let two_to_one_params = <CompressH as TwoToOneCRHScheme>::setup(rng).unwrap();

    (leaf_hash_params, two_to_one_params)
}

#[derive(Default, Clone)]
pub struct Sponge(Hasher);

impl CryptographicSponge for Sponge {
    type Config = Blake3Config;
    fn new(_config: &Self::Config) -> Self {
        Self::default()
    }
    fn absorb(&mut self, input: &impl Absorb) {
        self.0.update(&input.to_sponge_bytes_as_vec());
    }
    fn squeeze_bytes(&mut self, num_bytes: usize) -> Vec<u8> {
        let mut xof = self.0.finalize_xof();
        let mut output = vec![0u8; num_bytes];
        xof.fill(&mut output);
        self.0.update(&output);
        output
    }
    fn squeeze_bits(&mut self, num_bits: usize) -> Vec<bool> {
        let mut xof = self.0.finalize_xof();
        let mut output = vec![0u8; (num_bits + 7) / 8];
        xof.fill(&mut output);
        self.0.update(&output);
        output
            .iter()
            .flat_map(|byte| (0..8).rev().map(move |i| (byte >> i) & 1 == 1))
            .take(num_bits)
            .collect()
    }
}

impl TwoToOneCRHScheme for SHA3TwoToOneCRHScheme {
    type Input = SHA3Digest;
    type Output = SHA3Digest;
    type Parameters = ();

    fn setup<R: Rng>(_: &mut R) -> Result<Self::Parameters, Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _: &Self::Parameters,
        left_input: T,
        right_input: T,
    ) -> Result<Self::Output, Error> {
        let mut h = sha3::Sha3_256::new();
        h.update(&left_input.borrow().0);
        h.update(&right_input.borrow().0);
        let mut output = [0; 32];
        output.copy_from_slice(&h.finalize()[..]);
        HashCounter::add();
        Ok(SHA3Digest(output))
    }

    fn compress<T: Borrow<Self::Output>>(
        parameters: &Self::Parameters,
        left_input: T,
        right_input: T,
    ) -> Result<Self::Output, Error> {
        <Self as TwoToOneCRHScheme>::evaluate(parameters, left_input, right_input)
    }
}

#[derive(Debug, Default)]
pub struct HashCounter {
    counter: AtomicUsize,
}

lazy_static! {
    static ref HASH_COUNTER: HashCounter = HashCounter::default();
}

impl HashCounter {
    pub fn add() -> usize {
        HASH_COUNTER
            .counter
            .fetch_add(1, std::sync::atomic::Ordering::SeqCst)
    }
    pub fn reset() {
        HASH_COUNTER
            .counter
            .store(0, std::sync::atomic::Ordering::SeqCst)
    }
    pub fn get() -> usize {
        HASH_COUNTER
            .counter
            .load(std::sync::atomic::Ordering::SeqCst)
    }
}

impl Absorb for SHA3Digest {
    fn to_sponge_bytes(&self, dest: &mut Vec<u8>) {
        dest.extend_from_slice(&self.0);
    }

    fn to_sponge_field_elements<F: ark_ff::PrimeField>(&self, dest: &mut Vec<F>) {
        let mut buf = [0; 32];
        buf.copy_from_slice(&self.0);
        dest.push(F::from_be_bytes_mod_order(&buf));
    }
}

impl<F: CanonicalSerialize + Send> CRHScheme for SHA3LeafHash<F> {
    type Input = Vec<F>;
    type Output = SHA3Digest;
    type Parameters = ();

    fn setup<R: Rng>(_: &mut R) -> Result<Self::Parameters, Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _parameters: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, Error> {
        let mut buf = vec![];
        CanonicalSerialize::serialize_compressed(input.borrow(), &mut buf)?;
        let h = sha3::Sha3_256::new();
        let mut output = [0; 32];
        output.copy_from_slice(&h.finalize()[..]);
        HashCounter::add();
        Ok(SHA3Digest(output))
    }
}

#[derive(MontConfig)]
#[modulus = "4787605948707450321761805915146316350821882368518086721537"]
#[generator = "3"]
pub struct FrConfig192;
