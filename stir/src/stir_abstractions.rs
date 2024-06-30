use crate::utils::is_power_of_two;
use ark_crypto_primitives::merkle_tree::Config;
use ark_crypto_primitives::merkle_tree::{LeafParam, MerkleTree, MultiPath, TwoToOneParam};
use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use ark_ff::{FftField, PrimeField};
use ark_poly::univariate::DensePolynomial;
use ark_poly::GeneralEvaluationDomain;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use derivative::Derivative;
use std::fmt::Display;
use std::marker::PhantomData;
use std::ops::Deref;

#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub enum SoundnessType {
    Provable,
    Conjecture,
}

#[derive(Derivative)]
#[derivative(Debug, Clone(bound = ""))]
pub struct Parameters<F, MerkleConfig, FSConfig>
where
    F: FftField,
    MerkleConfig: Config,
    FSConfig: CryptographicSponge,
    FSConfig::Config: Clone,
{
    pub security_level: usize,
    pub protocol_security_level: usize,
    pub starting_degree: usize,
    pub stopping_degree: usize,
    pub folding_factor: usize,
    pub starting_rate: usize,
    pub soundness_type: SoundnessType,
    #[derivative(Debug = "ignore")]
    pub leaf_hash_params: LeafParam<MerkleConfig>,
    #[derivative(Debug = "ignore")]
    pub two_to_one_params: TwoToOneParam<MerkleConfig>,
    #[derivative(Debug = "ignore")]
    pub fiat_shamir_config: FSConfig::Config,
    #[derivative(Debug = "ignore")]
    pub _field: PhantomData<F>,
}

impl<F, MerkleConfig, FSConfig> Display for Parameters<F, MerkleConfig, FSConfig>
where
    F: FftField,
    MerkleConfig: Config,
    FSConfig: CryptographicSponge,
    FSConfig::Config: Clone,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let starting_degree_log = (self.starting_degree as f64).log2() as usize;
        let stopping_degree_log = (self.stopping_degree as f64).log2() as usize;

        writeln!(
            f,
            "Targeting {}-bits of security - protocol running at {}-bits - soundness: {:?}",
            self.security_level, self.protocol_security_level, self.soundness_type
        )?;
        writeln!(
            f,
            "Starting degree: 2^{}, stopping_degree: 2^{}",
            starting_degree_log, stopping_degree_log
        )?;
        writeln!(
            f,
            "Starting rate: 2^-{}, folding_factor: {}",
            self.starting_rate, self.folding_factor
        )
    }
}

impl<F, MerkleConfig, FSConfig> Parameters<F, MerkleConfig, FSConfig>
where
    F: FftField,
    MerkleConfig: Config,
    FSConfig: CryptographicSponge,
    FSConfig::Config: Clone,
{
    pub(crate) fn repetitions(&self, log_inv_rate: usize) -> usize {
        let constant = match self.soundness_type {
            SoundnessType::Provable => 2,
            SoundnessType::Conjecture => 1,
        };
        ((constant * self.protocol_security_level) as f64 / log_inv_rate as f64).ceil() as usize
    }

    pub(crate) fn pow_bits(&self, log_inv_rate: usize) -> usize {
        let repetitions = self.repetitions(log_inv_rate);
        // TODO: This will change with eta
        let scaling_factor = match self.soundness_type {
            SoundnessType::Provable => 2.,
            SoundnessType::Conjecture => 1.,
        };
        let achieved_security_bits = (log_inv_rate as f64 / scaling_factor) * repetitions as f64;
        let remaining_security_bits = self.security_level as f64 - achieved_security_bits;

        if remaining_security_bits <= 0. {
            0
        } else {
            remaining_security_bits.ceil() as usize
        }
    }
}

#[derive(Default)]
pub struct Stir<F, MerkleConfig, FSConfig> {
    _field: PhantomData<F>,
    _merkle_config: PhantomData<MerkleConfig>,
    _fs_config: PhantomData<FSConfig>,
}

pub trait LowDegreeTest<F, MerkleConfig, FSConfig>
where
    F: FftField,
    MerkleConfig: Config,
    FSConfig: CryptographicSponge,
    FSConfig::Config: Clone,
{
    type Prover: Prover<
        F,
        MerkleConfig,
        FSConfig,
        Commitment = <Self::Verifier as Verifier<F, MerkleConfig, FSConfig>>::Commitment,
        Proof = <Self::Verifier as Verifier<F, MerkleConfig, FSConfig>>::Proof,
    >;
    type Verifier: Verifier<F, MerkleConfig, FSConfig>;

    fn instantiate(
        parameters: Parameters<F, MerkleConfig, FSConfig>,
    ) -> (Self::Prover, Self::Verifier) {
        let prover = Self::Prover::new(parameters.clone());
        let verifier = Self::Verifier::new(parameters);

        (prover, verifier)
    }

    fn display(parameters: Parameters<F, MerkleConfig, FSConfig>);
}

impl<F, MerkleConfig, FSConfig> LowDegreeTest<F, MerkleConfig, FSConfig>
    for Stir<F, MerkleConfig, FSConfig>
where
    F: FftField + PrimeField + Absorb,
    MerkleConfig: Config<Leaf = Vec<F>>,
    MerkleConfig::InnerDigest: Absorb,
    FSConfig: CryptographicSponge,
    FSConfig::Config: Clone,
{
    type Prover = StirProver<F, MerkleConfig, FSConfig>;
    type Verifier = StirVerifier<F, MerkleConfig, FSConfig>;

    fn display(parameters: Parameters<F, MerkleConfig, FSConfig>) {
        println!("{:?}", FullParameters::from(parameters));
    }
}

pub struct StirProver<F, MerkleConfig, FSConfig>
where
    F: FftField,
    MerkleConfig: Config,
    FSConfig: CryptographicSponge,
    FSConfig::Config: Clone,
{
    pub parameters: FullParameters<F, MerkleConfig, FSConfig>,
}

#[derive(Debug)]
pub struct VerificationState<F: FftField> {
    pub oracle: OracleType<F>,
    pub domain_gen: F,
    pub domain_size: usize,
    pub domain_offset: F,
    pub root_of_unity: F,
    pub folding_randomness: F,
    pub num_round: usize,
}

#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Commitment<MerkleConfig: Config> {
    pub root: MerkleConfig::InnerDigest,
}

#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: FftField, MerkleConfig: Config> {
    pub round_proofs: Vec<RoundProof<F, MerkleConfig>>,
    pub final_polynomial: DensePolynomial<F>,
    pub queries_to_final: (Vec<Vec<F>>, MultiPath<MerkleConfig>),
    pub pow_nonce: Option<usize>,
}

#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct RoundProof<F: FftField, MerkleConfig: Config> {
    pub g_root: MerkleConfig::InnerDigest,
    pub betas: Vec<F>,
    pub ans_polynomial: DensePolynomial<F>,
    pub queries_to_prev: (Vec<Vec<F>>, MultiPath<MerkleConfig>),
    pub shake_polynomial: DensePolynomial<F>,
    pub pow_nonce: Option<usize>,
}

#[derive(Derivative)]
#[derivative(Debug, Clone(bound = ""))]
pub struct FullParameters<F, MerkleConfig, FSConfig>
where
    F: FftField,
    MerkleConfig: Config,
    FSConfig: CryptographicSponge,
    FSConfig::Config: Clone,
{
    #[derivative(Debug(bound = "F: std::fmt::Debug"))]
    pub(crate) parameters: Parameters<F, MerkleConfig, FSConfig>,
    pub(crate) num_rounds: usize,
    pub(crate) rates: Vec<usize>,
    pub(crate) repetitions: Vec<usize>,
    pub(crate) pow_bits: Vec<usize>,
    pub(crate) ood_samples: usize,
    pub(crate) degrees: Vec<usize>,
}

impl<F, MerkleConfig, FSConfig> From<Parameters<F, MerkleConfig, FSConfig>>
    for FullParameters<F, MerkleConfig, FSConfig>
where
    F: FftField,
    MerkleConfig: Config,
    FSConfig: CryptographicSponge,
    FSConfig::Config: Clone,
{
    fn from(parameters: Parameters<F, MerkleConfig, FSConfig>) -> Self {
        assert!(is_power_of_two(parameters.folding_factor));
        assert!(is_power_of_two(parameters.starting_degree));
        assert!(is_power_of_two(parameters.stopping_degree));

        // TODO: I don't even need to iterate but I don't pay for these cycles
        let mut d = parameters.starting_degree;
        let mut degrees = vec![d];
        let mut num_rounds = 0;
        while d > parameters.stopping_degree {
            assert!(d % parameters.folding_factor == 0);
            d /= parameters.folding_factor;
            degrees.push(d);
            num_rounds += 1;
        }

        num_rounds -= 1;
        degrees.pop();

        let mut rates = vec![parameters.starting_rate];
        let log_folding = parameters.folding_factor.ilog2() as usize;
        rates.extend((1..num_rounds + 1).map(|i| parameters.starting_rate + i * (log_folding - 1)));
        let pow_bits: Vec<_> = rates
            .iter()
            .map(|&log_inv_rate| parameters.pow_bits(log_inv_rate))
            .collect();
        let mut repetitions: Vec<_> = rates
            .iter()
            .map(|&log_inv_rate| parameters.repetitions(log_inv_rate))
            .collect();

        // Note, this skips the last repetition
        for i in 0..num_rounds {
            repetitions[i] = repetitions[i].min(degrees[i] / parameters.folding_factor);
        }

        assert_eq!(num_rounds + 1, rates.len());
        assert_eq!(num_rounds + 1, repetitions.len());

        Self {
            parameters,
            num_rounds,
            degrees,
            rates,
            pow_bits,
            ood_samples: 2,
            repetitions,
        }
    }
}

// This is necessary to be able to access fields from both FullParameters and Parameters
impl<F, MerkleConfig, FSConfig> Deref for FullParameters<F, MerkleConfig, FSConfig>
where
    F: FftField,
    MerkleConfig: Config,
    FSConfig: CryptographicSponge,
    FSConfig::Config: Clone,
{
    type Target = Parameters<F, MerkleConfig, FSConfig>;

    fn deref(&self) -> &Parameters<F, MerkleConfig, FSConfig> {
        &self.parameters
    }
}

pub trait Verifier<F, MerkleConfig, FSConfig>
where
    F: FftField,
    MerkleConfig: Config,
    FSConfig: CryptographicSponge,
    FSConfig::Config: Clone,
{
    type FullParameter;

    type Commitment;

    type Proof: CanonicalSerialize;

    fn new(parameters: Parameters<F, MerkleConfig, FSConfig>) -> Self;
    //fn new_full(full_parameters: Self::FullParameter) -> Self;

    fn verify(&self, commitment: &Self::Commitment, proof: &Self::Proof) -> bool;
}

pub trait Prover<F, MerkleConfig, FSConfig>
where
    F: FftField,
    MerkleConfig: Config,
    FSConfig: CryptographicSponge,
    FSConfig::Config: Clone,
{
    type FullParameter;

    type Commitment;

    type Witness: Clone;
    type Proof;

    fn new(parameters: Parameters<F, MerkleConfig, FSConfig>) -> Self;
    fn new_full(full_parameters: Self::FullParameter) -> Self;

    fn commit(&self, polynomial: DensePolynomial<F>) -> (Self::Commitment, Self::Witness);

    fn prove(&self, witness: Self::Witness) -> Self::Proof;
}

impl<F: FftField> Deref for Domain<F> {
    type Target = GeneralEvaluationDomain<F>;

    fn deref(&self) -> &Self::Target {
        &self.backing_domain
    }
}

#[derive(Debug)]
pub struct VirtualFunction<F: FftField> {
    pub comb_randomness: F,
    pub interpolating_polynomial: DensePolynomial<F>,
    pub quotient_set: Vec<F>,
}

#[derive(Debug)]
pub enum OracleType<F: FftField> {
    Initial,
    Virtual(VirtualFunction<F>),
}

#[derive(Derivative)]
#[derivative(Clone(bound = "F: Clone"))]
pub struct Witness<F: FftField, MerkleConfig: Config> {
    pub(crate) domain: Domain<F>,
    pub(crate) polynomial: DensePolynomial<F>,
    pub(crate) merkle_tree: MerkleTree<MerkleConfig>,
    pub(crate) folded_evals: Vec<Vec<F>>,
}

#[derive(Derivative)]
#[derivative(Debug)]
pub struct WitnessExtended<F: FftField, MerkleConfig: Config> {
    #[derivative(Debug = "ignore")]
    pub(crate) domain: Domain<F>,
    pub(crate) polynomial: DensePolynomial<F>,

    #[derivative(Debug = "ignore")]
    pub(crate) merkle_tree: MerkleTree<MerkleConfig>,
    #[derivative(Debug = "ignore")]
    pub(crate) folded_evals: Vec<Vec<F>>,
    pub(crate) num_round: usize,
    pub(crate) folding_randomness: F,
}

pub struct StirVerifier<F, MerkleConfig, FSConfig>
where
    F: FftField,
    MerkleConfig: Config,
    FSConfig: CryptographicSponge,
    FSConfig::Config: Clone,
{
    pub parameters: FullParameters<F, MerkleConfig, FSConfig>,
}

#[derive(Debug, Clone)]
pub struct Domain<F: FftField> {
    pub root_of_unity: F,
    pub root_of_unity_inv: F,
    pub backing_domain: GeneralEvaluationDomain<F>,
}
