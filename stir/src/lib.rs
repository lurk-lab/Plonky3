mod common;
mod hash;
mod prover;
mod stir_abstractions;
mod utils;
mod verifier;

#[cfg(test)]
mod test {
    use ark_ff::{Fp192, MontBackend};
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_serialize::CanonicalSerialize;
    use std::time::Instant;

    use crate::hash::{
        default_fs_config, merkle_tree_default_config, Blake3Config, FrConfig192, HashCounter,
        MerkleTreeParams, Sponge,
    };
    use crate::stir_abstractions::{
        LowDegreeTest, Parameters, Prover, SoundnessType, Stir, Verifier,
    };

    pub type Field192 = Fp192<MontBackend<FrConfig192, 3>>;

    #[test]
    fn test_stir_codebase_estimation() {
        let security_level = 128usize;
        let protocol_security_level = 106usize;
        let initial_degree = 20usize;
        let final_degree = 6usize;
        let rate = 2usize;
        let verifier_repetitions = 1000usize;
        let folding_factor = 16usize;
        let soundness_type = SoundnessType::Conjecture;
        let fiat_shamir_config: Blake3Config = default_fs_config();

        let mut rng = ark_std::test_rng();
        let poly = DensePolynomial::<Field192>::rand(initial_degree - 1, &mut rng);
        let (leaf_hash_params, two_to_one_params) =
            merkle_tree_default_config::<Field192>(&mut rng, folding_factor);

        let params: Parameters<Field192, MerkleTreeParams<Field192>, Sponge> = Parameters {
            security_level,
            protocol_security_level,
            starting_degree: 1 << initial_degree,
            stopping_degree: 1 << final_degree,
            folding_factor,
            starting_rate: rate,
            soundness_type,

            leaf_hash_params: leaf_hash_params.clone(),
            two_to_one_params: two_to_one_params.clone(),
            fiat_shamir_config,
            _field: Default::default(),
        };

        Stir::display(params.clone());

        let (prover, verifier) = Stir::instantiate(params);
        let (commitment, witness) = prover.commit(poly.clone());

        let stir_prover_time = Instant::now();
        let proof = prover.prove(witness);
        dbg!(stir_prover_time.elapsed());
        dbg!(proof.serialized_size(ark_serialize::Compress::Yes));
        let prover_hashes = HashCounter::get();
        dbg!(prover_hashes);
        HashCounter::reset();

        let stir_verifier_time = Instant::now();
        for _ in 0..verifier_repetitions {
            let result = verifier.verify(&commitment, &proof);
            assert!(result);
        }
        dbg!(stir_verifier_time.elapsed());
        let verifier_hashes = HashCounter::get() / verifier_repetitions;
        dbg!(verifier_hashes);
        HashCounter::reset();
    }
}
