use crate::stir_abstractions::{Domain, OracleType, VerificationState};
use crate::utils::quotient_with_hint;
use ark_ff::FftField;
use ark_poly::{
    EvaluationDomain, GeneralEvaluationDomain, MixedRadixEvaluationDomain, Radix2EvaluationDomain,
};

impl<F: FftField> VerificationState<F> {
    // Now, I need to query f_i at a given point.
    // This induces some query to the previous oracle, whose answer I get
    pub fn query(
        &self,
        evaluation_point: F,
        value_of_prev_oracle: F,
        common_factors_inverse: F,
        denom_hint: F,
        ans_eval: F,
    ) -> F {
        match &self.oracle {
            OracleType::Initial => value_of_prev_oracle, // In case this is the initial function, we just return the value of the previous oracle
            OracleType::Virtual(virtual_function) => {
                let num_terms = virtual_function.quotient_set.len();
                let quotient_evaluation = quotient_with_hint(
                    value_of_prev_oracle,
                    evaluation_point,
                    &virtual_function.quotient_set,
                    denom_hint,
                    ans_eval,
                );

                let common_factor = evaluation_point * virtual_function.comb_randomness;

                let scale_factor = if common_factor != F::ONE {
                    (F::ONE - common_factor.pow([(num_terms + 1) as u64])) * common_factors_inverse
                } else {
                    F::from((num_terms + 1) as u64)
                };

                quotient_evaluation * scale_factor
            }
        }
    }
}

impl<F: FftField> Domain<F> {
    pub fn new(degree: usize, log_rho_inv: usize) -> Option<Self> {
        let size = degree * (1 << log_rho_inv);
        let backing_domain = GeneralEvaluationDomain::new(size)?;
        let root_of_unity: F = match backing_domain {
            GeneralEvaluationDomain::Radix2(r2) => r2.group_gen,
            GeneralEvaluationDomain::MixedRadix(mr) => mr.group_gen,
        };
        let root_of_unity_inv = match backing_domain {
            GeneralEvaluationDomain::Radix2(r2) => r2.group_gen_inv,
            GeneralEvaluationDomain::MixedRadix(mr) => mr.group_gen_inv,
        };
        Some(Self {
            backing_domain,
            root_of_unity,
            root_of_unity_inv,
        })
    }

    pub fn size(&self) -> usize {
        self.backing_domain.size()
    }

    // Takes the underlying backing_domain = <w>, and computes the new domain
    // <w^power> (note this will have size |L| / power)
    // NOTE: This should not be mixed with scale_offset
    fn scale_generator_by(&self, power: usize) -> GeneralEvaluationDomain<F> {
        let starting_size = self.size();
        assert_eq!(starting_size % power, 0);
        let new_size = starting_size / power;
        let log_size_of_group = new_size.trailing_zeros();
        let size_as_field_element = F::from(new_size as u64);

        match self.backing_domain {
            GeneralEvaluationDomain::Radix2(r2) => {
                let group_gen = r2.group_gen.pow([power as u64]);
                let group_gen_inv = group_gen.inverse().unwrap();

                let offset = r2.offset.pow([power as u64]);
                let offset_inv = r2.offset_inv.pow([power as u64]);
                let offset_pow_size = offset.pow([new_size as u64]);

                GeneralEvaluationDomain::Radix2(Radix2EvaluationDomain {
                    size: new_size as u64,
                    log_size_of_group,
                    size_as_field_element,
                    size_inv: size_as_field_element.inverse().unwrap(),
                    group_gen,
                    group_gen_inv,
                    offset,
                    offset_inv,
                    offset_pow_size,
                })
            }
            GeneralEvaluationDomain::MixedRadix(mr) => {
                let group_gen = mr.group_gen.pow([power as u64]);
                let group_gen_inv = mr.group_gen_inv.pow([power as u64]);

                let offset = mr.offset.pow([power as u64]);
                let offset_inv = mr.offset_inv.pow([power as u64]);
                let offset_pow_size = offset.pow([new_size as u64]);

                GeneralEvaluationDomain::MixedRadix(MixedRadixEvaluationDomain {
                    size: new_size as u64,
                    log_size_of_group,
                    size_as_field_element,
                    size_inv: size_as_field_element.inverse().unwrap(),
                    group_gen,
                    group_gen_inv,
                    offset,
                    offset_inv,
                    offset_pow_size,
                })
            }
        }
    }

    // Take a domain L_0 = o * <w> and compute a new domain L_1 = w * o^power * <w^power>.
    // Note that L_0^k \cap L_1 = \emptyset for k > power.
    fn scale_with_offset(&self, power: usize) -> GeneralEvaluationDomain<F> {
        let starting_size = self.size();
        assert_eq!(starting_size % power, 0);
        let new_size = starting_size / power;
        let log_size_of_group = new_size.trailing_zeros();
        let size_as_field_element = F::from(new_size as u64);
        match self.backing_domain {
            GeneralEvaluationDomain::Radix2(r2) => {
                let group_gen = r2.group_gen.pow([power as u64]);
                let group_gen_inv = r2.group_gen_inv.pow([power as u64]);

                let offset = r2.offset.pow([power as u64]) * self.root_of_unity;
                let offset_inv = r2.offset_inv.pow([power as u64]) * self.root_of_unity_inv;

                GeneralEvaluationDomain::Radix2(Radix2EvaluationDomain {
                    size: new_size as u64,
                    log_size_of_group,
                    size_as_field_element,
                    size_inv: size_as_field_element.inverse().unwrap(),
                    group_gen,
                    group_gen_inv,
                    offset,
                    offset_inv,
                    offset_pow_size: offset.pow([new_size as u64]),
                })
            }
            GeneralEvaluationDomain::MixedRadix(mr) => {
                let group_gen = mr.group_gen.pow([power as u64]);
                let group_gen_inv = mr.group_gen_inv.pow([power as u64]);

                let offset = mr.offset.pow([power as u64]) * self.root_of_unity;
                let offset_inv = mr.offset_inv.pow([power as u64]) * self.root_of_unity_inv;

                GeneralEvaluationDomain::MixedRadix(MixedRadixEvaluationDomain {
                    size: new_size as u64,
                    log_size_of_group,
                    size_as_field_element,
                    size_inv: size_as_field_element.inverse().unwrap(),
                    group_gen,
                    group_gen_inv,
                    offset,
                    offset_inv,
                    offset_pow_size: offset.pow([new_size as u64]),
                })
            }
        }
    }

    pub fn scale(&self, power: usize) -> Self {
        Self {
            backing_domain: self.scale_generator_by(power),
            ..*self
        }
    }

    pub fn scale_offset(&self, power: usize) -> Self {
        Self {
            backing_domain: self.scale_with_offset(power),
            ..*self
        }
    }
}
