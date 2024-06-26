use alloc::vec::Vec;

use icicle_babybear::field::ScalarField;
use icicle_core::ntt::{initialize_domain, ntt_inplace, release_domain, NTTConfig, NTTDir};
use icicle_cuda_runtime::device_context::DeviceContext;
use icicle_cuda_runtime::memory::HostSlice;
use p3_dft::TwoAdicSubgroupDft;
use p3_field::{AbstractField, PrimeField32, TwoAdicField};
use p3_matrix::dense::RowMajorMatrix;
use p3_matrix::Matrix;
use p3_util::log2_strict_usize;

use crate::BabyBear;

type F = BabyBear;

/// GPU DFT/NTT for BabyBear with Icicle
#[derive(Debug, Default, Clone)]
pub struct BabyBearIcicleDft;

impl TwoAdicSubgroupDft<F> for BabyBearIcicleDft {
    type Evaluations = RowMajorMatrix<F>;
    fn dft_batch(&self, mat: RowMajorMatrix<F>) -> RowMajorMatrix<F> {
        let ctx = DeviceContext::default();

        let h = mat.height();
        let w = mat.width();
        let log_h = log2_strict_usize(h);
        let root = F::two_adic_generator(log_h);

        // To compute FFTs using icicle, we first need to initialize it using plonky3's "two adic generator"
        initialize_domain(ScalarField::from([root.as_canonical_u32()]), &ctx, false).unwrap();
        let mut scalars: Vec<ScalarField> = mat
            .values
            .iter()
            .map(|x| ScalarField::from([x.as_canonical_u32()]))
            .collect();
        // let scalars_p3: Vec<BabyBear> = scalars
        //     .iter()
        //     .map(|x| BabyBear::from_wrapped_u32(Into::<[u32; 1]>::into(*x)[0]))
        //     .collect();
        // let matrix_p3 = RowMajorMatrix::new(scalars_p3, w);
        // assert_eq!(matrix_p3, mat);

        let mut ntt_cfg: NTTConfig<'_, ScalarField> = NTTConfig::default();
        // Next two lines signalize that we want to compute `nrows` FFTs in column-ordered fashion
        ntt_cfg.batch_size = w as i32;
        ntt_cfg.columns_batch = true;
        ntt_inplace(
            HostSlice::from_mut_slice(&mut scalars[..]),
            NTTDir::kForward,
            &ntt_cfg,
        )
        .unwrap();

        let result_p3: Vec<BabyBear> = scalars
            .iter()
            .map(|x| BabyBear::from_wrapped_u32(Into::<[u32; 1]>::into(*x)[0]))
            .collect();
        let matrix_p3 = RowMajorMatrix::new(result_p3, w);
        // let result_p3_test = Radix2DitParallel.dft_batch(mat);
        // assert_eq!(matrix_p3, result_p3_test.to_row_major_matrix());

        release_domain::<ScalarField>(&ctx).unwrap();
        matrix_p3
    }
}