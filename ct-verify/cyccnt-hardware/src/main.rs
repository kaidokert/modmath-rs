#![no_main]
#![no_std]

use core::hint::black_box;
use cortex_m_rt::entry;
use krabi_caliper::cortex_m::DwtMeasurementPlatform;
use krabi_caliper::report::Field;
use krabi_caliper::suite::{PairedSuite, PairedSuiteConfig, PairedSuiteFields};

const TRIALS: usize = 4;
const BATCHES: usize = 1;
const MAX_POSITIVE_SPREAD: u32 = 32;

const M64_A: u64 = 0xffff_ffff_ffff_ff43;
const NP64_A: u64 = 0xa53f_a94f_ea53_fa95;
const M64_B: u64 = 0xffff_ffff_ffff_ffc5;
const NP64_B: u64 = 0xcbee_a4e1_a08a_d8f3;

const _: () = assert!(M64_A.wrapping_mul(NP64_A) == u64::MAX);
const _: () = assert!(M64_B.wrapping_mul(NP64_B) == u64::MAX);

const M256_A: [u32; 8] = [
    0xffff_ff43,
    0xffff_ffff,
    0xffff_ffff,
    0xffff_ffff,
    0xffff_ffff,
    0xffff_ffff,
    0xffff_ffff,
    0xffff_ffff,
];
const NP256_A: [u32; 8] = [
    0xea53_fa95,
    0xa53f_a94f,
    0x53fa_94fe,
    0x3fa9_4fea,
    0xfa94_fea5,
    0xa94f_ea53,
    0x94fe_a53f,
    0x4fea_53fa,
];
const M256_B: [u32; 8] = [
    0xffff_ffc5,
    0xffff_ffff,
    0xffff_ffff,
    0xffff_ffff,
    0xffff_ffff,
    0xffff_ffff,
    0xffff_ffff,
    0xffff_ffff,
];
const NP256_B: [u32; 8] = [
    0xa08a_d8f3,
    0xcbee_a4e1,
    0x8682_2b63,
    0x8f2f_ba93,
    0x4e1a_08ad,
    0xb63c_beea,
    0xa938_6822,
    0x8ad8_f2fb,
];

// CIOS consumes n_prime's low word; both pairs satisfy m * n_prime = -1
// modulo the 32-bit word radix. The full arrays above are the corresponding
// -m^-1 values modulo 2^256.
const _: () = assert!(M256_A[0].wrapping_mul(NP256_A[0]) == u32::MAX);
const _: () = assert!(M256_B[0].wrapping_mul(NP256_B[0]) == u32::MAX);

unsafe extern "C" {
    fn ct_fix__wide_mont_mul__u64__N1(
        a: *const u64,
        b: *const u64,
        m: *const u64,
        np: *const u64,
        out: *mut u64,
    );
    fn ct_fix__wide_redc__u64__N1(
        lo: *const u64,
        hi: *const u64,
        m: *const u64,
        np: *const u64,
        out: *mut u64,
    );
    fn ct_fix__cios_mont_mul__u64__N1(
        a: *const u64,
        b: *const u64,
        m: *const u64,
        np: *const u64,
        out: *mut u64,
    );
    fn ct_fix__cios_mont_mul__fb32__N8(
        a: *const [u32; 8],
        b: *const [u32; 8],
        m: *const [u32; 8],
        np: *const [u32; 8],
        out: *mut [u32; 8],
    );
    fn ct_fix__field_exp__u64__N1(base: *const u64, e: *const u64, m: *const u64, out: *mut u64);
    fn ct_fix__field_exp__fb32__N8(
        base: *const [u32; 8],
        e: *const [u32; 8],
        m: *const [u32; 8],
        out: *mut [u32; 8],
    );
    fn ct_fix__field_inv_safegcd__u64__N1(m: *const u64, v: *const u64) -> u8;
    fn ct_fix__field_inv_safegcd__fb32__N8(m: *const [u32; 8], v: *const [u32; 8]) -> u8;
    fn ct_fix__field_blind_path__fb32__N8(
        m: *const [u32; 8],
        x: *const [u32; 8],
        e: *const [u32; 8],
        out: *mut [u32; 8],
    ) -> u8;
    fn ct_fix__field_cswap_eq__fb32__N8(
        m: *const [u32; 8],
        a: *const [u32; 8],
        b: *const [u32; 8],
        choice: *const u8,
        out: *mut [u32; 8],
    ) -> u8;
    fn ct_fix__ASYM__field_inv_safegcd_public_m__u64__N1(v: *const u64) -> u8;
    fn ct_fix__ASYM__field_exp_secret_e__fb32__N8(e: *const [u32; 8], out: *mut [u32; 8]);
    fn ct_fix__ASYM__field_cswap_choice__fb32__N8(choice: *const u8, out: *mut [u32; 8]);
    fn nct_fix__neg__eea_inv__u64__N1(a: *const u64, m: *const u64, out: *mut u64);
    fn nct_fix__neg__schoolbook_exp__u64__N1(
        base: *const u64,
        e: *const u64,
        m: *const u64,
        out: *mut u64,
    );
    fn nct_fix__neg__table_lookup__u64__N1(a: *const u64, out: *mut u64);
}

#[derive(Clone, Copy)]
struct Mont64 {
    a: u64,
    b: u64,
    m: u64,
    np: u64,
}

#[derive(Clone, Copy)]
struct Mont256 {
    a: [u32; 8],
    b: [u32; 8],
    m: [u32; 8],
    np: [u32; 8],
}

#[derive(Clone, Copy)]
struct Exp64 {
    base: u64,
    exp: u64,
    modulus: u64,
}

#[derive(Clone, Copy)]
struct Exp256 {
    base: [u32; 8],
    exp: [u32; 8],
    modulus: [u32; 8],
}

#[derive(Clone, Copy)]
struct Field64 {
    modulus: u64,
    value: u64,
}

#[derive(Clone, Copy)]
struct Field256 {
    modulus: [u32; 8],
    value: [u32; 8],
}

#[derive(Clone, Copy)]
struct Blind256 {
    modulus: [u32; 8],
    value: [u32; 8],
    exponent: [u32; 8],
}

#[derive(Clone, Copy)]
struct Swap256 {
    modulus: [u32; 8],
    a: [u32; 8],
    b: [u32; 8],
    choice: u8,
}

#[inline(never)]
fn fixture_wide_mul(v: &Mont64) -> bool {
    let mut out = 0;
    unsafe { ct_fix__wide_mont_mul__u64__N1(&v.a, &v.b, &v.m, &v.np, &mut out) };
    let _ = black_box(out);
    true
}

#[inline(never)]
fn fixture_wide_redc(v: &Mont64) -> bool {
    let mut out = 0;
    unsafe { ct_fix__wide_redc__u64__N1(&v.a, &v.b, &v.m, &v.np, &mut out) };
    let _ = black_box(out);
    true
}

#[inline(never)]
fn fixture_cios64(v: &Mont64) -> bool {
    let mut out = 0;
    unsafe { ct_fix__cios_mont_mul__u64__N1(&v.a, &v.b, &v.m, &v.np, &mut out) };
    let _ = black_box(out);
    true
}

#[inline(never)]
fn fixture_cios256(v: &Mont256) -> bool {
    let mut out = [0; 8];
    unsafe { ct_fix__cios_mont_mul__fb32__N8(&v.a, &v.b, &v.m, &v.np, &mut out) };
    let _ = black_box(out);
    true
}

#[inline(never)]
fn fixture_modexp64(v: &Exp64) -> bool {
    let mut out = 0;
    unsafe { ct_fix__field_exp__u64__N1(&v.base, &v.exp, &v.modulus, &mut out) };
    let _ = black_box(out);
    true
}

#[inline(never)]
fn fixture_modexp256(v: &Exp256) -> bool {
    let mut out = [0; 8];
    unsafe { ct_fix__field_exp__fb32__N8(&v.base, &v.exp, &v.modulus, &mut out) };
    let _ = black_box(out);
    true
}

#[inline(never)]
fn fixture_field_inv64(v: &Field64) -> bool {
    let result = unsafe { ct_fix__field_inv_safegcd__u64__N1(&v.modulus, &v.value) };
    let _ = black_box(result);
    true
}

#[inline(never)]
fn fixture_field_inv256(v: &Field256) -> bool {
    let result = unsafe { ct_fix__field_inv_safegcd__fb32__N8(&v.modulus, &v.value) };
    let _ = black_box(result);
    true
}

#[inline(never)]
fn fixture_blind256(v: &Blind256) -> bool {
    let mut out = [0; 8];
    let result =
        unsafe { ct_fix__field_blind_path__fb32__N8(&v.modulus, &v.value, &v.exponent, &mut out) };
    let _ = black_box((out, result));
    true
}

#[inline(never)]
fn fixture_cswap_eq256(v: &Swap256) -> bool {
    let mut out = [0; 8];
    let result =
        unsafe { ct_fix__field_cswap_eq__fb32__N8(&v.modulus, &v.a, &v.b, &v.choice, &mut out) };
    let _ = black_box((out, result));
    true
}

#[inline(never)]
fn fixture_asym_inv64(v: &u64) -> bool {
    let result = unsafe { ct_fix__ASYM__field_inv_safegcd_public_m__u64__N1(v) };
    let _ = black_box(result);
    true
}

#[inline(never)]
fn fixture_asym_exp256(v: &[u32; 8]) -> bool {
    let mut out = [0; 8];
    unsafe { ct_fix__ASYM__field_exp_secret_e__fb32__N8(v, &mut out) };
    let _ = black_box(out);
    true
}

#[inline(never)]
fn fixture_asym_cswap(v: &u8) -> bool {
    let mut out = [0; 8];
    unsafe { ct_fix__ASYM__field_cswap_choice__fb32__N8(v, &mut out) };
    let _ = black_box(out);
    true
}

#[inline(never)]
fn fixture_neg_eea(v: &Field64) -> bool {
    let mut out = 0;
    unsafe { nct_fix__neg__eea_inv__u64__N1(&v.value, &v.modulus, &mut out) };
    let _ = black_box(out);
    true
}

#[inline(never)]
fn fixture_neg_exp(v: &Exp64) -> bool {
    let mut out = 0;
    unsafe { nct_fix__neg__schoolbook_exp__u64__N1(&v.base, &v.exp, &v.modulus, &mut out) };
    let _ = black_box(out);
    true
}

#[inline(never)]
fn fixture_table(v: &u64) -> bool {
    let mut out = 0;
    unsafe { nct_fix__neg__table_lookup__u64__N1(v, &mut out) };
    let _ = black_box(out);
    true
}

#[entry]
fn main() -> ! {
    let mut reporter = krabi_caliper::protocol::rtt::init_ct_compatible();
    ct_fixtures::link_anchor();
    let mut peripherals = cortex_m::Peripherals::take().unwrap();
    let mut platform = DwtMeasurementPlatform::enable(
        &mut peripherals.DCB,
        &mut peripherals.DWT,
        Some(16_000_000),
    )
    .unwrap();

    let mont64_a = Mont64 {
        a: 5,
        b: 7,
        m: M64_A,
        np: NP64_A,
    };
    let mont64_b = Mont64 {
        a: 0x1234_5678_9abc_def0,
        b: 0xfedc_ba98_7654_3210,
        m: M64_B,
        np: NP64_B,
    };
    let mont256_a = Mont256 {
        a: [5, 0, 0, 0, 0, 0, 0, 0],
        b: [7, 0, 0, 0, 0, 0, 0, 0],
        m: M256_A,
        np: NP256_A,
    };
    let mont256_b = Mont256 {
        a: [0x7654_3210, 0xfedc_ba98, 3, 4, 5, 6, 7, 8],
        b: [0x89ab_cdef, 0x0123_4567, 8, 7, 6, 5, 4, 3],
        m: M256_B,
        np: NP256_B,
    };
    let exp64_a = Exp64 {
        base: 3,
        exp: 3,
        modulus: M64_B,
    };
    let exp64_b = Exp64 {
        base: 0x1234_5678_9abc_def0,
        exp: 0xa55a_3cc3_9669_f00f,
        modulus: M64_B,
    };
    let exp256_a = Exp256 {
        base: [3, 0, 0, 0, 0, 0, 0, 0],
        exp: [3, 0, 0, 0, 0, 0, 0, 0],
        modulus: M256_B,
    };
    let exp256_b = Exp256 {
        base: [0x7654_3210, 0xfedc_ba98, 3, 4, 5, 6, 7, 8],
        exp: [0xa55a_3cc3, 0x9669_f00f, 2, 4, 8, 16, 32, 64],
        modulus: M256_B,
    };
    let field64_a = Field64 {
        modulus: M64_A,
        value: 5,
    };
    let field64_b = Field64 {
        modulus: M64_B,
        value: 0x1234_5678_9abc_def1,
    };
    let field256_a = Field256 {
        modulus: M256_A,
        value: mont256_a.a,
    };
    let field256_b = Field256 {
        modulus: M256_B,
        value: mont256_b.a,
    };
    let blind_a = Blind256 {
        modulus: M256_A,
        value: mont256_a.a,
        exponent: exp256_a.exp,
    };
    let blind_b = Blind256 {
        modulus: M256_B,
        value: mont256_b.a,
        exponent: exp256_b.exp,
    };
    let swap_a = Swap256 {
        modulus: M256_A,
        a: mont256_a.a,
        b: mont256_a.b,
        choice: 0,
    };
    let swap_b = Swap256 {
        modulus: M256_B,
        a: mont256_b.a,
        b: mont256_b.b,
        choice: 1,
    };
    let asym_e_a = exp256_a.exp;
    let asym_e_b = exp256_b.exp;

    let run_fields = [
        Field::u64("trials", TRIALS as u64),
        Field::u64("max_positive_spread", MAX_POSITIVE_SPREAD as u64),
    ];
    let summary_fields = [Field::u64("diagnostics", 1)];
    let mut suite = PairedSuite::<_, _, TRIALS>::start(
        &mut platform,
        &mut reporter,
        PairedSuiteConfig {
            suite: "modmath-cyccnt",
            target: "thumbv7em-none-eabihf",
            board: Some("stm32f407vg"),
            unit: krabi_caliper::Unit::CoreCycles,
            frequency_hz: Some(16_000_000),
            warmup_blocks: 1,
            batches: BATCHES,
            positive_max_spread: MAX_POSITIVE_SPREAD as u64,
            positive_require_overlap: false,
            fields: PairedSuiteFields {
                run: &run_fields,
                fixture: &[],
                summary: &summary_fields,
            },
        },
    )
    .unwrap();
    macro_rules! positive {
        ($name:literal, $a:expr, $b:expr, $fixture:expr) => {
            suite.positive($name, $a, $b, $fixture).unwrap()
        };
    }
    positive!("wide_mont_mul_u64", &mont64_a, &mont64_b, fixture_wide_mul);
    positive!("wide_redc_u64", &mont64_a, &mont64_b, fixture_wide_redc);
    positive!("cios_mont_mul_u64", &mont64_a, &mont64_b, fixture_cios64);
    positive!(
        "cios_mont_mul_fb32",
        &mont256_a,
        &mont256_b,
        fixture_cios256
    );
    positive!("mod_exp_public_m_u64", &exp64_a, &exp64_b, fixture_modexp64);
    positive!(
        "mod_exp_public_m_fb32",
        &exp256_a,
        &exp256_b,
        fixture_modexp256
    );
    positive!(
        "field_inv_secret_m_u64",
        &field64_a,
        &field64_b,
        fixture_field_inv64
    );
    positive!(
        "field_inv_secret_m_fb32",
        &field256_a,
        &field256_b,
        fixture_field_inv256
    );
    positive!(
        "field_blind_path_fb32",
        &blind_a,
        &blind_b,
        fixture_blind256
    );
    positive!("field_cswap_eq_fb32", &swap_a, &swap_b, fixture_cswap_eq256);
    positive!(
        "asym_inv_public_m_u64",
        &5,
        &0x1234_5678_9abc_def1,
        fixture_asym_inv64
    );
    positive!(
        "asym_exp_secret_e_fb32",
        &asym_e_a,
        &asym_e_b,
        fixture_asym_exp256
    );
    positive!("asym_cswap_choice_fb32", &0, &1, fixture_asym_cswap);
    suite
        .negative(
            "negative_eea_inverse",
            &Field64 {
                modulus: M64_B,
                value: 1,
            },
            &field64_b,
            fixture_neg_eea,
        )
        .unwrap();
    suite
        .negative(
            "negative_schoolbook_exp",
            &Exp64 {
                base: 3,
                exp: 1,
                modulus: M64_B,
            },
            &exp64_b,
            fixture_neg_exp,
        )
        .unwrap();
    suite
        .diagnostic(
            "negative_table_lookup",
            "address-only",
            &0,
            &15,
            fixture_table,
        )
        .unwrap();
    suite.finish().unwrap();
    loop {
        cortex_m::asm::nop();
    }
}
