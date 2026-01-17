
pub const FPSMOD: u32 = 998_244_353;
pub const PRIMITIVE_ROOT: u32 = 3;
const MOD2: u32 = FPSMOD * 2;

const fn mont_r() -> u32 {
    // r ≡ mod^{-1} (mod 2^32)
    let mut ret = FPSMOD;
    let mut i = 0;
    while i < 4 {
        ret = ret.wrapping_mul(2u32.wrapping_sub(FPSMOD.wrapping_mul(ret)));
        i += 1;
    }
    ret
}
const R: u32 = mont_r();
const NEG_R: u32 = 0u32.wrapping_sub(R);
const N2: u32 = ((0u64.wrapping_sub(FPSMOD as u64)) % (FPSMOD as u64)) as u32; // 2^64 mod FPSMOD

#[inline(always)]
const fn reduce_const(b: u64) -> u32 {
    let t = (b as u32).wrapping_mul(NEG_R) as u64;
    let u = (b + t * (FPSMOD as u64)) >> 32;
    u as u32
}
#[inline(always)]
fn reduce(b: u64) -> u32 {
    reduce_const(b)
}

use std::cmp::min;
use std::ops::*;

#[derive(Copy, Clone, Eq, PartialEq, Debug, Default)]
pub struct Mint(pub u32); // Montgomery form, lazy range ~ [0,2*mod)

impl Mint {
    pub const ZERO: Mint = Mint(0);
    pub const ONE: Mint = Mint(reduce_const(N2 as u64)); // new(1)

    #[inline(always)]
    pub fn new<T: Into<i64>>(v: T) -> Self {
        let mut x = v.into() % (FPSMOD as i64);
        if x < 0 {
            x += FPSMOD as i64;
        }
        Mint(reduce((x as u64) * (N2 as u64)))
    }

    #[inline(always)]
    pub fn get(self) -> u32 {
        let mut ret = reduce(self.0 as u64);
        if ret >= FPSMOD {
            ret -= FPSMOD;
        }
        ret
    }

    #[inline(always)]
    pub fn pow(self, mut e: u64) -> Self {
        let mut a = self;
        let mut r = Mint::ONE;
        while e > 0 {
            if e & 1 == 1 {
                r *= a;
            }
            a *= a;
            e >>= 1;
        }
        r
    }

    #[inline(always)]
    pub fn inv(self) -> Self {
        // FPSMOD is prime
        self.pow((FPSMOD as u64) - 2)
    }
}

impl std::fmt::Display for Mint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.get())
    }
}

impl std::ops::Add for Mint {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self::Output {
        let mut s = self.0 as u64 + rhs.0 as u64;
        if s >= MOD2 as u64 {
            s -= MOD2 as u64;
        }
        Mint(s as u32)
    }
}
impl std::ops::Sub for Mint {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self::Output {
        let mut s = self.0 as i64 - rhs.0 as i64;
        if s < 0 {
            s += MOD2 as i64;
        }
        Mint(s as u32)
    }
}
impl std::ops::Mul for Mint {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self::Output {
        Mint(reduce((self.0 as u64) * (rhs.0 as u64)))
    }
}
impl std::ops::Div for Mint {
    type Output = Self;
    #[inline(always)]
    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inv()
    }
}
impl std::ops::AddAssign for Mint {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}
impl std::ops::SubAssign for Mint {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}
impl std::ops::MulAssign for Mint {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}
impl std::ops::DivAssign for Mint {
    #[inline(always)]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}
impl std::ops::Neg for Mint {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        if self.0 == 0 {
            self
        } else {
            Mint(MOD2 - self.0)
        }
    }
}

// =======================
//  NTT (ACL style butterfly)
// =======================

#[inline(always)]
fn is_pow2(n: usize) -> bool {
    n != 0 && (n & (n - 1)) == 0
}

#[derive(Clone)]
struct NttTables {
    rate2: Vec<Mint>,
    irate2: Vec<Mint>,
    rate3: Vec<Mint>,
    irate3: Vec<Mint>,
    imag: Mint,
    iimag: Mint,
    rank2: usize,
}

use std::sync::OnceLock;
fn ntt_tables() -> &'static NttTables {
    static ONCE: OnceLock<NttTables> = OnceLock::new();
    ONCE.get_or_init(|| {
        let rank2 = (FPSMOD - 1).trailing_zeros() as usize; // 23 for 998244353
        let root = Mint::new(PRIMITIVE_ROOT as i64).pow(((FPSMOD - 1) >> rank2) as u64); // primitive 2^rank2-th root
        let iroot = root.inv();

        // imag = root^(2^(rank2-2)) is primitive 4th root
        let imag = root.pow(1u64 << (rank2 - 2));
        let iimag = imag.inv();

        let mut rate2 = vec![Mint::ZERO; rank2 - 1];
        let mut irate2 = vec![Mint::ZERO; rank2 - 1];
        {
            let mut prod = Mint::ONE;
            let mut iprod = Mint::ONE;
            for i in 0..(rank2 - 1) {
                let e = 1u64 << (rank2 - i - 2);
                let w = root.pow(e);
                let iw = iroot.pow(e);
                rate2[i] = w * prod;
                irate2[i] = iw * iprod;
                prod *= iw;
                iprod *= w;
            }
        }

        let mut rate3 = vec![Mint::ZERO; rank2 - 2];
        let mut irate3 = vec![Mint::ZERO; rank2 - 2];
        {
            let mut prod = Mint::ONE;
            let mut iprod = Mint::ONE;
            for i in 0..(rank2 - 2) {
                let e = 1u64 << (rank2 - i - 3);
                let w = root.pow(e);
                let iw = iroot.pow(e);
                rate3[i] = w * prod;
                irate3[i] = iw * iprod;
                prod *= iw;
                iprod *= w;
            }
        }

        NttTables {
            rate2,
            irate2,
            rate3,
            irate3,
            imag,
            iimag,
            rank2,
        }
    })
}

#[inline(always)]
fn butterfly(a: &mut [Mint]) {
    let n = a.len();
    if n <= 1 {
        return;
    }
    debug_assert!(is_pow2(n));
    let h = n.trailing_zeros() as usize;
    let tab = ntt_tables();
    debug_assert!(h <= tab.rank2);

    let mut len = 0usize;
    while len < h {
        if h - len == 1 {
            let p = 1usize << len;
            let mut rot = Mint::ONE;
            for s in 0..(1usize << (h - len - 1)) {
                let offset = s << (len + 1);
                for i in 0..p {
                    let l = a[offset + i];
                    let r = a[offset + i + p] * rot;
                    a[offset + i] = l + r;
                    a[offset + i + p] = l - r;
                }
                let idx = (!s).trailing_zeros() as usize;
                rot *= tab.rate2[idx];
            }
            len += 1;
        } else {
            let p = 1usize << len;
            let mut rot = Mint::ONE;
            for s in 0..(1usize << (h - len - 2)) {
                let rot2 = rot * rot;
                let rot3 = rot2 * rot;
                let offset = s << (len + 2);
                for i in 0..p {
                    let a0 = a[offset + i];
                    let a1 = a[offset + i + p] * rot;
                    let a2 = a[offset + i + 2 * p] * rot2;
                    let a3 = a[offset + i + 3 * p] * rot3;

                    let a1na3 = a1 - a3;
                    let a1pa3 = a1 + a3;
                    let a0pa2 = a0 + a2;
                    let a0ma2 = a0 - a2;

                    a[offset + i] = a0pa2 + a1pa3;
                    a[offset + i + p] = a0pa2 - a1pa3;
                    a[offset + i + 2 * p] = a0ma2 + a1na3 * tab.imag;
                    a[offset + i + 3 * p] = a0ma2 - a1na3 * tab.imag;
                }
                let idx = (!s).trailing_zeros() as usize;
                rot *= tab.rate3[idx];
            }
            len += 2;
        }
    }
}

#[inline(always)]
fn butterfly_inv(a: &mut [Mint]) {
    let n = a.len();
    if n <= 1 {
        return;
    }
    debug_assert!(is_pow2(n));
    let h = n.trailing_zeros() as usize;
    let tab = ntt_tables();
    debug_assert!(h <= tab.rank2);

    let mut len = h;
    while len > 0 {
        if len == 1 {
            let p = 1usize << (h - len);
            let mut irot = Mint::ONE;
            for s in 0..(1usize << (len - 1)) {
                let offset = s << (h - len + 1);
                for i in 0..p {
                    let l = a[offset + i];
                    let r = a[offset + i + p];
                    a[offset + i] = l + r;
                    a[offset + i + p] = (l - r) * irot;
                }
                let idx = (!s).trailing_zeros() as usize;
                irot *= tab.irate2[idx];
            }
            len -= 1;
        } else {
            let p = 1usize << (h - len);
            let mut irot = Mint::ONE;
            for s in 0..(1usize << (len - 2)) {
                let irot2 = irot * irot;
                let irot3 = irot2 * irot;
                let offset = s << (h - len + 2);
                for i in 0..p {
                    let y0 = a[offset + i];
                    let y1 = a[offset + i + p];
                    let y2 = a[offset + i + 2 * p];
                    let y3 = a[offset + i + 3 * p];

                    let t0 = y0 + y1;
                    let t1 = y0 - y1;
                    let t2 = y2 + y3;
                    let t3 = (y2 - y3) * tab.iimag;

                    a[offset + i] = t0 + t2;
                    a[offset + i + p] = (t1 + t3) * irot;
                    a[offset + i + 2 * p] = (t0 - t2) * irot2;
                    a[offset + i + 3 * p] = (t1 - t3) * irot3;
                }
                let idx = (!s).trailing_zeros() as usize;
                irot *= tab.irate3[idx];
            }
            len -= 2;
        }
    }
}

#[inline(always)]
pub fn ntt_inplace(a: &mut [Mint]) {
    butterfly(a);
}

#[inline(always)]
pub fn intt_inplace(a: &mut [Mint]) {
    butterfly_inv(a);
    let inv_n = Mint::new(a.len() as i64).inv();
    for x in a.iter_mut() {
        *x *= inv_n;
    }
}

#[inline]
fn ceil_pow2(mut n: usize) -> usize {
    let mut p = 1usize;
    while p < n {
        p <<= 1;
    }
    p
}

pub fn convolution(a: &[Mint], b: &[Mint]) -> Vec<Mint> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    if min(a.len(), b.len()) <= 40 {
        let mut res = vec![Mint::ZERO; a.len() + b.len() - 1];
        for i in 0..a.len() {
            for j in 0..b.len() {
                res[i + j] += a[i] * b[j];
            }
        }
        return res;
    }
    let need = a.len() + b.len() - 1;
    let n = ceil_pow2(need);
    let mut fa = vec![Mint::ZERO; n];
    let mut fb = vec![Mint::ZERO; n];
    fa[..a.len()].copy_from_slice(a);
    fb[..b.len()].copy_from_slice(b);

    ntt_inplace(&mut fa);
    ntt_inplace(&mut fb);
    for i in 0..n {
        fa[i] *= fb[i];
    }
    butterfly_inv(&mut fa);
    let inv_n = Mint::new(n as i64).inv();
    for i in 0..need {
        fa[i] *= inv_n;
    }
    fa.truncate(need);
    fa
}

// =======================
//  Binomial (fac/finv/inv/C)
// =======================

pub struct Binomial {
    f: Vec<Mint>, // fac
    g: Vec<Mint>, // finv
    h: Vec<Mint>, // inv
}

impl Binomial {
    pub fn new(max_n: usize) -> Self {
        let mut b = Binomial {
            f: vec![Mint::ONE],
            g: vec![Mint::ONE],
            h: vec![Mint::ZERO], // inv[0] unused
        };
        if max_n > 0 {
            b.extend(max_n + 1);
        }
        b
    }

    pub fn extend(&mut self, mut m: usize) {
        m = min(m, FPSMOD as usize);
        let n = self.f.len();
        if n >= m {
            return;
        }
        self.f.resize(m, Mint::ZERO);
        self.g.resize(m, Mint::ZERO);
        self.h.resize(m, Mint::ZERO);

        for i in n..m {
            self.f[i] = self.f[i - 1] * Mint::new(i as i64);
        }
        self.g[m - 1] = self.f[m - 1].inv();
        for i in (n..=(m - 2)).rev() {
            self.g[i] = self.g[i + 1] * Mint::new((i + 1) as i64);
        }
        // inv[i] = fac[i-1]*finv[i]  (for i>=1)
        self.h[1] = Mint::ONE;
        for i in max(2, n)..m {
            self.h[i] = self.f[i - 1] * self.g[i];
        }
    }

    #[inline(always)]
    pub fn fac(&mut self, i: isize) -> Mint {
        if i < 0 {
            return Mint::ZERO;
        }
        let i = i as usize;
        if i >= FPSMOD as usize {
            return Mint::ZERO;
        }
        if i >= self.f.len() {
            self.extend(i + 1);
        }
        self.f[i]
    }

    #[inline(always)]
    pub fn finv(&mut self, i: isize) -> Mint {
        if i < 0 {
            return Mint::ZERO;
        }
        let i = i as usize;
        if i >= FPSMOD as usize {
            return Mint::ZERO;
        }
        if i >= self.g.len() {
            self.extend(i + 1);
        }
        self.g[i]
    }

    #[inline(always)]
    pub fn inv(&mut self, i: isize) -> Mint {
        if i < 0 {
            return -self.inv(-i);
        }
        let i = i as usize;
        if i == 0 {
            return Mint::ZERO;
        }
        if i >= FPSMOD as usize {
            return Mint::ZERO;
        }
        if i >= self.h.len() {
            self.extend(i + 1);
        }
        self.h[i]
    }

    #[inline(always)]
    pub fn c(&mut self, n: isize, r: isize) -> Mint {
        if n < 0 || r < 0 || r > n {
            return Mint::ZERO;
        }
        self.fac(n) * self.finv(n - r) * self.finv(r)
    }
}

// =======================
//  FPS helpers
// =======================

pub type Fps = Vec<Mint>;
pub type Matrix = Vec<Vec<Mint>>;

#[inline(always)]
pub fn shrink(a: &mut Fps) {
    while let Some(&x) = a.last() {
        if x == Mint::ZERO {
            a.pop();
        } else {
            break;
        }
    }
}

#[inline(always)]
pub fn pre(a: &[Mint], n: usize) -> Fps {
    let mut r = a[..min(a.len(), n)].to_vec();
    if r.len() < n {
        r.resize(n, Mint::ZERO);
    }
    r
}

#[inline(always)]
pub fn rev(a: &[Mint]) -> Fps {
    let mut r = a.to_vec();
    r.reverse();
    r
}

#[inline(always)]
pub fn shl(a: &[Mint], k: usize) -> Fps {
    if a.is_empty() {
        return vec![];
    }
    let mut r = vec![Mint::ZERO; k];
    r.extend_from_slice(a);
    r
}

#[inline(always)]
pub fn shr(a: &[Mint], k: usize) -> Fps {
    if a.len() <= k {
        vec![]
    } else {
        a[k..].to_vec()
    }
}

#[inline(always)]
pub fn add_fps(a: &[Mint], b: &[Mint]) -> Fps {
    let mut r = a.to_vec();
    if b.len() > r.len() {
        r.resize(b.len(), Mint::ZERO);
    }
    for i in 0..b.len() {
        r[i] += b[i];
    }
    r
}

#[inline(always)]
pub fn sub_fps(a: &[Mint], b: &[Mint]) -> Fps {
    let mut r = a.to_vec();
    if b.len() > r.len() {
        r.resize(b.len(), Mint::ZERO);
    }
    for i in 0..b.len() {
        r[i] -= b[i];
    }
    r
}

#[inline(always)]
pub fn mul_scalar(a: &[Mint], v: Mint) -> Fps {
    let mut r = a.to_vec();
    for x in r.iter_mut() {
        *x *= v;
    }
    r
}

#[inline(always)]
pub fn diff(a: &[Mint]) -> Fps {
    if a.len() <= 1 {
        return vec![];
    }
    let mut r = vec![Mint::ZERO; a.len() - 1];
    for i in 1..a.len() {
        r[i - 1] = a[i] * Mint::new(i as i64);
    }
    r
}

#[inline(always)]
pub fn integral(a: &[Mint]) -> Fps {
    let mut r = vec![Mint::ZERO; a.len() + 1];
    for i in 0..a.len() {
        r[i + 1] = a[i] / Mint::new((i + 1) as i64);
    }
    r
}

#[inline(always)]
pub fn eval(a: &[Mint], x: Mint) -> Mint {
    let mut r = Mint::ZERO;
    let mut w = Mint::ONE;
    for &c in a.iter() {
        r += w * c;
        w *= x;
    }
    r
}

// series inverse (Newton)
pub fn inv_series(a: &[Mint], deg: usize) -> Fps {
    assert!(deg >= 1);
    assert!(!a.is_empty() && a[0] != Mint::ZERO);

    let mut res = vec![Mint::ZERO; deg];
    res[0] = a[0].inv();
    let mut m = 1usize;
    while m < deg {
        let m2 = min(deg, m * 2);
        let f = pre(a, m2);
        let g = res[..m].to_vec();
        let fg = convolution(&f, &g);
        let mut t = vec![Mint::ZERO; m2];
        t[0] = Mint::new(2);
        for i in 0..min(m2, fg.len()) {
            t[i] -= fg[i];
        }
        let ng = convolution(&g, &t);
        for i in 0..m2 {
            res[i] = ng.get(i).copied().unwrap_or(Mint::ZERO);
        }
        m = m2;
    }
    res
}

pub fn log_series(a: &[Mint], deg: usize) -> Fps {
    assert!(deg >= 1);
    assert!(!a.is_empty() && a[0] == Mint::ONE);
    let da = diff(a);
    let ia = inv_series(a, deg);
    let mut r = convolution(&da, &ia);
    r.truncate(deg.saturating_sub(1));
    let mut r = integral(&r);
    r.truncate(deg);
    r
}

pub fn exp_series(a: &[Mint], deg: usize) -> Fps {
    assert!(deg >= 1);
    assert!(a.is_empty() || a[0] == Mint::ZERO);

    let mut g = vec![Mint::ONE];
    let mut m = 1usize;
    while m < deg {
        let m2 = min(deg, m * 2);
        let lg = log_series(&pre(&g, m2), m2);
        let mut diffv = pre(a, m2);
        for i in 0..m2 {
            diffv[i] -= lg.get(i).copied().unwrap_or(Mint::ZERO);
        }
        diffv[0] += Mint::ONE;
        let ng = convolution(&g, &diffv);
        g = pre(&ng, m2);
        m = m2;
    }
    g.truncate(deg);
    g
}

pub fn pow_series(a: &[Mint], k: i64, deg: usize) -> Fps {
    if deg == 0 {
        return vec![];
    }
    if k == 0 {
        let mut r = vec![Mint::ZERO; deg];
        r[0] = Mint::ONE;
        return r;
    }
    let n = a.len();
    let mut i0 = 0usize;
    while i0 < n && a[i0] == Mint::ZERO {
        i0 += 1;
    }
    if i0 == n {
        return vec![Mint::ZERO; deg];
    }
    let shift = (i0 as i128) * (k as i128);
    if shift >= deg as i128 {
        return vec![Mint::ZERO; deg];
    }

    let a0 = a[i0];
    let inv_a0 = a0.inv();
    let a0k = a0.pow(k.unsigned_abs() as u64);
    let scale = if k > 0 { a0k } else { a0k.inv() };

    let mut base = mul_scalar(&shr(a, i0), inv_a0);
    // base[0] should be 1
    assert!(!base.is_empty() && base[0] == Mint::ONE);
    base.resize(deg, Mint::ZERO);

    let mut lg = log_series(&base, deg);
    let kk = Mint::new(k as i64);
    for x in lg.iter_mut() {
        *x *= kk;
    }
    let mut r = exp_series(&lg, deg);
    r = mul_scalar(&r, scale);
    let sh = shift as usize;
    let mut r = shl(&r, sh);
    r.truncate(deg);
    r
}

// poly division / remainder (fast via reverse+inv; small via long division)
pub fn div_rem(a: &[Mint], b: &[Mint]) -> (Fps, Fps) {
    let mut aa = a.to_vec();
    let mut bb = b.to_vec();
    shrink(&mut aa);
    shrink(&mut bb);
    if bb.is_empty() {
        panic!("division by zero polynomial");
    }
    if aa.len() < bb.len() {
        return (vec![], aa);
    }

    if bb.len() <= 64 {
        let n = aa.len();
        let m = bb.len();
        let inv_lead = bb[m - 1].inv();
        let mut q = vec![Mint::ZERO; n - m + 1];
        for k in (0..=n - m).rev() {
            let coef = aa[k + m - 1] * inv_lead;
            q[k] = coef;
            for j in 0..m {
                aa[k + j] -= coef * bb[j];
            }
        }
        aa.truncate(m - 1);
        shrink(&mut aa);
        shrink(&mut q);
        return (q, aa);
    }

    let n = aa.len();
    let m = bb.len();
    let deg_q = n - m + 1;

    let mut ar = rev(&aa);
    let mut br = rev(&bb);
    ar.truncate(deg_q);
    br.truncate(deg_q);

    let inv_br = inv_series(&br, deg_q);
    let mut qr = convolution(&ar, &inv_br);
    qr.truncate(deg_q);
    let mut q = rev(&qr);
    shrink(&mut q);

    let mut prod = convolution(&bb, &q);
    prod.truncate(n);
    let mut r = a.to_vec();
    if prod.len() > r.len() {
        r.resize(prod.len(), Mint::ZERO);
    }
    for i in 0..prod.len() {
        r[i] -= prod[i];
    }
    r.truncate(m - 1);
    shrink(&mut r);
    (q, r)
}

// =======================
//  multiply2d (same style)
// =======================

fn transpose(mat: &Matrix) -> Matrix {
    let h = mat.len();
    let w = mat[0].len();
    let mut t = vec![vec![Mint::ZERO; h]; w];
    for i in 0..h {
        for j in 0..w {
            t[j][i] = mat[i][j];
        }
    }
    t
}

fn add_shifted_row(dst: &mut Vec<Mint>, src: &[Mint], shift: usize) {
    if dst.len() < src.len() + shift {
        dst.resize(src.len() + shift, Mint::ZERO);
    }
    for j in 0..src.len() {
        dst[j + shift] += src[j];
    }
}

fn multiply2d_naive(a: &Matrix, b: &Matrix) -> Matrix {
    if a.is_empty() || b.is_empty() || a[0].is_empty() || b[0].is_empty() {
        return vec![];
    }
    let ha = a.len();
    let wa = a[0].len();
    let hb = b.len();
    let wb = b[0].len();
    let mut c = vec![vec![Mint::ZERO; wa + wb - 1]; ha + hb - 1];
    for ia in 0..ha {
        for ja in 0..wa {
            for ib in 0..hb {
                for jb in 0..wb {
                    c[ia + ib][ja + jb] += a[ia][ja] * b[ib][jb];
                }
            }
        }
    }
    c
}

fn multiply2d_partially_naive(mut a: Matrix, mut b: Matrix) -> Matrix {
    if a.is_empty() || b.is_empty() || a[0].is_empty() || b[0].is_empty() {
        return vec![];
    }
    let mut ha = a.len();
    let mut wa = a[0].len();
    let mut hb = b.len();
    let mut wb = b[0].len();

    if min(ha, hb) * min(wa, wb) <= 40 {
        return multiply2d_naive(&a, &b);
    }

    let mut w = 1usize;
    while w < wa + wb - 1 {
        w <<= 1;
    }

    // split trick (same idea)
    if w >= 64 && wa + wb - 1 <= w / 2 + 20 {
        if wa <= 20 {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut ha, &mut hb);
            std::mem::swap(&mut wa, &mut wb);
        }
        let d = wa + wb - 1 - w / 2;
        let mut a1 = vec![vec![Mint::ZERO; wa - d]; ha];
        let mut a2 = vec![vec![Mint::ZERO; d]; ha];
        for i in 0..ha {
            a1[i].copy_from_slice(&a[i][..wa - d]);
            a2[i].copy_from_slice(&a[i][wa - d..]);
        }
        let mut c1 = multiply2d_partially_naive(a1, b.clone());
        let c2 = multiply2d_partially_naive(a2, b);
        for i in 0..(ha + hb - 1) {
            add_shifted_row(&mut c1[i], &c2[i], wa - d);
            c1[i].truncate(wa + wb - 1);
        }
        return c1;
    }

    for row in a.iter_mut() {
        row.resize(w, Mint::ZERO);
        ntt_inplace(row);
    }
    for row in b.iter_mut() {
        row.resize(w, Mint::ZERO);
        ntt_inplace(row);
    }

    // for each frequency column, convolve heights
    let mut c_t: Matrix = Vec::with_capacity(w);
    for j in 0..w {
        let mut col_a = vec![Mint::ZERO; ha];
        let mut col_b = vec![Mint::ZERO; hb];
        for i in 0..ha {
            col_a[i] = a[i][j];
        }
        for i in 0..hb {
            col_b[i] = b[i][j];
        }
        c_t.push(convolution(&col_a, &col_b)); // length ha+hb-1
    }

    let mut c = transpose(&c_t); // (ha+hb-1) x w
    for row in c.iter_mut() {
        intt_inplace(row);
        row.truncate(wa + wb - 1);
    }
    c
}

pub fn multiply2d(mut a: Matrix, mut b: Matrix) -> Matrix {
    if a.is_empty() || b.is_empty() || a[0].is_empty() || b[0].is_empty() {
        return vec![];
    }
    let ha = a.len();
    let wa = a[0].len();
    let hb = b.len();
    let wb = b[0].len();

    if min(ha, hb) * min(wa, wb) <= 40 {
        return multiply2d_naive(&a, &b);
    }
    if min(ha, hb) <= 40 {
        return multiply2d_partially_naive(a, b);
    }
    if min(wa, wb) <= 40 {
        let at = transpose(&a);
        let bt = transpose(&b);
        let ct = multiply2d_partially_naive(at, bt);
        return transpose(&ct);
    }

    let mut h = 1usize;
    let mut w = 1usize;
    while h < ha + hb - 1 {
        h <<= 1;
    }
    while w < wa + wb - 1 {
        w <<= 1;
    }

    // split width trick
    if wa + wb - 1 < w / 2 + 20 {
        let d = wa + wb - 1 - w / 2;
        let mut a1 = vec![vec![Mint::ZERO; wa - d]; ha];
        let mut a2 = vec![vec![Mint::ZERO; d]; ha];
        for i in 0..ha {
            a1[i].copy_from_slice(&a[i][..wa - d]);
            a2[i].copy_from_slice(&a[i][wa - d..]);
        }
        let mut c1 = multiply2d(a1, b.clone());
        let c2 = multiply2d(a2, b);
        for i in 0..(ha + hb - 1) {
            add_shifted_row(&mut c1[i], &c2[i], wa - d);
            c1[i].truncate(wa + wb - 1);
        }
        return c1;
    }
    // split height trick
    if ha + hb - 1 < h / 2 + 20 {
        let at = transpose(&a);
        let bt = transpose(&b);
        let ct = multiply2d(at, bt);
        return transpose(&ct);
    }

    // pad to HxW
    a.resize(h, vec![Mint::ZERO; wa]);
    for row in a.iter_mut() {
        row.resize(w, Mint::ZERO);
    }
    b.resize(h, vec![Mint::ZERO; wb]);
    for row in b.iter_mut() {
        row.resize(w, Mint::ZERO);
    }

    // 2D NTT: rows then cols
    for i in 0..h {
        ntt_inplace(&mut a[i]);
        ntt_inplace(&mut b[i]);
    }
    // columns
    let mut col = vec![Mint::ZERO; h];
    for j in 0..w {
        for i in 0..h {
            col[i] = a[i][j];
        }
        ntt_inplace(&mut col);
        for i in 0..h {
            a[i][j] = col[i];
        }

        for i in 0..h {
            col[i] = b[i][j];
        }
        ntt_inplace(&mut col);
        for i in 0..h {
            b[i][j] = col[i];
        }
    }

    // pointwise
    for i in 0..h {
        for j in 0..w {
            a[i][j] *= b[i][j];
        }
    }

    // inverse cols then rows
    for j in 0..w {
        for i in 0..h {
            col[i] = a[i][j];
        }
        intt_inplace(&mut col);
        for i in 0..h {
            a[i][j] = col[i];
        }
    }
    for i in 0..h {
        intt_inplace(&mut a[i]);
    }

    a.truncate(ha + hb - 1);
    for row in a.iter_mut() {
        row.truncate(wa + wb - 1);
    }
    a
}

// =======================
//  middle_product / pow_enumerate / inner / composite / compositional_inverse
// =======================

pub fn middle_product(a: &[Mint], c: &[Mint]) -> Fps {
    let s = a.len();
    let t = c.len();
    assert!(s > 0 && s <= t);
    let bsz = ceil_pow2(t);

    let mut ar = a.to_vec();
    ar.reverse();
    ar.resize(bsz, Mint::ZERO);
    let mut cr = c.to_vec();
    cr.resize(bsz, Mint::ZERO);

    ntt_inplace(&mut ar);
    ntt_inplace(&mut cr);
    for i in 0..bsz {
        ar[i] *= cr[i];
    }
    butterfly_inv(&mut ar);
    let inv_n = Mint::new(bsz as i64).inv();
    for x in ar.iter_mut() {
        *x *= inv_n;
    }
    ar[(s - 1)..t].to_vec()
}

// [x^n] f(x)^i g(x) を i=0..m で列挙 (C++ と同じ)
pub fn pow_enumerate(mut f: Fps, mut g: Fps, mut m: isize) -> Fps {
    let mut n = f.len() - 1;
    let mut k = 1usize;
    g.resize(n + 1, Mint::ZERO);
    if m == -1 {
        m = n as isize;
    }

    // P, Q: (n+1) x (k+1)
    let mut p: Matrix = vec![vec![Mint::ZERO; k + 1]; n + 1];
    let mut q: Matrix = vec![vec![Mint::ZERO; k + 1]; n + 1];
    q[0][0] = Mint::ONE;
    for i in 0..=n {
        p[i][0] = g[i];
        if k >= 1 {
            q[i][1] = -f[i];
        }
    }

    while n > 0 {
        let mut r = q.clone();
        for i in (1..=n).step_by(2) {
            for j in 0..=k {
                r[i][j] = -r[i][j];
            }
        }
        let s = multiply2d(p, r.clone());
        let t = multiply2d(q, r);

        let mut u: Matrix = vec![vec![Mint::ZERO; k * 2 + 1]; n / 2 + 1];
        let mut v: Matrix = vec![vec![Mint::ZERO; k * 2 + 1]; n / 2 + 1];

        for i in 0..=n / 2 {
            u[i] = s[i * 2 + (n & 1)].clone();
            v[i] = t[i * 2].clone();
        }
        p = u;
        q = v;
        n /= 2;
        k *= 2;
    }

    // return (P[0] * Q[0].inv(m+1)).pre(m+1)
    let deg = (m as usize) + 1;
    let invq = inv_series(&q[0], deg);
    let mut res = convolution(&p[0], &invq);
    res.truncate(deg);
    res
}

fn inner(g: &[Mint], q: Matrix, n: usize, k: usize) -> Matrix {
    if n == 0 {
        // h = g * Q[0].inv().rev()
        let invq = inv_series(&q[0], q[0].len());
        let mut invq_rev = invq;
        invq_rev.reverse();
        let h = convolution(g, &invq_rev);

        let mut p = vec![vec![Mint::ZERO; k + 1]; 1];
        // P[0][i] = h[i + Q[0].size() - 1]
        let off = q[0].len() - 1;
        for i in 0..=k {
            p[0][i] = h[i + off];
        }
        return p;
    }

    let mut r = q.clone();
    for i in (1..=n).step_by(2) {
        for j in 0..r[i].len() {
            r[i][j] = -r[i][j];
        }
    }
    let t = multiply2d(q, r.clone());

    let mut v: Matrix = vec![vec![Mint::ZERO; k * 2 + 1]; n / 2 + 1];
    for i in 0..=n / 2 {
        v[i] = t[i * 2].clone();
    }
    let u = inner(g, v, n / 2, k * 2);

    let mut s: Matrix = vec![vec![Mint::ZERO; k * 2 + 1]; n * 2 + 1];
    for i in 0..=n / 2 {
        s[i * 2 + (n & 1)] = u[i].clone();
    }

    // flatten S and R
    let w2 = k * 2 + 1;
    let mut s2 = vec![Mint::ZERO; (n * 2 + 1) * w2];
    for i in 0..=n * 2 {
        for j in 0..=k * 2 {
            s2[i * w2 + j] = s[i][j];
        }
    }

    let mut r2 = vec![Mint::ZERO; n * w2 + (k + 1)];
    for i in 0..=n {
        for j in 0..=k {
            r2[i * w2 + j] = r[i][j];
        }
    }

    let p2 = middle_product(&r2, &s2);
    let mut p: Matrix = vec![vec![Mint::ZERO; k + 1]; n + 1];
    for i in 0..=n {
        for j in 0..=k {
            p[i][j] = p2[i * w2 + j];
        }
    }
    p
}

// g(f(x)) を計算 (C++ の composite と同じ)
pub fn composite(f: &[Mint], g: &[Mint]) -> Fps {
    assert_eq!(f.len(), g.len());
    let n = f.len() - 1;
    let k = 1usize;

    let mut q: Matrix = vec![vec![Mint::ZERO; 2]; n + 1];
    q[0][0] = Mint::ONE;
    for i in 0..=n {
        q[i][1] = -f[i];
    }

    let p = inner(g, q, n, k);
    let mut h = vec![Mint::ZERO; n + 1];
    for i in 0..=n {
        h[i] = p[i][0];
    }
    h.reverse();
    h
}

// 高速 compositional inverse (C++ と同じ式)
pub fn compositional_inverse(f: &[Mint], deg: usize) -> Fps {
    assert!(deg >= 1);
    assert!(f.len() >= 2);
    assert!(f[1] != Mint::ZERO);

    if deg < 2 {
        return pre(&[Mint::ZERO], deg);
    }
    let n = deg - 1;

    let mut h = pow_enumerate(f.to_vec(), vec![Mint::ONE], -1); // size n+1
                                                                // h *= n
    let nn = Mint::new(n as i64);
    for x in h.iter_mut() {
        *x *= nn;
    }
    // for k=1..=n: h[k]/=k
    for k in 1..=n {
        h[k] /= Mint::new(k as i64);
    }
    h.reverse();
    // normalize h[0] to 1
    let inv0 = h[0].inv();
    for x in h.iter_mut() {
        *x *= inv0;
    }
    // g = exp( log(h) * (-n)^{-1} )
    let logh = log_series(&h, n + 1);
    let inv_minus_n = (-Mint::new(n as i64)).inv();
    let mut t = mul_scalar(&logh, inv_minus_n);
    let mut g = exp_series(&t, n + 1);

    // * f[1]^{-1}
    g = mul_scalar(&g, f[1].inv());
    // <<1, pre(deg)
    let mut ans = shl(&g, 1);
    ans.truncate(deg);
    ans
}

// =======================
//  Fast Composition Q(P(x))  (C++ Composition port)
// =======================

pub fn composition(mut p: Fps, mut q: Fps, c: &mut Binomial, deg: Option<usize>) -> Fps {
    // deg default: min(|p|,|q|)
    let mut n = deg.unwrap_or_else(|| min(p.len(), q.len()));
    if n == 0 {
        return vec![];
    }
    shrink(&mut p);
    if p.is_empty() {
        let mut r = vec![Mint::ZERO; n];
        r[0] = if q.is_empty() { Mint::ZERO } else { q[0] };
        return r;
    }
    if n == 1 {
        return vec![eval(&q, p[0])];
    }

    p.resize(n, Mint::ZERO);
    q.resize(n, Mint::ZERO);

    // M = max(2, sqrt(N / log2(N)))
    let nf = n as f64;
    let m = max(2usize, ((nf / nf.log2()).sqrt().floor() as usize));
    let l = (n + m - 1) / m;

    let pm = p[..min(m, n)].to_vec();
    let pr = if n > m { p[m..].to_vec() } else { vec![] };

    // J = ceil(log2 N)
    let j = (usize::BITS as usize) - (n - 1).leading_zeros() as usize;

    let mut pms: Vec<Fps> = vec![vec![]; j];
    pms[0] = pm.clone();
    for i in 1..j {
        let sq = convolution(&pms[i - 1], &pms[i - 1]);
        pms[i] = pre(&sq, n);
    }

    fn comp_rec(pms: &Vec<Fps>, q: &Vec<Mint>, n: usize, left: usize, j: usize) -> Fps {
        if j == 1 {
            if left >= n {
                return vec![];
            }
            let q1 = q[left];
            let q2 = if left + 1 < n {
                q[left + 1]
            } else {
                Mint::ZERO
            };
            let mut r = mul_scalar(&pms[0], q2);
            if r.is_empty() {
                r.resize(1, Mint::ZERO);
            }
            r[0] += q1;
            pre(&r, n)
        } else {
            if left >= n {
                return vec![];
            }
            let q1 = comp_rec(pms, q, n, left, j - 1);
            let q2 = comp_rec(pms, q, n, left + (1usize << (j - 1)), j - 1);
            let mut t = convolution(&pre(&pms[j - 1], n), &q2);
            t = pre(&t, n);
            let mut r = add_fps(&q1, &t);
            r = pre(&r, n);
            r
        }
    }

    let mut qpm = comp_rec(&pms, &q, n, 0, j);
    let mut r = qpm.clone();

    let mut pw_pr: Fps = vec![Mint::ONE];

    let mut dpm = diff(&pm);
    shrink(&mut dpm);

    let mut deg_dpm = 0usize;
    while deg_dpm < dpm.len() && dpm[deg_dpm] == Mint::ZERO {
        deg_dpm += 1;
    }
    let mut idpm: Fps = if dpm.is_empty() {
        vec![]
    } else {
        let dpm2 = shr(&dpm, deg_dpm);
        inv_series(&dpm2, n)
    };

    let mut d = m;
    for block in 1..=l {
        if d >= n {
            break;
        }
        // pw_Pr *= Pr
        pw_pr = pre(&convolution(&pw_pr, &pr), n - d);

        if dpm.is_empty() {
            // R += (pw_Pr * Q[block]) << d
            let term = mul_scalar(&pw_pr, q.get(block).copied().unwrap_or(Mint::ZERO));
            let term = shl(&term, d);
            r = add_fps(&r, &term);
        } else {
            // QPm = (QPm' / Pm')  (with shift for leading zeros)
            let mut dq = diff(&qpm);
            if deg_dpm > 0 {
                dq = shr(&dq, deg_dpm);
            }
            idpm.resize(n - d, Mint::ZERO);
            let mut next_qpm = convolution(&dq, &idpm);
            next_qpm = pre(&next_qpm, n - d);
            qpm = next_qpm;

            let mut term = convolution(&qpm, &pw_pr);
            term = pre(&term, n - d);
            term = mul_scalar(&term, c.finv(block as isize)); // / block!
            term = shl(&term, d);
            r = add_fps(&r, &term);
        }

        d += m;
    }

    r.resize(n, Mint::ZERO);
    r
}
