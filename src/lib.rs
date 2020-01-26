extern crate num_bigint as bigint;
extern crate num_integer as integer;
extern crate num_traits as traits;
extern crate rand;
extern crate rayon;

use {
    bigint::{BigUint, RandBigInt, ToBigUint},
    rayon::prelude::*,
    std::iter::repeat_with,
    traits::{One, Zero},
};

macro_rules! biguint {
    ($e:expr) => {
        ($e).to_biguint().unwrap()
    };
}

fn decompose(n: &BigUint) -> (BigUint, BigUint) {
    let one = One::one();
    let ref two = biguint!(2);
    let ref n_minus_one: BigUint = n - 1u8;
    let mut d: BigUint = n_minus_one.clone();
    let mut r: BigUint = Zero::zero();

    while &d % two == one {
        d /= two;
        r += 1u8;
    }

    (d, r)
}

fn __miller_rabin(a: &BigUint, n: &BigUint, d: &BigUint, r: &BigUint) -> bool {
    let n_minus_one = n - 1u8;
    let mut x = a.modpow(d, n);
    let mut count: BigUint = One::one();
    let ref two = biguint!(2);

    if x == One::one() || x == n_minus_one {
        return false;
    }
    while &count < r {
        x = x.modpow(two, n);
        if x == n_minus_one {
            return false;
        }
        count += 1u8;
    }

    true
}

/// Test whether an integer `a` is a witness for the compositeness of `n`.
///
/// NOTE: This function fails if `a < 2` or `n < 3`.
///
/// # Examples
///
/// ```
/// extern crate rand;
/// use rand::distributions::{Distribution, Uniform};
/// use miller_rabin::is_witness;
///
/// let n: u64 = 27;
/// let dist = Uniform::new(2, n);
/// let mut rng = rand::thread_rng();
/// // A random integer in [2..n]
/// let a: u64 = dist.sample(&mut rng);
/// assert!(is_witness(&a, &n).unwrap());
/// ```
pub fn is_witness<T: ToBigUint>(a: &T, n: &T) -> Option<bool> {
    let (ref a, ref n) = (biguint!(a), biguint!(n));

    if a < &biguint!(2) || n < &biguint!(3) {
        return None;
    }

    let (ref d, ref r) = decompose(n);
    Some(__miller_rabin(a, n, d, r))
}

/// Test whether an integer `n` is likely prime using the Miller-Rabin primality test.
///
/// # Examples
///
/// ```
/// use miller_rabin::is_prime;
///
/// // Mersenne Prime (2^31 - 1)
/// let n: u64 = 0x7FFF_FFFF;
/// // Try the miller-rabin test 100 times in parallel
/// // (or, if `n` is less than `u64::max()`, test only the numbers known to be sufficient).
/// // In general, the algorithm should fail at a rate of at most `4^{-k}`.
/// assert!(is_prime(&n, 100));
/// ```
pub fn is_prime<T: ToBigUint>(n: &T, k: usize) -> bool {
    let ref n = biguint!(n);
    let (ref d, ref r) = decompose(n);

    if n == &Zero::zero() {
        return false;
    } else if n < &biguint!(3) {
        return true;
    } else if n < &biguint!(0xFFFF_FFFF_FFFF_FFFFu64) {
        let samples: Vec<u8> = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
        return samples
            .par_iter()
            .find_any(|&&a| __miller_rabin(&biguint!(a), n, d, r))
            .is_none();
    }

    let mut rng = rand::thread_rng();
    let samples: Vec<BigUint> = repeat_with(|| rng.gen_biguint(n.bits()))
        .filter(|m| m < &n)
        .take(k)
        .collect();

    samples
        .par_iter()
        .find_any(|&a| __miller_rabin(a, n, d, r))
        .is_none()
}

#[cfg(test)]
mod tests {
    const K: usize = 16;

    use super::*;
    use std::io;

    #[test]
    fn test_prime() -> io::Result<()> {
        let prime: u64 = 0x7FFF_FFFF;
        assert!(is_prime(&prime, K));
        Ok(())
    }

    #[test]
    fn test_composite() -> io::Result<()> {
        let composite: u64 = 0x7FFF_FFFE;
        assert!(!is_prime(&composite, K));
        Ok(())
    }

    #[test]
    fn test_big_mersenne_prime() -> io::Result<()> {
        let prime: BigUint =
            BigUint::parse_bytes(b"170141183460469231731687303715884105727", 10).unwrap();
        assert!(is_prime(&prime, K));
        Ok(())
    }

    #[test]
    fn test_big_wagstaff_prime() -> io::Result<()> {
        let prime: BigUint =
            BigUint::parse_bytes(b"56713727820156410577229101238628035243", 10).unwrap();
        assert!(is_prime(&prime, K));
        Ok(())
    }

    #[test]
    fn test_big_composite() -> io::Result<()> {
        let prime: BigUint =
            BigUint::parse_bytes(b"170141183460469231731687303715884105725", 10).unwrap();
        assert!(!is_prime(&prime, K));
        Ok(())
    }
}
