use std::mem;
use std::path::Path;

use anyhow::Result;
use bio::stats::LogProb;
use itertools_num::linspace;
use statrs::function::beta::ln_beta;

use crate::kallisto::KallistoQuant;
use crate::prior::Prior;

pub(crate) fn sample_expression(
    preprocessing: &Path,
    kallisto_quant: &Path,
    epsilon: LogProb,
    prior: &Prior,
) -> Result<()> {
    let quant = KallistoQuant::new(kallisto_quant)?;

    

    let feature_ids = quant.feature_ids()?;
    let estimates = Estimates::new(&quant)?;

    for i in 0..feature_ids.len() {
        let d_ij = estimates.means[i];
        if let Some(t_ij) = estimates.dispersions[i] {
            let (left_window, right_window) = window(d_ij);

            let mut likelihoods = Vec::new();

            let mut likelihood_mu_ik = |mu_ik| {
                let mut max_prob = LogProb::ln_zero();
                for theta_i in prior.window() {
                    let prob = likelihood_mu_ik_theta_i(d_ij, mu_ik, t_ij, theta_i, epsilon);

                    if prob > max_prob {
                        max_prob = prob;
                    }

                    likelihoods.push(Likelihood {
                        mu_ik: mu_ik,
                        theta_i: theta_i,
                        prob,
                    });
                }

                max_prob
            };

            for mu_ik in left_window {
                if likelihood_mu_ik(mu_ik) < epsilon {
                    break;
                }
            }

            for mu_ik in right_window {
                if likelihood_mu_ik(mu_ik) < epsilon {
                    break;
                }
            }
        } else {
        }
    }

    Ok(())
}

struct Likelihood {
    mu_ik: f64,
    theta_i: f64,
    prob: LogProb,
}

fn likelihood_mu_ik_theta_i(
    d_ij: f64,
    mu_ik: f64,
    t_ij: f64,
    theta_i: f64,
    epsilon: LogProb,
) -> LogProb {
    // TODO determine whether this is the best window given that we also have access to mu_ik here.
    let (left_window, right_window) = window(d_ij);

    let prob = |x| LogProb((neg_binom(d_ij, x, t_ij) * neg_binom(x, mu_ik, theta_i)).ln());
    let is_greater_epsilon = |prob: &LogProb| *prob >= epsilon;

    LogProb::ln_sum_exp(
        &left_window
            .map(&prob)
            .take_while(&is_greater_epsilon)
            .chain(right_window.map(&prob).take_while(&is_greater_epsilon))
            .collect::<Vec<_>>(),
    )
}

fn window(d_ij: f64) -> (impl Iterator<Item = f64>, impl Iterator<Item = f64>) {
    (
        linspace(d_ij, 5.0 * d_ij, 100),
        linspace(d_ij / 5.0, d_ij, 100).rev(),
    )
}

fn neg_binom(x: f64, mu: f64, theta: f64) -> f64 {
    let n = 1.0 / theta;
    let p = n / (n + mu);
    let mut p1 = if n > 0.0 { n * p.ln() } else { 0.0 };
    let mut p2 = if x > 0.0 { x * (1.0 - p).ln() } else { 0.0 };
    let b = ln_beta(x + 1.0, n);

    if p1 < p2 {
        mem::swap(&mut p1, &mut p2);
    }
    (p1 - b + p2).exp() / (x + n)
    // (p1 - b + p2) - (x + n).ln() // TODO is this the correct form for returning LogProb?
}

// def neg_binom(x, mu, theta):
//     """ Own implementation of the negative negative binomial distribution using betaln
//     """
//     n = 1. / theta
//     p = n / (n + mu)
//     p1 = n*np.log(p) if n > 0 else 0
//     p2 = x*np.log(1-p) if x > 0 else 0
//     b = betaln(x + 1, n)
//     if (p1 < p2):
//         return exp(p2 - b + p1) / (x+n)
//     else:
//         return exp(p1 - b + p2) / (x+n)
