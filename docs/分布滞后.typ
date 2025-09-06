#import "@local/modern-cug-report:0.1.2": *
#show: doc => template(doc, size: 15pt, footer: "CUG水文气象学2025", header: "")

// #let ell = sym.ell

= 1 非线性滞后 

== 1.1 交叉基展开

$ f(x, tau) = sum_{v,ell} beta_{v ell} thin s_v(x) h_ell(tau) $

- `basisvar`: $[n, v_x]$，`basislag`: $[L, v_ell]$

- $v_x$: 非线性影响基函数个数；
- $v_ell$: 滞后基函数个数；

#box-red[
  如果是bspline，需要指定knots节点，以及基函数的个数。
]

== 1.2 考虑滞后时间

$
  eta_t & = sum_{tau=0}^L f(x_{t - tau}, tau), tau in [0, 1, ..., L] \
  eta_t & = sum_{tau=0}^L sum_{v=0}^v_x sum_{ell=0}^v_ell beta_{v_ell} s_v (x_{t - tau}) h_ell (tau)
$

```R
x = basisvar[, v] # 第i个基函数
mat <- as.matrix(Lag(basisvar[, v], seqlag(lag), group = group))
# mat: {X[n, t], X[n, t-1], ..., X[n, t-L]}, [n, L+1]

for (l in seq(length = vl)) {
  ck <- basislag[, l] # [L+1]
  crossbasis[, vl * (v - 1) + l] <- mat %*% ck # [n, L+1] %*% [L+1, 1] = [n]
}
```

= 2 矩阵形式

上式写成，矩阵形式

$ S(x) = [s_1 (x), ..., s_{v x} (x)]^T, in RR^[n, v_x] $

// $ H = [H_1 (tau), ..., H_{v ell}] in RR^[L+1, v_ell] $

$ H(tau) = [h_1 (tau), ..., h_{v ell}] in RR^{v_ell}, tau in [0, ..., L] $

$
  B in R^{v_x, v_ell}, beta = "vec"(B) in RR^{v_x times v_ell}
$

$ f(x, tau) = S(x)^T B H(tau) =[ H(tau) times.circle S(x) ]^T beta $

$
  eta_t = sum_{tau=0}^{L} [ H(tau) times.circle S(x) ]^T beta
$
