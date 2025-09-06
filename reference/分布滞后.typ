#import "@local/modern-cug-report:0.1.2": *
#show: doc => template(doc, size: 15pt, footer: "CUG水文气象学2025", header: "")

#set par(leading: 1.24em, spacing: 1.24em)

// #let ell = sym.ell
= 1 *非线性滞后 *

== 1.1 交叉基展开

$ f(x, tau) = sum_{v,ell} beta_{v ell} thin s_v (x) h_ell(tau) $

- `basisvar`: $[n, v_x]$，`basislag`: $[L, v_ell]$

- $v_x$: 非线性影响基函数个数；$v_ell$: 滞后基函数个数。

#box-red[
  如果是bspline，需要指定knots节点，以及基函数的个数。
]

#v(-0.4em)
*bspline为何能实现压缩模型参数个数？* 类似于曲线拟合，只需要简单的几个节点，便可描述非线性过程。

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

= 2 *矩阵形式*

上式写成，矩阵形式

$ S(x) = [s_1 (x), ..., s_{v x} (x)]^T, in RR^{n, v_x} $

// $ H = [H_1 (tau), ..., H_{v ell}] in RR^[L+1, v_ell] $

$ H(tau) = [h_1 (tau), ..., h_{v ell}] in RR^{v_ell}, tau in [0, ..., L] $

$ B in R^{v_x, v_ell}, beta = "vec"(B) in RR^{v_x times v_ell} $

$ f(x, tau) = S(x)^T B H(tau) =[ H(tau) times.circle S(x) ]^T beta $

$ eta_t = sum_{tau=0}^{L} [ H(tau) times.circle S(x) ]^T beta $


= 3 *B样条滞后基的构造（log 间距结点）*

- 选择滞后变换
  $u = g(tau)$（常用 $g(tau) = log(tau + tau_0), tau_0 > 0$）。
- 在 u 轴给定结点序列
  $kappa = { kappa_0, ..., kappa_K }$（对数间距），样条次数 q（如 q = 3）。
- 定义
  $h_ell(tau) = N_(ell, q)( g(tau) ; kappa ), quad ell = 1, ..., v_ell = K + q$，
  其中 $N_(ell, q)$ 为对应结点与次数的 B 样条基。


$ f(x, tau) = sum_(ell=1)^(v_ell) beta_ell thin x thin h_ell(tau) $

$ eta_t = sum_(tau=0)^(L) f(x_(t - tau), tau), tau in [0, 1, ..., L] $

$
  eta_t = sum_{tau=0}^L sum_{ell=0}^ell beta_{ell} x_{t - tau} h_ell (tau)
  = sum_(ell=1)^(v_ell) beta_ell thin sum_(tau=0)^(L) h_ell(tau) thin x_(t - tau)
$

*另一种写法*

$
  s_j (t) = sum_{tau = 0}^L x_{t - tau} h_j (tau), quad
  eta_t = sum_{j=1}^{v ell} beta_j s_j (t)
$

$ X = [x_{t}, ..., x_{t-l}, ..., x_{t - L}]in RR^{n, L+1}, $

$ H_j = [h_j (0), ..., h_j (l), ..., h_j (L)]^T in RR^{L+1,1}, H = [H_1, ..., H_{v ell}] in RR^{L+1, v_ell} $ 

$ S = X H in RR^{n, v_ell}, quad eta = S beta in RR^n, beta = [beta_1, ..., beta_{v ell}]^T $

其中$H_j$为bspline生成的以logknots为输入的基函数。

// #box-red[
//   直觉：小滞后区（秒–分钟）用密集结点刻画陡峭“早期”峰，大滞后区（小时–天）用稀疏结点刻画“长尾”衰减，在相同参数数目下取得更均衡的时域分辨率。
// ]

= 5 *从时域核到频域传递函数*

先定义时域加载核
$rho(tau) := sum_(ell=1)^(v_ell) beta_ell thin h_ell(tau)$。

其离散傅里叶变换给出传递函数
$R(omega) = sum_(tau=0)^(L_b) rho(tau) e^(- i omega tau)$。

- 幅值：$abs(R(omega))$（频率依赖幅值比）
- 相位：$arg R(omega)$
- 步响应（近似静态 BE）：$S(T) = sum_(tau=0)^T rho(tau)$

#show "---": doc => { pagebreak() }
