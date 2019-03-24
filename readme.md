# Assignment: Interpolation

数值计算方法第一次作业。

## 实现功能

设在区间 $[0, 6]$ 上有函数 $f(x)=e^{-2x}$。把区间 $n$ 等分，给定各结点上的函数值。利用三种不同的插值算法：

1. 全区间上的 $n​$ 阶 Lagrange 插值多项式；
2. 分段线性插值式；
3. 三次样条插值式（采用第 1 种边界条件）；

对 $n=20$ 和 $n=30$，分别计算

$$
x_{k-\frac{1}{2}}=\frac{6(k-\frac{1}{2})}{n} (1\le k\le n)
$$

除此以外，本人还采用了 $n$ 阶 Newton 插值公式进行计算。



## 算法分析

### Lagrange 插值多项式

根据 $n$ 次 Lagrange 插值公式：
$$
\varphi_n(x_j)=\sum _{i=0}^n \frac{w_n(x_j)y_i}{(x_j-x_i)w'_n(x_i)}
$$

定义
$$
\begin{aligned}
frac(n)&=n\cdot(n-1)\cdots 1\\
frac2(n)&=(n-0.5)\cdot(n-1.5)\cdots 0.5\\
\end{aligned}
$$

我们可以在 $O(n)$ 时间复杂度内完成上述两个函数的全部计算，再带入上式中，即可得到对应的结果。总体算法复杂度为 $O(n^2)$。

### Newton 插值多项式

Newton 插值采用递推的思路简化了计算，因此不需要做额外的优化就可以实现 $O(n^2)$ 复杂度的求值。做法完全同课本。总体算法复杂度为 $O(n^2)$。

### 分段线性插值

直接两两取平均即可。总体算法复杂度为 $O(n)$。

### 样条插值

首先计算出 $\alpha_i, \beta_i$ 的结果：
``` c++
a[0] = 0; b[0] = -4;
a[N] = 1; b[N] = -4 * target[N];
long double coefb = 1.5 / unit;
for (int i = 1; i < N; ++i) {
    a[i] = 0.5;
    b[i] = (source[i + 1] - source[i - 1]) * coefb;
}
```

再依次求出 $A_i, B_i$：
``` c++
A[0] = 0; B[0] = -2;
for (int i = 1; i <= N; ++i) {
long double deno = 2 + (1 - a[i]) * A[i - 1];
    A[i] = -a[i] / deno;
    B[i] = (b[i] - (1 - a[i]) * B[i - 1]) / deno;
}
```

之后反向求出 $m_i$，并利用
$$
\varphi_i=\frac{y_i+y_{i+1}}{2}+\frac{d(m_i-m_{i+1})}{8}
$$

得到最终的值。总体算法复杂度为 $O(n^2)$。

## 实验结果

环境：Windows 10，Sublime Text 3，TDM-GCC 6.8.1 64-bit

可以参考下面的三个文件中的数据：

1. `out_20.csv`: $n=20$ 时四种方法在每个测试点处的绝对误差；
2. `out_30.csv`: $n=30$ 时四种方法在每个测试点处的绝对误差；
3. `errors.csv`: $n=10, 20, \cdots, 90$ 时四种方法的最大误差比较。

## 数据分析

### out_20.csv

分析：线性插值误差最大，样条插值其次，Newton/Lagrange 插值最小且它们几乎相等（因为它们本来就应该相等）。可见数据量较小时 Newton/Lagrange 插值体现出了非常大的优势。

### out_30.csv

分析：线性插值略微优化，样条插值显著优化，Newton/Lagrange 插值显著优化。

### errors.csv

分析：数据量继续增大时，Newton/Lagrange 插值的误差开始增大，直至线性插值和样条插值；线性插值优化不明显，样条插值的效果得到了不错的提升。

