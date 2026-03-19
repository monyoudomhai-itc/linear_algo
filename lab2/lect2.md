# Chapter II: Determinant — Summary
**Course:** Linear Algebra | **Institute:** ITC | **Instructor:** PHAUK Sokkhey, Ph.D

---

## 1. Definition of the Determinant

### 2×2 Matrix
For a matrix $A = \begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{bmatrix}$, the determinant is:

$$\det(A) = |A| = a_{11}a_{22} - a_{21}a_{12}$$

**Examples:**
- $|A| = \begin{vmatrix} 2 & -3 \\ 1 & 2 \end{vmatrix} = 4 + 3 = 7$
- $|B| = \begin{vmatrix} 2 & 1 \\ 4 & 2 \end{vmatrix} = 4 - 4 = 0$
- $|C| = \begin{vmatrix} 0 & 3/2 \\ 2 & 4 \end{vmatrix} = 0 - 3 = -3$

---

## 2. Minor and Cofactor

### Definitions
Let $A = (a_{ij})_n$ and let $M_{ij}$ be the $(n-1) \times (n-1)$ submatrix obtained by **deleting the $i$-th row and $j$-th column** of $A$.

- **Minor** of $a_{ij}$: defined as $M_{ij}$ (the determinant of the submatrix)
- **Cofactor** of $a_{ij}$: defined as $C_{ij} = (-1)^{i+j} |M_{ij}|$

### Sign Pattern for Cofactors
The sign $(-1)^{i+j}$ follows a checkerboard pattern:

$$3 \times 3: \begin{bmatrix} + & - & + \\ - & + & - \\ + & - & + \end{bmatrix} \qquad 4 \times 4: \begin{bmatrix} + & - & + & - \\ - & + & - & + \\ + & - & + & - \\ - & + & - & + \end{bmatrix}$$

### Laplace Expansion (Cofactor Expansion)
The determinant of an $n \times n$ matrix can be computed by expanding along **any row or column**:

- **Row $i$ expansion:** $\det(A) = \sum_{j=1}^{n} a_{ij} C_{ij}$
- **Column $j$ expansion:** $\det(A) = \sum_{i=1}^{n} a_{ij} C_{ij}$

**Example:** For $A = \begin{bmatrix} 0 & 2 & 1 \\ 3 & -1 & 2 \\ 4 & 0 & 1 \end{bmatrix}$, expanding along row 1:

$$|A| = 0(-1) + 2(5) + 1(4) = 14$$

> The same result is obtained by expanding along any other row or column.

### Adjoint of a Matrix
The **adjoint** of $A$ is the transpose of the cofactor matrix:

$$\text{adj}(A) = (C_{ij})^t$$

**Key theorem:** $A \cdot \text{adj}(A) = \text{adj}(A) \cdot A = |A| \cdot I$

If $|A| \neq 0$, then: $A^{-1} = \dfrac{1}{|A|} \text{adj}(A)$

**Properties of adj:**
| Property | Formula |
|---|---|
| adj of inverse | $\text{adj}(A^{-1}) = (\text{adj}(A))^{-1}$ |
| adj of transpose | $\text{adj}(A^t) = (\text{adj}(A))^t$ |
| adj of product | $\text{adj}(AB) = \text{adj}(B)\,\text{adj}(A)$ |
| det of adj | $|\text{adj}(A)| = |A|^{n-1}$ |
| adj of scalar mult | $\text{adj}(\alpha A) = \alpha^{n-1}\,\text{adj}(A)$ |

---

## 3. Properties of Determinants

### Elementary Row Operations
If $B$ is obtained from $A$ by:

| Operation | Effect on det |
|---|---|
| Swap two rows | $|B| = -|A|$ |
| Multiply a row by scalar $\lambda$ | $|B| = \lambda|A|$ |
| Add $\lambda \times$ (row $j$) to row $i$ | $|B| = |A|$ (unchanged) |

### Special Cases (Theorem 2)
- If any **row or column is all zeros** → $|A| = 0$
- If **two rows (or columns) are equal** → $|A| = 0$
- If $A$ is a **triangular matrix** → $|A| = a_{11} \cdot a_{22} \cdots a_{nn}$ (product of diagonal)

### General Properties (Theorem 3)
For $A, B \in M_n(\mathbb{R})$ and scalar $\lambda$:

| Property | Formula |
|---|---|
| Transpose | $|A^t| = |A|$ |
| Product | $|AB| = |A| \cdot |B|$ |
| Inverse | $|A^{-1}| = \dfrac{1}{|A|}$ |
| Scalar multiple | $|\lambda A| = \lambda^n |A|$ |

**Example:** If $|A| = 7$ for a $3 \times 3$ matrix:
- $|2A| = 2^3 \cdot 7 = 56$
- $|3A^{-1}| = 3^3 \cdot \frac{1}{7} = \frac{27}{7}$
- $|(3A)^{-1}| = \frac{1}{|3A|} = \frac{1}{189}$

---

## 4. Cramer's Rule

### Setup
For the system $AX = b$ with $n$ equations and $n$ unknowns, where:

$$A = \begin{bmatrix} a_{11} & \cdots & a_{1n} \\ \vdots & & \vdots \\ a_{n1} & \cdots & a_{nn} \end{bmatrix}, \quad X = \begin{bmatrix} x_1 \\ \vdots \\ x_n \end{bmatrix}, \quad b = \begin{bmatrix} b_1 \\ \vdots \\ b_n \end{bmatrix}$$

### Theorem (Cramer's Rule)
The system has a **unique solution if and only if $\Delta = |A| \neq 0$**. The solution is:

$$x_1 = \frac{\Delta_1}{\Delta}, \quad x_2 = \frac{\Delta_2}{\Delta}, \quad \ldots, \quad x_n = \frac{\Delta_n}{\Delta}$$

where $\Delta_i$ = determinant of the matrix formed by **replacing the $i$-th column of $A$ with $b$**.

### Example
Solve: $-x + 2y - 3z = 1$, $2x + z = 0$, $3x - 4y + 4z = 2$

$$\Delta = |A| = 10 \neq 0 \implies \text{unique solution exists}$$

$$x = \frac{\Delta_1}{\Delta} = \frac{8}{10} = \frac{4}{5}$$

### Invertibility & Homogeneous Systems
A square matrix $A$ is **invertible** ⟺ $AX = 0$ has only the trivial solution ⟺ $|A| \neq 0$

The homogeneous system $AX = 0$ has a **nontrivial solution** if and only if $|A| = 0$.

---

## Quick Reference Summary

```
det(2×2): ad - bc

Cofactor:  C_ij = (-1)^(i+j) * M_ij

Laplace:   det(A) = Σ a_ij * C_ij  (any row/column)

Cramer:    x_i = Δ_i / Δ,  where Δ = |A| ≠ 0

Inverse:   A⁻¹ = (1/|A|) * adj(A)

Triangular: |A| = product of diagonal entries
```
