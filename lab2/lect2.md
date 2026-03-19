# Chapter II: Determinant — Summary

**Course:** Linear Algebra | **Institute:** ITC | **Instructor:** PHAUK Sokkhey, Ph.D

---

## 1. Definition of the Determinant

### 2×2 Matrix

For a matrix A with entries a₁₁, a₁₂, a₂₁, a₂₂:

```
det(A) = |A| = a₁₁·a₂₂ − a₂₁·a₁₂
```

**Examples:**

| Matrix | Calculation | Result |
|--------|-------------|--------|
| A = \|2  −3 / 1   2\| | 2(2) − 1(−3) = 4 + 3 | **7** |
| B = \|2   1 / 4   2\| | 2(2) − 4(1)  = 4 − 4 | **0** |
| C = \|0  3/2/ 2   4\| | 0(4) − 2(3/2) = 0 − 3 | **−3** |

---

## 2. Minor and Cofactor

### Definitions

Let A be an n×n matrix.

- **Minor Mᵢⱼ** — the determinant of the submatrix formed by **deleting row i and column j** from A.
- **Cofactor Cᵢⱼ** — defined as:

```
Cᵢⱼ = (−1)^(i+j) × Mᵢⱼ
```

### Sign Pattern (Checkerboard)

The sign `(−1)^(i+j)` follows this pattern:

```
3×3:          4×4:
+ − +         + − + −
− + −         − + − +
+ − +         + − + −
              − + − +
```

- If (i+j) is **even** → sign is **+**
- If (i+j) is **odd**  → sign is **−**

### Laplace Expansion

The determinant can be expanded along **any row or column**:

```
Row i expansion:    det(A) = Σⱼ aᵢⱼ · Cᵢⱼ
Column j expansion: det(A) = Σᵢ aᵢⱼ · Cᵢⱼ
```

> The result is the same regardless of which row or column you choose.

**Example** — for the 3×3 matrix A = [[0,2,1],[3,−1,2],[4,0,1]]:

```
Expanding along row 1:
|A| = 0·C₁₁ + 2·C₁₂ + 1·C₁₃
    = 0(−1) + 2(5) + 1(4)
    = 14
```

### Adjoint of A

The **adjoint** (adj) of A is the transpose of the cofactor matrix:

```
adj(A) = (Cᵢⱼ)ᵀ
```

Key relations:

```
A · adj(A) = adj(A) · A = |A| · I

A⁻¹ = (1 / |A|) · adj(A)     [only if |A| ≠ 0]
```

**Properties of adj:**

| Property | Formula |
|----------|---------|
| adj of inverse | adj(A⁻¹) = (adj(A))⁻¹ |
| adj of transpose | adj(Aᵀ) = (adj(A))ᵀ |
| adj of product | adj(AB) = adj(B) · adj(A) |
| det of adj | \|adj(A)\| = \|A\|ⁿ⁻¹ |
| adj of scalar multiple | adj(αA) = αⁿ⁻¹ · adj(A) |

---

## 3. Properties of Determinants

### Effect of Elementary Row Operations

| Operation | Effect on det(A) |
|-----------|-----------------|
| Swap two rows | det(B) = −det(A) |
| Multiply a row by scalar λ | det(B) = λ · det(A) |
| Add λ × (row j) to row i | det(B) = det(A) ← unchanged |

### Special Cases

| Condition | Result |
|-----------|--------|
| Any row or column is all zeros | \|A\| = 0 |
| Two rows (or columns) are equal | \|A\| = 0 |
| A is a triangular matrix | \|A\| = a₁₁ × a₂₂ × … × aₙₙ |

### General Properties

| Property | Formula |
|----------|---------|
| Transpose | \|Aᵀ\| = \|A\| |
| Product | \|AB\| = \|A\| · \|B\| |
| Inverse | \|A⁻¹\| = 1 / \|A\| |
| Scalar multiple | \|λA\| = λⁿ · \|A\| |
| det of adj | \|adj(A)\| = \|A\|ⁿ⁻¹ |

**Example** — if |A| = 7 for a 3×3 matrix:

```
|2A|       = 2³ × 7        = 56
|3A⁻¹|     = 3³ × (1/7)   = 27/7
|(3A)⁻¹|   = 1 / |3A|     = 1/189
```

---

## 4. Cramer's Rule

### System AX = b

For a system of n equations with n unknowns:

```
a₁₁x₁ + a₁₂x₂ + … + a₁ₙxₙ = b₁
a₂₁x₁ + a₂₂x₂ + … + a₂ₙxₙ = b₂
⋮
aₙ₁x₁ + aₙ₂x₂ + … + aₙₙxₙ = bₙ
```

### Formula

```
Δ  = |A|
Δᵢ = |A with column i replaced by b|

xᵢ = Δᵢ / Δ     [requires Δ ≠ 0]
```

> The system has a **unique solution** if and only if Δ ≠ 0.

### Example

System: −x + 2y − 3z = 1,  2x + z = 0,  3x − 4y + 4z = 2

```
Δ = |A| = 10  →  unique solution exists

x = Δ₁ / Δ = 8/10 = 4/5
```

### Invertibility Equivalence

For a square matrix A, the following are all equivalent:

```
(1) A is invertible
(2) AX = 0 has only the trivial solution (X = 0)
(3) |A| ≠ 0
```

> AX = 0 has a **nontrivial** solution ⟺ |A| = 0

---

## Quick Reference

```
2×2 det:      ad − bc

Cofactor:     Cᵢⱼ = (−1)^(i+j) × Mᵢⱼ

Laplace:      det(A) = Σ aᵢⱼ · Cᵢⱼ  (any row or column)

Inverse:      A⁻¹ = (1/|A|) · adj(A)

Cramer:       xᵢ = Δᵢ / Δ,   Δ = |A| ≠ 0

Triangular:   |A| = product of diagonal entries
```