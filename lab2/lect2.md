# Chapter II: Determinant — Full Study Notes

**Course:** Linear Algebra | **Institute:** ITC | **Instructor:** PHAUK Sokkhey, Ph.D

---

## Table of Contents

1. [Definition of the Determinant](#1-definition-of-the-determinant)
2. [Minor and Cofactor](#2-minor-and-cofactor)
3. [Properties of Determinants](#3-properties-of-determinants)
4. [Cramer's Rule](#4-cramers-rule)

---

## 1. Definition of the Determinant

### 1.1 — 2×2 Matrix

For a 2×2 matrix A:

```
    | a₁₁  a₁₂ |
A = |           |
    | a₂₁  a₂₂ |
```

The determinant is defined as:

```
det(A) = |A| = a₁₁·a₂₂ − a₂₁·a₁₂
```

> **In plain words:** multiply the top-left by bottom-right, then subtract (bottom-left × top-right).

---

### 1.2 — Worked Examples (2×2)

**Example 1:**
```
    | 2  -3 |
A = |       |     →   |A| = 2(2) − 1(−3) = 4 + 3 = 7
    | 1   2 |
```

**Example 2:**
```
    | 2  1 |
B = |      |     →   |B| = 2(2) − 4(1) = 4 − 4 = 0
    | 4  2 |
```
> Note: |B| = 0 means the matrix is **singular** (not invertible). The two rows are proportional.

**Example 3:**
```
    | 0   3/2 |
C = |         |     →   |C| = 0(4) − 2(3/2) = 0 − 3 = −3
    | 2   4   |
```

---

## 2. Minor and Cofactor

### 2.1 — Definitions

Let A = (aᵢⱼ) be an n×n matrix.

**Minor Mᵢⱼ:**
- Delete **row i** and **column j** from A.
- The minor Mᵢⱼ is the **determinant** of the remaining (n−1)×(n−1) submatrix.

**Cofactor Cᵢⱼ:**
- The cofactor is the minor with a sign attached:

```
Cᵢⱼ = (−1)^(i+j) × Mᵢⱼ
```

- If (i+j) is **even** → Cᵢⱼ = +Mᵢⱼ
- If (i+j) is **odd**  → Cᵢⱼ = −Mᵢⱼ

---

### 2.2 — Sign Pattern (Checkerboard)

The signs `(−1)^(i+j)` form a checkerboard pattern across the matrix:

```
For a 3×3 matrix:         For a 4×4 matrix:

  + − +                     + − + −
  − + −                     − + − +
  + − +                     + − + −
                            − + − +
```

**General rule for any size n×n:**
```
Position (i,j) has sign (+) if i+j is even
Position (i,j) has sign (−) if i+j is odd
```

---

### 2.3 — Finding Minors and Cofactors (Example)

Given:
```
    | 0   2   1 |
A = | 3  -1   2 |
    | 4   0   1 |
```

**Finding M₁₁** (delete row 1, column 1):
```
Remaining matrix:
    | -1   2 |
    |  0   1 |

M₁₁ = (−1)(1) − (0)(2) = −1 − 0 = −1
C₁₁ = (−1)^(1+1) × (−1) = (+1)(−1) = −1
```

**Finding M₁₂** (delete row 1, column 2):
```
Remaining matrix:
    | 3   2 |
    | 4   1 |

M₁₂ = (3)(1) − (4)(2) = 3 − 8 = −5
C₁₂ = (−1)^(1+2) × (−5) = (−1)(−5) = 5
```

**Finding M₁₃** (delete row 1, column 3):
```
Remaining matrix:
    | 3  -1 |
    | 4   0 |

M₁₃ = (3)(0) − (4)(−1) = 0 + 4 = 4
C₁₃ = (−1)^(1+3) × 4 = (+1)(4) = 4
```

**All minors and cofactors for matrix A above:**
```
Minors:               Cofactors:
M₁₁ = −1             C₁₁ = −1
M₁₂ = −5             C₁₂ =  5
M₁₃ =  4             C₁₃ =  4
M₂₁ =  2             C₂₁ = −2
M₂₂ = −4             C₂₂ = −4
M₂₃ = −8             C₂₃ =  8
M₃₁ =  5             C₃₁ =  5
M₃₂ = −3             C₃₂ =  3
M₃₃ = −6             C₃₃ = −6
```

---

### 2.4 — Laplace Expansion (Cofactor Expansion)

For any n×n matrix, the determinant can be computed by expanding along **any row or column**.

**Row i expansion:**
```
det(A) = aᵢ₁·Cᵢ₁ + aᵢ₂·Cᵢ₂ + … + aᵢₙ·Cᵢₙ
```

**Column j expansion:**
```
det(A) = a₁ⱼ·C₁ⱼ + a₂ⱼ·C₂ⱼ + … + aₙⱼ·Cₙⱼ
```

> **Important:** The result is **always the same** no matter which row or column you pick. Choose the one with the most zeros to save work!

---

### 2.5 — Computing det of a 3×3 Matrix (Example)

Using the same matrix:
```
    | 0   2   1 |
A = | 3  -1   2 |
    | 4   0   1 |
```

**Expand along Row 1:**
```
|A| = a₁₁·C₁₁ + a₁₂·C₁₂ + a₁₃·C₁₃
    = 0·(−1) + 2·(5) + 1·(4)
    = 0 + 10 + 4
    = 14
```

**Verify by expanding along Row 2:**
```
|A| = a₂₁·C₂₁ + a₂₂·C₂₂ + a₂₃·C₂₃
    = 3·(−2) + (−1)·(−4) + 2·(8)
    = −6 + 4 + 16
    = 14  ✓
```

**Verify by expanding along Column 1:**
```
|A| = a₁₁·C₁₁ + a₂₁·C₂₁ + a₃₁·C₃₁
    = 0·(−1) + 3·(−2) + 4·(5)
    = 0 − 6 + 20
    = 14  ✓
```

---

### 2.6 — Adjoint of a Matrix

**Definition:**
The adjoint (or adjugate) of A is the **transpose of the cofactor matrix**:

```
adj(A) = (Cᵢⱼ)ᵀ
```

In expanded form, the cofactor matrix is:

```
           | C₁₁  C₁₂  …  C₁ₙ |ᵀ
adj(A)  =  | C₂₁  C₂₂  …  C₂ₙ |
           |  ⋮    ⋮       ⋮   |
           | Cₙ₁  Cₙ₂  …  Cₙₙ |
```

> Note: The (i,j) entry of adj(A) is Cⱼᵢ (transposed — row and column are swapped).

**Key theorem:**
```
A · adj(A) = adj(A) · A = |A| · I
```

**Inverse formula using adjoint:**
```
A⁻¹ = (1 / |A|) · adj(A)       [only valid when |A| ≠ 0]
```

---

### 2.7 — Adjoint Example

Given:
```
     | −1   3   2 |
A =  |  0  −2   1 |
     |  1   0  −2 |
```

Step 1 — compute all cofactors:
```
C₁₁ = (−1)^2 · |(−2)(−2)−(0)(1)|  = +(4−0)  =  4
C₁₂ = (−1)^3 · |(0)(−2)−(1)(1)|   = −(0−1)  =  1
C₁₃ = (−1)^4 · |(0)(0)−(1)(−2)|   = +(0+2)  =  2
C₂₁ = (−1)^3 · |(3)(−2)−(0)(2)|   = −(−6−0) =  6
C₂₂ = (−1)^4 · |(−1)(−2)−(1)(2)|  = +(2−2)  =  0
C₂₃ = (−1)^5 · |(−1)(0)−(1)(3)|   = −(0−3)  =  3
C₃₁ = (−1)^4 · |(3)(1)−(−2)(2)|   = +(3+4)  =  7
C₃₂ = (−1)^5 · |(−1)(1)−(0)(2)|   = −(−1−0) =  1
C₃₃ = (−1)^6 · |(−1)(−2)−(0)(3)|  = +(2−0)  =  2
```

Step 2 — cofactor matrix:
```
    | 4  1  2 |
    | 6  0  3 |
    | 7  1  2 |
```

Step 3 — transpose to get adj(A):
```
           | 4  6  7 |
adj(A)  =  | 1  0  1 |
           | 2  3  2 |
```

---

### 2.8 — Properties of adj

| Property | Formula |
|----------|---------|
| adj of inverse | adj(A⁻¹) = (adj(A))⁻¹ |
| adj of transpose | adj(Aᵀ) = (adj(A))ᵀ |
| adj of product | adj(AB) = adj(B) · adj(A) |
| determinant of adj | \|adj(A)\| = \|A\|ⁿ⁻¹ |
| adj of scalar multiple | adj(αA) = αⁿ⁻¹ · adj(A) |

---

## 3. Properties of Determinants

### 3.1 — Effect of Elementary Row Operations

These are the three types of row operations and how each affects the determinant:

**Operation 1 — Swap two rows:**
```
If B is obtained from A by swapping any two rows:
    det(B) = −det(A)
```
> Each row swap flips the sign of the determinant.

**Operation 2 — Multiply a row by a scalar λ:**
```
If B is obtained from A by multiplying one row by λ:
    det(B) = λ · det(A)
```
> Factoring out λ from a row scales the determinant by λ.

**Operation 3 — Add a multiple of one row to another:**
```
If B is obtained from A by: (row i) ← (row i) + λ·(row j),  i ≠ j:
    det(B) = det(A)
```
> This operation does NOT change the determinant.

---

### 3.2 — Worked Example (Row Operations)

Starting with |A| = 11:
```
    | 2  -3 |
A = |       |     →   |A| = 2(4)−1(−3) = 11
    | 1   4 |
```

**Swap rows → sign flips:**
```
    | 1   4 |
B = |       |     →   |B| = −11
    | 2  -3 |
```

**Multiply row 1 by 1/2:**
```
    | 1  -4 |
B = |       |     →   |B| = (1/2) · 11 = 11/2 = 5.5  (not 11!)
    | 1   4 |
```

**Add −2×(row 1) to row 2:**
```
    | 1  -3 |
B = |       |     →   |B| = |A| = 2   (unchanged)
    | 0   2 |
```

---

### 3.3 — Using Row Reduction to Find det

**Strategy:** Use row operations to reduce A to triangular form, then multiply diagonal entries.

**Example:** Find det of:
```
    |  0  -7  14 |
A = |  1   2  -2 |
    |  0   3  -8 |
```

Step 1 — Swap row 1 and row 2 (sign flips):
```
    |  1   2  -2 |
= − |  0  -7  14 |
    |  0   3  -8 |
```

Step 2 — Factor −7 out of row 2:
```
    |  1   2  -2 |
= 7 |  0   1  -2 |
    |  0   3  -8 |
```

Step 3 — Add −3×(row 2) to row 3:
```
    |  1   2  -2 |
= 7 |  0   1  -2 |   ← now triangular!
    |  0   0  -2 |
```

Step 4 — det of triangular matrix = product of diagonal:
```
|A| = 7 × (1)(1)(−2) = 7 × (−2) = −14
```

---

### 3.4 — Special Cases

| Condition | Result | Reason |
|-----------|--------|--------|
| A row (or column) is all zeros | \|A\| = 0 | Expanding along that row gives all zero terms |
| Two rows (or columns) are identical | \|A\| = 0 | Swapping them changes sign but nothing changes, so must be 0 |
| A is upper or lower triangular | \|A\| = a₁₁ × a₂₂ × … × aₙₙ | Off-diagonal terms vanish in expansion |
| A is diagonal | \|A\| = a₁₁ × a₂₂ × … × aₙₙ | Special case of triangular |
| A is the identity matrix | \|I\| = 1 | Product of diagonal = 1×1×…×1 |

---

### 3.5 — General Algebraic Properties

| Property | Formula | Notes |
|----------|---------|-------|
| Transpose | \|Aᵀ\| = \|A\| | Rows and columns are interchangeable |
| Product | \|AB\| = \|A\| · \|B\| | Multiplicative property |
| Inverse | \|A⁻¹\| = 1 / \|A\| | Follows from \|A·A⁻¹\| = \|I\| = 1 |
| Scalar multiple (n×n) | \|λA\| = λⁿ · \|A\| | λ is factored from each of n rows |
| Powers | \|Aᵏ\| = \|A\|ᵏ | Follows from the product rule |
| det of adj | \|adj(A)\| = \|A\|ⁿ⁻¹ | — |

---

### 3.6 — Worked Example (General Properties)

Suppose A is a 3×3 matrix with **|A| = 7**. Compute:

**(a) |2A|**
```
|2A| = 2³ · |A| = 8 × 7 = 56
```

**(b) |3A⁻¹|**
```
|3A⁻¹| = 3³ · |A⁻¹| = 27 · (1/7) = 27/7
```

**(c) |(3A)⁻¹|**
```
|(3A)⁻¹| = 1 / |3A| = 1 / (3³ · 7) = 1 / (27 × 7) = 1/189
```

**(d) det of Aᵀ with columns rearranged**

If two columns of A are swapped, the determinant changes sign:
```
|A with columns 2 and 3 swapped| = −|A| = −7
```

---

### 3.7 — Determinant of a Matrix Product (Example)

Given:
```
    |  1  -2   2 |           |  2   0   1 |
A = |  0   3   2 |     B =   |  0  -1  -2 |
    |  1   0   1 |           |  3   1  -2 |
```

Compute separately:
```
|A| = −7      |B| = 11
```

Product AB:
```
     |  8   4   1  |
AB = |  6  -1  -10 |
     |  5   1  -1  |

|AB| = −77
```

Verify: |AB| = |A| · |B| = (−7)(11) = −77 ✓

---

## 4. Cramer's Rule

### 4.1 — The System AX = b

Consider a system of n linear equations with n unknowns:

```
a₁₁x₁ + a₁₂x₂ + a₁₃x₃ + … + a₁ₙxₙ = b₁
a₂₁x₁ + a₂₂x₂ + a₂₃x₃ + … + a₂ₙxₙ = b₂
a₃₁x₁ + a₃₂x₂ + a₃₃x₃ + … + a₃ₙxₙ = b₃
⋮
aₙ₁x₁ + aₙ₂x₂ + aₙ₃x₃ + … + aₙₙxₙ = bₙ
```

In matrix form:
```
A · X = b

where:
A = coefficient matrix (n×n)
X = [x₁, x₂, …, xₙ]ᵀ  (unknowns)
b = [b₁, b₂, …, bₙ]ᵀ  (right-hand side)
```

---

### 4.2 — Cramer's Rule Formula

**Step 1:** Compute Δ = |A| (determinant of coefficient matrix).

**Step 2:** For each unknown xᵢ, form matrix Aᵢ by replacing the **i-th column** of A with b.

**Step 3:** Compute Δᵢ = |Aᵢ|.

**Step 4:** The solution is:
```
x₁ = Δ₁/Δ,   x₂ = Δ₂/Δ,   …,   xₙ = Δₙ/Δ
```

> **Condition:** The system has a **unique solution** if and only if Δ = |A| ≠ 0.
> If Δ = 0, Cramer's Rule cannot be applied.

---

### 4.3 — Worked Example (Cramer's Rule)

**Solve for x only** in the system:
```
−x + 2y − 3z = 1
 2x      +  z = 0
 3x − 4y + 4z = 2
```

**Coefficient matrix A and vector b:**
```
     | −1   2  −3 |       | 1 |
A =  |  2   0   1 |   b = | 0 |
     |  3  −4   4 |       | 2 |
```

**Step 1 — Compute Δ = |A|:**
```
Δ = |A| = −1[(0)(4)−(1)(−4)] − 2[(2)(4)−(1)(3)] + (−3)[(2)(−4)−(0)(3)]
        = −1[0+4] − 2[8−3] − 3[−8−0]
        = −4 − 10 + 24
        = 10
```

**Step 2 — Form A₁** (replace column 1 with b):
```
      |  1   2  −3 |
A₁ =  |  0   0   1 |
      |  2  −4   4 |
```

**Step 3 — Compute Δ₁ = |A₁|:**
```
Δ₁ = 1[(0)(4)−(1)(−4)] − 2[(0)(4)−(1)(2)] + (−3)[(0)(−4)−(0)(2)]
   = 1[0+4] − 2[0−2] − 3[0−0]
   = 4 + 4 + 0
   = 8
```

**Step 4 — Solve:**
```
x = Δ₁/Δ = 8/10 = 4/5
```

---

### 4.4 — Invertibility and Homogeneous Systems

**Theorem — Equivalent Conditions:**

For a square matrix A, the following statements are all equivalent:

```
(1) A is invertible (A⁻¹ exists)
(2) AX = 0 has only the trivial solution X = 0
(3) det(A) ≠ 0
```

**Homogeneous System AX = 0:**

```
AX = 0 has a NONTRIVIAL solution  ⟺  |A| = 0
AX = 0 has ONLY the trivial solution  ⟺  |A| ≠ 0
```

> A homogeneous system always has at least one solution (X = 0). The question is whether there are others.

---

### 4.5 — Finding A⁻¹ Using Cramer's Rule / Adjoint

Since `A⁻¹ = (1/|A|) · adj(A)`, we can find the inverse without row reduction.

**Example:** Find A⁻¹ for:
```
    |  1   2  −4 |
A = |  0   2   3 |
    |  1   1  −1 |
```

Step 1 — Compute |A|:
```
|A| = 1[(2)(−1)−(3)(1)] − 2[(0)(−1)−(3)(1)] + (−4)[(0)(1)−(2)(1)]
    = 1[−2−3] − 2[0−3] − 4[0−2]
    = −5 + 6 + 8
    = 9
```

Step 2 — Compute all cofactors and form adj(A):
```
C₁₁ = +(2·(−1)−3·1)   = −5      C₁₂ = −(0·(−1)−3·1) =  3     C₁₃ = +(0·1−2·1)  = −2
C₂₁ = −(2·(−1)−(−4)·1)= 6       C₂₂ = +(1·(−1)−(−4)·1)= 3    C₂₃ = −(1·1−2·1)  =  1
C₃₁ = +(2·3−(−4)·2)   = 14      C₃₂ = −(1·3−(−4)·0) = −3     C₃₃ = +(1·2−2·0)  =  2

           | −5   6   14 |
adj(A)  =  |  3   3   −3 |   (transpose of cofactor matrix)
           | −2   1    2 |
```

Step 3 — Compute A⁻¹:
```
           1          | −5   6   14 |
A⁻¹  =   ─── · adj = |  3   3   −3 |
           9          | −2   1    2 |
```

---

## Quick Reference Summary

```
═══════════════════════════════════════════════════════
  DETERMINANT QUICK REFERENCE
═══════════════════════════════════════════════════════

  2×2 det:
      |a  b|
      |c  d|  =  ad − bc

  ───────────────────────────────────────────────────
  Cofactor:
      Cᵢⱼ = (−1)^(i+j) × Mᵢⱼ

  Sign pattern:   + − + − …
                  − + − + …
                  + − + − …

  ───────────────────────────────────────────────────
  Laplace Expansion:
      det(A) = Σⱼ aᵢⱼ · Cᵢⱼ   (expand along any row i)
      det(A) = Σᵢ aᵢⱼ · Cᵢⱼ   (expand along any column j)

  ───────────────────────────────────────────────────
  Row Operations:
      Swap rows          →  det changes sign
      Multiply row by λ  →  det multiplied by λ
      Add multiple       →  det unchanged

  ───────────────────────────────────────────────────
  Key Properties:
      |Aᵀ|    =  |A|
      |AB|    =  |A| · |B|
      |A⁻¹|   =  1/|A|
      |λA|    =  λⁿ · |A|     (n×n matrix)
      |adj(A)| = |A|ⁿ⁻¹

  Triangular matrix → |A| = product of diagonal entries

  ───────────────────────────────────────────────────
  Adjoint:
      adj(A) = transpose of cofactor matrix
      A⁻¹    = (1/|A|) · adj(A)

  ───────────────────────────────────────────────────
  Cramer's Rule:
      Δ  = |A|
      Δᵢ = |A with column i replaced by b|
      xᵢ = Δᵢ / Δ       [requires Δ ≠ 0]

  ───────────────────────────────────────────────────
  Invertibility:
      A invertible  ⟺  AX=0 has only trivial solution  ⟺  |A| ≠ 0
      AX=0 has nontrivial solution  ⟺  |A| = 0

═══════════════════════════════════════════════════════
```