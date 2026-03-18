# 📐 Linear Algebra — Complete Review Guide

> Covers: Systems of Linear Equations · Elementary Row Operations · Elementary Matrices ·
> Gaussian & Gauss-Jordan Elimination · REF · RREF · Homogeneous Systems · Solution Vectors ·
> Inverse of a Matrix · Invertibility · Rank · Nullity

---

## 📌 Table of Contents

1. [What is a Matrix?](#1-what-is-a-matrix)
2. [System of Linear Equations](#2-system-of-linear-equations)
3. [Elementary Row Operations](#3-elementary-row-operations)
4. [Elementary Matrices](#4-elementary-matrices)
5. [Row Echelon Form (REF)](#5-row-echelon-form-ref)
6. [Reduced Row Echelon Form (RREF)](#6-reduced-row-echelon-form-rref)
7. [REF vs RREF — Side-by-Side](#7-ref-vs-rref--side-by-side)
8. [Gaussian Elimination (→ REF)](#8-gaussian-elimination--ref)
9. [Gauss-Jordan Elimination (→ RREF)](#9-gauss-jordan-elimination--rref)
10. [Homogeneous Systems](#10-homogeneous-systems)
11. [Solution Vector](#11-solution-vector)
12. [Inverse of a Matrix (A⁻¹)](#12-inverse-of-a-matrix-a)
13. [Invertible (Non-singular) Matrix](#13-invertible-non-singular-matrix)
14. [Rank of a Matrix](#14-rank-of-a-matrix)
15. [Nullity of a Matrix](#15-nullity-of-a-matrix)
16. [Rank–Nullity Theorem](#16-ranknullity-theorem)
17. [Common Mistakes](#17-common-mistakes)
18. [Quick Checklists](#18-quick-checklists)
19. [Summary Map](#19-summary-map)

---

## 1. What is a Matrix?

A **matrix** is a rectangular array of numbers with `m` rows and `n` columns — called an **m × n matrix**.

```
A =  [ a₁₁  a₁₂  a₁₃ ]    (2 × 3 matrix)
     [ a₂₁  a₂₂  a₂₃ ]
```

### Special matrices

| Name | Description | Example |
|------|-------------|---------|
| **Square matrix** | m = n | 3×3 |
| **Identity matrix Iₙ** | 1s on diagonal, 0s elsewhere | diag(1,1,1) |
| **Zero matrix** | All entries are 0 | |
| **Augmented matrix [A\|b]** | A with solution column appended | Used for solving Ax = b |

---

## 2. System of Linear Equations

### Definition

A **system of m linear equations in n unknowns** has the form:

```
a₁₁x₁ + a₁₂x₂ + ··· + a₁ₙxₙ = b₁
a₂₁x₁ + a₂₂x₂ + ··· + a₂ₙxₙ = b₂
              ⋮
aₘ₁x₁ + aₘ₂x₂ + ··· + aₘₙxₙ = bₘ
```

Written compactly as **Ax = b**, where:

```
A = coefficient matrix (m×n)
x = unknown vector (n×1)
b = constant vector (m×1)
```

### Three Possible Solution Types

| Case | Description | Geometric (2D) |
|------|-------------|----------------|
| **Unique solution** | Exactly one solution | Lines intersect at one point |
| **Infinitely many solutions** | A free variable exists | Lines are the same (coincident) |
| **No solution** | Inconsistent system | Lines are parallel |

### How to check via augmented matrix

```
[ A | b ] → row reduce → check result

Row [ 0  0  ···  0 | c ] with c ≠ 0  →  NO SOLUTION
All variables are pivot variables        →  UNIQUE SOLUTION
At least one free variable exists        →  INFINITELY MANY SOLUTIONS
```

### Example

```
x  + 2y  =  5
3x -  y  =  1
```

Augmented matrix:
```
[ 1   2 | 5 ]
[ 3  -1 | 1 ]
```

R2 → R2 - 3·R1:
```
[ 1   2 |  5 ]
[ 0  -7 | -14]
```

R2 → (-1/7)·R2:
```
[ 1   2 | 5 ]
[ 0   1 | 2 ]
```

R1 → R1 - 2·R2:
```
[ 1   0 | 1 ]
[ 0   1 | 2 ]
```

**Solution: x = 1, y = 2** (unique solution ✓)

---

## 3. Elementary Row Operations

These are the **only three** operations allowed when manipulating a matrix. They preserve the solution set of the system.

| # | Operation | Symbol | Description |
|---|-----------|--------|-------------|
| **Type 1** | Row Swap | `Rᵢ ↔ Rⱼ` | Swap row i and row j |
| **Type 2** | Row Scaling | `Rᵢ → k·Rᵢ` | Multiply row i by nonzero scalar k (k ≠ 0) |
| **Type 3** | Row Addition | `Rᵢ → Rᵢ + k·Rⱼ` | Add k times row j to row i |

> 💡 These operations are **reversible** — you can always undo them, which means the original and reduced matrices have the **same solution set**. Two matrices related by row operations are called **row equivalent**.

### Key Properties
- Row swap does **not** change the solution set, but **changes the sign** of the determinant.
- Row scaling by `k` **multiplies** the determinant by `k`.
- Row addition does **not** change the determinant.

### Example of Each Operation

Starting with:
```
A = [ 1   2   3 ]
    [ 4   5   6 ]
    [ 7   8   9 ]
```

**Type 1 — R1 ↔ R2:**
```
[ 4   5   6 ]
[ 1   2   3 ]
[ 7   8   9 ]
```

**Type 2 — R1 → 2·R1:**
```
[ 2   4   6 ]
[ 4   5   6 ]
[ 7   8   9 ]
```

**Type 3 — R2 → R2 - 4·R1:**
```
[ 1   2   3 ]
[ 0  -3  -6 ]
[ 7   8   9 ]
```

---

## 4. Elementary Matrices

### Definition

An **elementary matrix** is a matrix obtained by performing **exactly one** elementary row operation on the **identity matrix Iₙ**.

> 💡 Every elementary row operation on matrix A can be achieved by **left-multiplying** A by a corresponding elementary matrix E:
> `E · A = (result of the row operation on A)`

### The Three Types of Elementary Matrices

#### Type 1 — Row Swap (Swap Rᵢ ↔ Rⱼ)

Swap rows 1 and 2 in I₃:
```
I₃ = [ 1  0  0 ]    →  E = [ 0  1  0 ]
     [ 0  1  0 ]            [ 1  0  0 ]
     [ 0  0  1 ]            [ 0  0  1 ]
```

Effect: `E · A` swaps rows 1 and 2 of A.

#### Type 2 — Row Scaling (Rᵢ → k·Rᵢ)

Multiply row 2 by k = 3 in I₃:
```
I₃ = [ 1  0  0 ]    →  E = [ 1  0  0 ]
     [ 0  1  0 ]            [ 0  3  0 ]
     [ 0  0  1 ]            [ 0  0  1 ]
```

Effect: `E · A` multiplies row 2 of A by 3.

#### Type 3 — Row Addition (Rᵢ → Rᵢ + k·Rⱼ)

Add 2·R1 to R3 (k=2, i=3, j=1) in I₃:
```
I₃ = [ 1  0  0 ]    →  E = [ 1  0  0 ]
     [ 0  1  0 ]            [ 0  1  0 ]
     [ 0  0  1 ]            [ 2  0  1 ]
```

Effect: `E · A` adds 2 times row 1 of A to row 3.

### Inverses of Elementary Matrices

Every elementary matrix is **invertible**, and its inverse is also an elementary matrix:

| Type | Operation E | Inverse E⁻¹ |
|------|-------------|-------------|
| Row Swap | `Rᵢ ↔ Rⱼ` | `Rᵢ ↔ Rⱼ` (swap again) |
| Row Scale by k | `Rᵢ → k·Rᵢ` | `Rᵢ → (1/k)·Rᵢ` |
| Row Addition by k | `Rᵢ → Rᵢ + k·Rⱼ` | `Rᵢ → Rᵢ - k·Rⱼ` |

### Connecting to Inverse of A

If A is reduced to RREF by operations E₁, E₂, ..., Eₖ (applied left to right), then:

```
Eₖ · ··· · E₂ · E₁ · A = I

Therefore:  A⁻¹ = Eₖ · ··· · E₂ · E₁
```

This is the theoretical basis for the `[A | I] → [I | A⁻¹]` method.

---

## 5. Row Echelon Form (REF)

### Definition

A matrix is in **Row Echelon Form (REF)** if:

1. **All-zero rows are at the bottom.**
2. **Each leading entry (pivot) is strictly to the right of the pivot in the row above.**
3. **All entries below a pivot are zero.**

> 💡 The pivot does **not** need to equal 1 in REF.

### Example — REF ✅

```
[ 2   1  -1 |  8 ]   ← pivot at column 1
[ 0   3  -5 |  4 ]   ← pivot at column 2 (right of above ✓)
[ 0   0   7 |  9 ]   ← pivot at column 3 (right of above ✓)
```

### Example — NOT REF ❌

```
[ 0   1  -5 |  4 ]   ← starts with 0, but row below has pivot at col 1 ❌
[ 2   1  -1 |  8 ]
[ 0   0   3 |  9 ]
```

---

## 6. Reduced Row Echelon Form (RREF)

### Definition

A matrix is in **Reduced Row Echelon Form (RREF)** if it satisfies all REF conditions **plus**:

4. **Each pivot (leading entry) is exactly `1`** — called a **leading 1**.
5. **All entries above AND below each pivot are zero.**

> 💡 Every matrix has **exactly one** RREF — it is unique. REF is NOT unique.

### Example — RREF ✅

```
[ 1   0   0 |  2 ]   ← pivot = 1, entire column 1 is zero except here
[ 0   1   0 | -1 ]   ← pivot = 1, entire column 2 is zero except here
[ 0   0   1 |  3 ]   ← pivot = 1, entire column 3 is zero except here
```

Solution reads directly: **x₁ = 2, x₂ = -1, x₃ = 3**

### Example — NOT RREF ❌

```
[ 1   2   0 |  5 ]   ← the 2 above pivot in column 2 is nonzero ❌
[ 0   1   0 |  3 ]
[ 0   0   1 |  2 ]
```

---

## 7. REF vs RREF — Side-by-Side

| Feature | REF | RREF |
|---------|-----|------|
| Zero rows at bottom | ✅ | ✅ |
| Staircase (pivots go right) | ✅ | ✅ |
| Zeros **below** each pivot | ✅ | ✅ |
| Each pivot = **1** | ❌ optional | ✅ required |
| Zeros **above** each pivot | ❌ | ✅ |
| Uniqueness | ❌ not unique | ✅ always unique |
| Method | Gaussian Elimination | Gauss-Jordan Elimination |
| Next step needed | Back-substitution | Read solution directly |

---

## 8. Gaussian Elimination (→ REF)

### Algorithm

```
1. Find the leftmost nonzero column (pivot column)
2. Swap rows if needed: put a nonzero entry at the top of the pivot column
3. Use row addition to make all entries BELOW the pivot = 0
4. Ignore the current row, move to the submatrix below
5. Repeat until in REF
6. Solve using back-substitution (bottom → top)
```

### Full Worked Example

**System:**
```
2x  +  y  -  z  =  8
-3x -  y  + 2z  = -11
-2x +  y  + 2z  = -3
```

**Augmented matrix:**
```
[ 2   1  -1 |  8  ]
[-3  -1   2 | -11 ]
[-2   1   2 | -3  ]
```

**Step 1:** Pivot column 1. Eliminate below pivot (2).
```
R2 → R2 + (3/2)·R1       [multiplier = 3/2]
R3 → R3 + (1)·R1         [multiplier = 1]
```
```
[ 2    1    -1  |   8  ]
[ 0   1/2  1/2  |   1  ]
[ 0    2    1   |   5  ]
```

**Step 2:** Pivot column 2. Eliminate below pivot (1/2).
```
R3 → R3 + (-4)·R2        [multiplier = -4]
```
```
[ 2    1    -1  |   8  ]
[ 0   1/2  1/2  |   1  ]
[ 0    0   -1   |   1  ]
```

✅ **Now in REF!**

**Back-substitution:**
```
Row 3:  -z = 1                          →  z = -1
Row 2:  (1/2)y + (1/2)(-1) = 1         →  y = 3
Row 1:  2x + (3) - (-1) = 8            →  x = 2
```

**✅ Solution: x = 2, y = 3, z = -1**

---

## 9. Gauss-Jordan Elimination (→ RREF)

### Algorithm

```
1. Perform Gaussian Elimination to reach REF
2. Scale each pivot row so the pivot = 1
3. Eliminate ALL entries ABOVE each pivot (work bottom-up)
4. → Result: RREF
5. Read solution directly from the augmented column
```

### Full Worked Example (continuing from REF above)

**From REF:**
```
[ 2    1    -1  |   8  ]
[ 0   1/2  1/2  |   1  ]
[ 0    0   -1   |   1  ]
```

**Step 3:** Scale rows so all pivots = 1.
```
R1 → (1/2)·R1
R2 → 2·R2
R3 → (-1)·R3
```
```
[ 1   1/2  -1/2  |   4  ]
[ 0    1    1    |   2  ]
[ 0    0    1    |  -1  ]
```

**Step 4a:** Eliminate above pivot in column 3.
```
R2 → R2 - (1)·R3
R1 → R1 + (1/2)·R3
```
```
[ 1   1/2   0   |  7/2 ]
[ 0    1    0   |   3  ]
[ 0    0    1   |  -1  ]
```

**Step 4b:** Eliminate above pivot in column 2.
```
R1 → R1 - (1/2)·R2
```
```
[ 1    0    0   |   2  ]
[ 0    1    0   |   3  ]
[ 0    0    1   |  -1  ]
```

✅ **Now in RREF!**

**✅ Solution reads directly: x = 2, y = 3, z = -1**

### Gaussian vs Gauss-Jordan Comparison

| | Gaussian Elimination | Gauss-Jordan Elimination |
|-|----------------------|--------------------------|
| **Target form** | REF | RREF |
| **Pivot requirement** | Any nonzero | Must be 1 |
| **Zero above pivot?** | No | Yes |
| **How to solve** | Back-substitution | Read directly |
| **Work required** | Less upfront | More upfront, no back-sub |

---

## 10. Homogeneous Systems

### Definition

A system **Ax = 0** (where the right-hand side is the zero vector) is called a **homogeneous system**.

```
a₁₁x₁ + a₁₂x₂ + ··· + a₁ₙxₙ = 0
a₂₁x₁ + a₂₂x₂ + ··· + a₂ₙxₙ = 0
                 ⋮
aₘ₁x₁ + aₘ₂x₂ + ··· + aₘₙxₙ = 0
```

### Key Facts

> 💡 A homogeneous system **ALWAYS** has at least one solution: the **trivial solution** x = 0 (all variables = 0).

| Condition | Solution |
|-----------|----------|
| No free variables (all pivot) | **Trivial solution only** (x = 0) |
| At least one free variable | **Infinitely many solutions** (nontrivial solutions exist) |
| Homogeneous system | **NEVER inconsistent** |

### Example

**System:**
```
x  + 2y  - z  = 0
2x -  y  + z  = 0
```

Augmented matrix (right side stays 0 throughout):
```
[ 1   2  -1 | 0 ]
[ 2  -1   1 | 0 ]
```

R2 → R2 - 2·R1:
```
[ 1   2  -1 | 0 ]
[ 0  -5   3 | 0 ]
```

R2 → (-1/5)·R2:
```
[ 1   2   -1  | 0 ]
[ 0   1  -3/5 | 0 ]
```

R1 → R1 - 2·R2:
```
[ 1   0   1/5 | 0 ]
[ 0   1  -3/5 | 0 ]
```

- **Pivot variables:** x₁, x₂
- **Free variable:** x₃ = t

**General solution:**
```
x₁ = -t/5
x₂ =  3t/5
x₃ =  t
```

→ Written as a **solution vector** (see Section 11).

---

## 11. Solution Vector

### Definition

A **solution vector** expresses the complete solution of a system in column-vector form. For systems with free variables, the solution is split into a **particular** part and a **homogeneous** part.

### Forms

#### Case 1: Unique Solution

Write as a single column vector:
```
x = [ x₁ ]   [  2 ]
    [ x₂ ] = [  3 ]
    [ x₃ ]   [ -1 ]
```

#### Case 2: Homogeneous System with Free Variables (Ax = 0)

```
x = t · v₁  +  s · v₂  +  ···      (one vector per free variable)
```

From Section 10 (one free variable t):
```
x = t · [ -1/5 ]
        [  3/5 ]
        [   1  ]
```

This vector spans the **null space** of A.

#### Case 3: Non-homogeneous with Free Variables (Ax = b)

```
x = x_p  +  t · x_h

where:
  x_p = particular solution (any one specific solution of Ax = b)
  x_h = null space vector(s) (solution(s) of Ax = 0)
```

**Example:** Suppose RREF gives:
```
x₁ =  3 - 2t
x₂ =  t         (free)
x₃ = -1
```

Then:
```
x = [  3 ]       [ -2 ]
    [  0 ]  + t · [  1 ]
    [ -1 ]       [  0 ]

    x_particular   x_homogeneous
```

### Solution Set Notation

For multiple free variables (t, s, ...):
```
Solution set = { x_p + t·v₁ + s·v₂ + ··· | t, s ∈ ℝ }
```

---

## 12. Inverse of a Matrix (A⁻¹)

### Definition

For a **square n×n matrix A**, the **inverse A⁻¹** is the unique matrix satisfying:

```
A · A⁻¹ = A⁻¹ · A = Iₙ
```

> 💡 Only **square** matrices can potentially have an inverse. Not all square matrices are invertible.

### Finding A⁻¹ by Gauss-Jordan Elimination

**Method:** Form `[A | I]` and row reduce until the left side becomes I.

```
[A | I]  →  row reduce  →  [I | A⁻¹]
```

If the left side **cannot** become I, then A is **not invertible**.

### Full Worked Example

Find A⁻¹ for:
```
A = [ 1   2 ]
    [ 3   4 ]
```

**Set up [A | I₂]:**
```
[ 1   2 | 1   0 ]
[ 3   4 | 0   1 ]
```

**R2 → R2 - 3·R1:**
```
[ 1   2 |  1   0 ]
[ 0  -2 | -3   1 ]
```

**R2 → (-1/2)·R2:**
```
[ 1   2 |  1    0  ]
[ 0   1 | 3/2  -1/2]
```

**R1 → R1 - 2·R2:**
```
[ 1   0 | -2    1  ]
[ 0   1 | 3/2  -1/2]
```

✅ Left side = I, so:
```
A⁻¹ = [  -2     1  ]
      [  3/2  -1/2  ]
```

**Verify A · A⁻¹ = I:**
```
[ 1  2 ]·[  -2    1  ] = [ 1·(-2)+2·(3/2)   1·1+2·(-1/2) ] = [ 1  0 ] ✓
[ 3  4 ]  [ 3/2  -1/2]   [ 3·(-2)+4·(3/2)   3·1+4·(-1/2) ]   [ 0  1 ]
```

### 2×2 Inverse Shortcut

For `A = [[a, b], [c, d]]`:
```
A⁻¹ = 1/(ad - bc) · [  d  -b ]
                     [ -c   a ]

Valid only when  det(A) = ad - bc ≠ 0
```

### Key Properties of Inverses

```
(A⁻¹)⁻¹ = A
(AB)⁻¹  = B⁻¹ · A⁻¹     ← order REVERSES
(Aᵀ)⁻¹  = (A⁻¹)ᵀ
(kA)⁻¹  = (1/k) · A⁻¹   for scalar k ≠ 0
```

---

## 13. Invertible (Non-singular) Matrix

### Definition

A square matrix A is **invertible** (also called **non-singular**) if there exists A⁻¹ such that `A·A⁻¹ = A⁻¹·A = I`.

A matrix that is **not invertible** is called **singular**.

### The Invertible Matrix Theorem

For an **n×n square matrix A**, the following statements are **all equivalent**:

| # | Statement |
|---|-----------|
| 1 | A is invertible |
| 2 | A⁻¹ exists |
| 3 | det(A) ≠ 0 |
| 4 | RREF of A is Iₙ |
| 5 | A has n pivot positions |
| 6 | rank(A) = n |
| 7 | Ax = 0 has only the trivial solution x = 0 |
| 8 | Ax = b has exactly one solution for every b |
| 9 | The columns of A are linearly independent |
| 10 | nullity(A) = 0 |

> 💡 If **any one** of these is true, then **all** of them are true. If any one is false, all are false.

### Quick Test: Is A Invertible?

```
Method 1 — Determinant:
  det(A) = 0   →  SINGULAR (not invertible) ❌
  det(A) ≠ 0  →  INVERTIBLE ✅

Method 2 — Row reduction:
  RREF(A) = I  →  INVERTIBLE ✅
  RREF(A) ≠ I  →  SINGULAR ❌
```

### Example

```
A = [ 1   2 ]    det(A) = 1·4 - 2·3 = -2  ≠ 0  →  INVERTIBLE ✅
    [ 3   4 ]

B = [ 1   2 ]    det(B) = 1·4 - 2·2 = 0   = 0  →  SINGULAR ❌
    [ 2   4 ]
```

---

## 14. Rank of a Matrix

### Definition

The **rank** of matrix A is the number of **pivot positions** in its RREF (= number of nonzero rows in REF).

```
rank(A) = number of pivot columns
        = number of nonzero rows in REF/RREF
```

### How to Find Rank

```
1. Row reduce A to REF (or RREF)
2. Count the number of nonzero rows
```

### Example

```
A = [ 1   2   3 ]
    [ 4   5   6 ]
    [ 7   8   9 ]
```

Row reduction:
```
R2 → R2 - 4·R1      →    [ 1   2   3 ]
R3 → R3 - 7·R1           [ 0  -3  -6 ]
                          [ 0  -6  -12]

R3 → R3 - 2·R2      →    [ 1   2   3 ]
                          [ 0  -3  -6 ]
                          [ 0   0   0 ]   ← zero row
```

**2 nonzero rows → rank(A) = 2**

### Rank and Consistency of Ax = b

Let A be m×n and [A|b] be the augmented matrix.

| Condition | Solution Type |
|-----------|---------------|
| rank(A) = rank([A\|b]) = n | **Unique solution** |
| rank(A) = rank([A\|b]) < n | **Infinitely many solutions** |
| rank(A) < rank([A\|b]) | **No solution (inconsistent)** |

### Rank Properties

```
rank(A)   ≤  min(m, n)            for any m×n matrix
rank(A)   =  rank(Aᵀ)
rank(AB)  ≤  min(rank(A), rank(B))
For n×n invertible A:  rank(A) = n   (full rank)
```

---

## 15. Nullity of a Matrix

### Null Space (Kernel)

The **null space** of A is the set of all vectors x satisfying Ax = 0:

```
Null(A) = { x ∈ ℝⁿ | Ax = 0 }
```

This is always a **subspace** of ℝⁿ (contains 0, closed under addition and scaling).

### Nullity

The **nullity** of A is the **dimension** of the null space:

```
nullity(A) = dim(Null(A))
           = number of free variables in solution of Ax = 0
           = number of non-pivot columns in RREF(A)
```

### How to Find Null Space and Nullity

```
1. Row reduce A to RREF
2. nullity = count of non-pivot (free) columns
3. Set each free variable = 1, others = 0 (one at a time)
4. Solve for pivot variables → each gives one null space basis vector
```

### Example

```
A = [ 1   2   3 ]     RREF:    [ 1   0  -1 ]
    [ 4   5   6 ]               [ 0   1   2 ]
    [ 7   8   9 ]               [ 0   0   0 ]
```

- **Pivot columns:** 1, 2 → pivot variables: x₁, x₂
- **Non-pivot column:** 3 → free variable: x₃ = t
- **nullity(A) = 1**

**Find null space vector:** Let x₃ = t.

From RREF:
```
x₁ - x₃ = 0   →   x₁ = t
x₂ + 2x₃ = 0  →   x₂ = -2t
x₃ = t
```

**Null space basis vector:**
```
x = t · [  1 ]
        [ -2 ]
        [  1 ]

Null(A) = span{ [1, -2, 1]ᵀ }    (a line through origin in ℝ³)
```

---

## 16. Rank–Nullity Theorem

### Statement

For any **m×n matrix A**:

```
rank(A)  +  nullity(A)  =  n

   (pivots)    (free vars)   (total columns)
```

Also called the **Dimension Theorem** or **Fundamental Theorem of Linear Maps**.

### Intuition

```
n columns = pivot columns + free (non-pivot) columns
          =    rank       +     nullity
```

### Example Verification

From Sections 14 and 15 (A is 3×3):
```
n = 3
rank(A)    = 2    (two pivot columns)
nullity(A) = 1    (one free column)

rank + nullity = 2 + 1 = 3 = n  ✅
```

### Summary Table

| Matrix A (m×n) | rank | nullity | Solution of Ax=b |
|----------------|------|---------|------------------|
| Full column rank (rank = n) | n | 0 | Unique (if consistent) |
| rank < n | < n | > 0 | Infinitely many (if consistent) |
| Full row rank (rank = m) | m | n−m | Always consistent |
| Square, invertible (rank = n = m) | n | 0 | Unique for every b |

---

## 17. Common Mistakes

### ❌ Mistake 1: Forgetting the augmented column in row operations
Apply row operations to the **entire** augmented row, including the `| b` part.

### ❌ Mistake 2: Scaling a row by zero
Row scaling requires k **≠ 0**. `Rᵢ → 0·Rᵢ` is not a valid operation.

### ❌ Mistake 3: Confusing REF and RREF
- REF: zeros only **below** pivots
- RREF: zeros **above AND below** + all pivots = 1

### ❌ Mistake 4: Not using back-substitution after REF
REF alone does not give the answer. Only RREF lets you read directly.

### ❌ Mistake 5: Thinking homogeneous = trivial only
Ax = 0 always has x = 0, but has **more solutions** whenever nullity > 0 (free variables exist).

### ❌ Mistake 6: Confusing rank and nullity
```
rank    = pivot columns      (dimensions the matrix "maps to")
nullity = free columns       (dimensions the matrix "collapses" to zero)
```

### ❌ Mistake 7: Assuming every square matrix is invertible
Square + det = 0 → **singular**, no inverse.

### ❌ Mistake 8: Wrong order for (AB)⁻¹
```
(AB)⁻¹ = B⁻¹A⁻¹    ← NOT A⁻¹B⁻¹
```

### ❌ Mistake 9: Building elementary matrix from A instead of I
Elementary matrix E is always built from **Iₙ**, not from A.

### ❌ Mistake 10: Forgetting rank-nullity when computing one of the two
If you know rank = r and the matrix is m×n, then nullity = n − r immediately.

---

## 18. Quick Checklists

### ✅ Is this matrix in REF?
- [ ] All-zero rows are at the bottom
- [ ] Each pivot is strictly right of the pivot in the row above
- [ ] Everything **below** each pivot is 0

### ✅ Is this matrix in RREF?
- [ ] All REF conditions above ✓
- [ ] Every pivot = exactly **1**
- [ ] Everything **above** each pivot is also 0

### ✅ Is matrix A invertible?
- [ ] A is **square** (n×n)
- [ ] det(A) ≠ 0  **OR**  RREF(A) = I  **OR**  rank(A) = n

### ✅ Solving Ax = b
- [ ] Write augmented matrix [A|b]
- [ ] Row reduce to REF or RREF
- [ ] Check for inconsistency: any row `[ 0 ··· 0 | c ]` with c ≠ 0 → no solution
- [ ] Identify pivot vs free variables
- [ ] Write solution in vector form if free variables exist

### ✅ Finding Rank and Nullity
- [ ] Row reduce A to RREF
- [ ] rank = number of nonzero rows (pivot columns)
- [ ] nullity = number of non-pivot columns
- [ ] Verify: rank + nullity = n (number of columns)

---

## 19. Summary Map

```
SYSTEM OF LINEAR EQUATIONS:  Ax = b
              │
              ▼
     AUGMENTED MATRIX  [A | b]
              │
              ▼  Elementary Row Operations
              │     Type 1: Rᵢ ↔ Rⱼ          (↔ Elementary Matrix: swap rows of I)
              │     Type 2: Rᵢ → k·Rᵢ         (↔ Elementary Matrix: scale row of I)
              │     Type 3: Rᵢ → Rᵢ + k·Rⱼ   (↔ Elementary Matrix: add row of I)
              │
        ┌─────┴──────┐
        ▼            ▼
       REF          RREF
  (Gaussian)   (Gauss-Jordan)
        │            │
        ▼            ▼
   Back-sub     Read directly
        └─────┬──────┘
              ▼
          SOLUTION
      ┌────────────────┐
      ▼                ▼
  Unique (rank=n)   Infinitely many (free vars)
                       │
                       ▼  Write as:
               x = x_p  +  t·x_h₁  +  s·x_h₂  +  ···
                         (null space vectors)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
RANK & NULLITY  (for m×n matrix A):

  rank(A)    = # pivot columns  =  # nonzero rows in REF
  nullity(A) = # free columns   =  dim(Null(A))
  rank(A) + nullity(A) = n              ← Rank–Nullity Theorem

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
INVERTIBILITY  (square n×n matrix A only):

  Invertible  ⟺  det(A) ≠ 0
              ⟺  rank(A) = n
              ⟺  nullity(A) = 0
              ⟺  RREF(A) = I
              ⟺  Ax = 0 has only trivial solution

  To find A⁻¹:
    [A | I]  →  row reduce  →  [I | A⁻¹]
    (theoretically: A⁻¹ = Eₖ·⋯·E₂·E₁ where Eᵢ are elementary matrices)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
HOMOGENEOUS  Ax = 0:

  Always consistent (trivial solution x = 0 always exists)
  Nontrivial solutions exist  ⟺  nullity(A) > 0  ⟺  free variable exists
```

---

*Linear Algebra Complete Review — Systems · Row Operations · Elementary Matrices · REF · RREF · Homogeneous Systems · Solution Vectors · Gaussian & Gauss-Jordan Elimination · Inverse · Invertibility · Rank · Nullity*