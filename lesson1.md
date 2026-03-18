# 📐 Row Echelon Form (REF) & Reduced Row Echelon Form (RREF)

> A review guide for matrix reduction techniques in Linear Algebra

---

## 📌 Table of Contents

1. [What is a Matrix?](#what-is-a-matrix)
2. [Elementary Row Operations](#elementary-row-operations)
3. [Row Echelon Form (REF)](#row-echelon-form-ref)
4. [Reduced Row Echelon Form (RREF)](#reduced-row-echelon-form-rref)
5. [REF vs RREF — Side-by-Side](#ref-vs-rref--side-by-side)
6. [Step-by-Step: Gaussian Elimination (→ REF)](#step-by-step-gaussian-elimination--ref)
7. [Step-by-Step: Gauss-Jordan Elimination (→ RREF)](#step-by-step-gauss-jordan-elimination--rref)
8. [Applications](#applications)
9. [Common Mistakes](#common-mistakes)
10. [Quick Checklist](#quick-checklist)

---

## What is a Matrix?

A **matrix** is a rectangular array of numbers arranged in rows and columns.

```
      col1  col2  col3
row1 [  2    1   -1  ]
row2 [ -3   -1    2  ]
row3 [ -2    1    2  ]
```

An **augmented matrix** `[A | b]` represents a system of linear equations:

```
2x  +  y  -  z  =  8       →    [ 2   1  -1 |  8 ]
-3x -  y  + 2z  = -11      →    [-3  -1   2 | -11]
-2x +  y  + 2z  = -3       →    [-2   1   2 | -3 ]
```

---

## Elementary Row Operations

These are the **only** operations allowed during row reduction. They do **not** change the solution set.

| Operation | Notation | Description |
|-----------|----------|-------------|
| Row Swap | `Rᵢ ↔ Rⱼ` | Swap row i and row j |
| Row Scale | `Rᵢ → k·Rᵢ` | Multiply row i by nonzero scalar k |
| Row Addition | `Rᵢ → Rᵢ + k·Rⱼ` | Add k times row j to row i |

---

## Row Echelon Form (REF)

### ✅ Definition

A matrix is in **Row Echelon Form** if it satisfies all three conditions:

1. **All zero rows are at the bottom.**
2. **Each leading entry (pivot) is to the right of the pivot in the row above.**
3. **All entries below a pivot are zero.**

> 💡 The leading entry in each nonzero row is called a **pivot**. It does not have to be 1.

### Example — REF ✅

```
[ 2   1  -1 |  8 ]
[ 0   2  -5 |  4 ]
[ 0   0   3 |  9 ]
```

- Pivots: `2`, `2`, `3` (staircase pattern ✓)
- Zeros below each pivot ✓
- No zero rows on top ✓

### Example — NOT REF ❌

```
[ 0   1  -5 |  4 ]    ← zero in first column but row is above a nonzero row
[ 2   1  -1 |  8 ]
[ 0   0   3 |  9 ]
```

---

## Reduced Row Echelon Form (RREF)

### ✅ Definition

A matrix is in **Reduced Row Echelon Form** if it satisfies all REF conditions **plus**:

4. **Each pivot is exactly `1`** (called a **leading 1**)
5. **All entries above AND below each pivot are zero.**

> 💡 RREF is **unique** — every matrix has exactly one RREF. REF is not unique.

### Example — RREF ✅

```
[ 1   0   0 |  2 ]
[ 0   1   0 | -1 ]
[ 0   0   1 |  3 ]
```

- All pivots are 1 ✓
- Zeros above and below each pivot ✓
- Solution reads directly: `x = 2`, `y = -1`, `z = 3`

### Example — NOT RREF ❌

```
[ 1   2   0 |  5 ]    ← entry above pivot in column 2 is nonzero (2 ≠ 0)
[ 0   1   0 |  3 ]
[ 0   0   1 |  2 ]
```

---

## REF vs RREF — Side-by-Side

| Feature | REF | RREF |
|---------|-----|------|
| Zero rows at bottom | ✅ | ✅ |
| Staircase pattern (pivots go right) | ✅ | ✅ |
| Zeros **below** each pivot | ✅ | ✅ |
| Each pivot equals **1** | ❌ (optional) | ✅ |
| Zeros **above** each pivot | ❌ | ✅ |
| Uniqueness | ❌ (not unique) | ✅ (always unique) |
| Used in | Gaussian Elimination | Gauss-Jordan Elimination |
| Easiest next step | Back-substitution | Read solution directly |

---

## Step-by-Step: Gaussian Elimination (→ REF)

**Goal:** Reduce to REF, then use **back-substitution**.

### Example: Solve the system

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

**Step 1:** Eliminate below pivot in column 1.

```
R2 → R2 + (3/2)·R1
R3 → R3 + (1)·R1
```

```
[ 2    1    -1  |   8  ]
[ 0   1/2  1/2  |   1  ]
[ 0    2    1   |   5  ]
```

**Step 2:** Eliminate below pivot in column 2.

```
R3 → R3 + (-4)·R2
```

```
[ 2    1    -1  |   8  ]
[ 0   1/2  1/2  |   1  ]
[ 0    0   -1   |   1  ]
```

✅ **Now in REF!**

**Back-substitution:**
- From R3: `-z = 1` → `z = -1`
- From R2: `y/2 + (-1)/2 = 1` → `y = 3`
- From R1: `2x + 3 - (-1) = 8` → `x = 2`

**Solution: `x = 2, y = 3, z = -1`**

---

## Step-by-Step: Gauss-Jordan Elimination (→ RREF)

**Goal:** Continue from REF until RREF — read solution directly.

**Continuing from REF above:**

```
[ 2    1    -1  |   8  ]
[ 0   1/2  1/2  |   1  ]
[ 0    0   -1   |   1  ]
```

**Step 3:** Make each pivot = 1.

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

**Step 4:** Eliminate **above** each pivot (back to zero).

```
R2 → R2 - (1)·R3
R1 → R1 + (1/2)·R3
```

```
[ 1   1/2   0   |  7/2 ]
[ 0    1    0   |   3  ]
[ 0    0    1   |  -1  ]
```

```
R1 → R1 - (1/2)·R2
```

```
[ 1    0    0   |   2  ]
[ 0    1    0   |   3  ]
[ 0    0    1   |  -1  ]
```

✅ **Now in RREF!**

**Solution reads directly: `x = 2, y = 3, z = -1`** ✓

---

## Applications

### 1. Solving Systems of Linear Equations
- REF + back-substitution → Gaussian Elimination
- RREF → Gauss-Jordan Elimination (direct read)

### 2. Finding the Rank of a Matrix
> **Rank** = number of nonzero rows in REF/RREF = number of pivots

```
rank(A) = number of pivot columns
```

### 3. Determining Solution Types

| Condition | Solution Type |
|-----------|--------------|
| Each variable is a pivot variable | **Unique solution** |
| Free variables exist (non-pivot columns) | **Infinitely many solutions** |
| Row `[ 0 0 0 | c ]` where `c ≠ 0` | **No solution (inconsistent)** |

### 4. Finding the Inverse of a Matrix
Augment with identity: `[A | I]` → row reduce → `[I | A⁻¹]`

### 5. Computing the Null Space / Column Space
Identify free variables from RREF pivot structure.

---

## Common Mistakes

### ❌ Mistake 1: Forgetting the pivot must move **right**
Every new pivot must be strictly to the right of the one above it.

### ❌ Mistake 2: Confusing REF and RREF
In REF, you only zero out **below** pivots. In RREF, you must also zero out **above**.

### ❌ Mistake 3: Pivots don't need to be 1 in REF
REF allows any nonzero pivot. Only RREF requires leading 1s.

### ❌ Mistake 4: Arithmetic errors during row operations
Always double-check: `Rᵢ → Rᵢ + k·Rⱼ` — the row being **added to** changes, not the one being multiplied.

### ❌ Mistake 5: Stopping back-substitution too early
In Gaussian Elimination, make sure to substitute all the way back up.

---

## Quick Checklist

### Is this matrix in REF? ✅ Check:
- [ ] All-zero rows are at the bottom
- [ ] Each pivot is strictly to the right of the pivot in the row above
- [ ] Everything **below** each pivot is zero

### Is this matrix in RREF? ✅ Check all REF conditions, plus:
- [ ] Every pivot equals exactly **1**
- [ ] Everything **above** each pivot is also zero

---

## Summary

```
Original Matrix
      ↓  (Elementary Row Operations)
   REF  ←  Gaussian Elimination
      ↓  (More row operations)
  RREF  ←  Gauss-Jordan Elimination
```

| | REF | RREF |
|--|-----|------|
| **Algorithm** | Gaussian Elimination | Gauss-Jordan Elimination |
| **To solve** | + Back-substitution needed | Read directly |
| **Uniqueness** | Not unique | Always unique |
| **Pivot requirement** | Any nonzero value | Must be 1 |
| **Zeros required** | Below pivot only | Above AND below pivot |

---

*Linear Algebra Review — REF & RREF*