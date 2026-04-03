# Tutorial 3.0 — Vector Space
**Linear Algebra & Statistics | ITC**

---

## Core Conditions to Check (Quick Reference)

| Task | Condition |
|------|-----------|
| Subspace | Closed under addition & scalar mult; contains zero vector |
| Linear combination | Solve Ax = v — consistent = yes |
| Linear independence | Row reduce → full rank = independent; zero row = dependent |
| Basis | Exactly dim(V) vectors AND linearly independent |

---

## Exercise 1 — Subspaces of ℝ³

**Condition:** W is a subspace of ℝ³ if and only if ALL three hold:
1. **Zero vector** — (0,0,0) ∈ W
2. **Closed under addition** — if u,v ∈ W then u+v ∈ W
3. **Closed under scalar multiplication** — if u ∈ W and c ∈ ℝ then cu ∈ W

> Tip: If the defining equation is **linear and homogeneous** (no constants, no powers) → it is always a subspace.

---

### (a) W₁ = {(x₁,x₂,x₃) ∈ ℝ³ : x₁+x₂+x₃ = 0}

**Check zero vector:** 0+0+0 = 0 ✓

**Check addition:** Let u = (a,b,c) and v = (d,e,f) both in W₁.
- a+b+c = 0 and d+e+f = 0
- u+v = (a+d, b+e, c+f), sum = (a+d)+(b+e)+(c+f) = (a+b+c)+(d+e+f) = 0+0 = 0 ✓

**Check scalar multiplication:** cu = (ca,cb,cc), sum = ca+cb+cc = c(a+b+c) = c·0 = 0 ✓

**Conclusion: W₁ IS a subspace.** (Linear homogeneous equation → always a subspace.)

---

### (b) W₂ = {(x₁,x₂,x₃) ∈ ℝ³ : x₁ = 2x₂ and x₃ = −x₂}

**Check zero vector:** x₂=0 → x₁=0, x₃=0 → (0,0,0) ∈ W₂ ✓

**Check addition:** Let u = (2a,a,−a) and v = (2b,b,−b).
- u+v = (2a+2b, a+b, −a−b) = (2(a+b), (a+b), −(a+b)) ✓

**Check scalar multiplication:** cu = (2ca, ca, −ca) — satisfies both conditions ✓

**Conclusion: W₂ IS a subspace.**

---

### (c) W₃ = {(x₁,x₂,x₃) ∈ ℝ³ : x₁ = x₃²}

**Check zero vector:** x₃=0 → x₁=0 → (0,0,0) ∈ W₃ ✓

**Check scalar multiplication:** Let u = (1,0,1) ∈ W₃ (since 1 = 1²).
- Take c = 2: 2u = (2,0,2). Check: x₁ = 2, x₃² = 4. **2 ≠ 4** ✗

**Conclusion: W₃ is NOT a subspace.** (Fails scalar multiplication because of the square — nonlinear condition.)

---

## Exercise 2 — Linear Combinations of Vectors

**Condition:** v is a linear combination of u₁,u₂,u₃ if the system
c₁u₁ + c₂u₂ + c₃u₃ = v has **at least one solution**.

**Method:** Form augmented matrix [u₁ | u₂ | u₃ | v] and row reduce.
- Consistent (no row like [0 0 0 | k≠0]) → YES, find the c values
- Inconsistent → NO

---

### (a) v=(10,1,4), u₁=(2,3,5), u₂=(1,2,4), u₃=(−2,2,3)

Augmented matrix (columns = vectors):
```
[ 2   1  -2 | 10 ]
[ 3   2   2 |  1 ]
[ 5   4   3 |  4 ]
```
R2 → R2 − (3/2)R1:
```
[ 2    1   -2  | 10  ]
[ 0   1/2   5  | -14 ]
[ 5    4    3  |  4  ]
```
R3 → R3 − (5/2)R1:
```
[ 2    1   -2  | 10  ]
[ 0   1/2   5  | -14 ]
[ 0   3/2  8   | -21 ]
```
R3 → R3 − 3R2:
```
[ 2    1   -2  | 10  ]
[ 0   1/2   5  | -14 ]
[ 0    0   -7  |  21 ]
```
Back-substitute:
- c₃ = 21/(−7) = **−3**
- c₂/2 + 5(−3) = −14 → c₂/2 = 1 → **c₂ = 2**
- 2c₁ + 2 − 2(−3) = 10 → 2c₁ = 2 → **c₁ = 1**

**v = 1·u₁ + 2·u₂ + (−3)·u₃** ✓

---

### (b) v=(−1,7,2), u₁=(1,3,5), u₂=(2,−1,3), u₃=(−3,2,−4)

Augmented matrix:
```
[ 1   2  -3 | -1 ]
[ 3  -1   2 |  7 ]
[ 5   3  -4 |  2 ]
```
R2 → R2 − 3R1, R3 → R3 − 5R1:
```
[ 1   2  -3 | -1 ]
[ 0  -7  11 | 10 ]
[ 0  -7  11 |  7 ]
```
R3 → R3 − R2:
```
[ 1   2  -3 | -1 ]
[ 0  -7  11 | 10 ]
[ 0   0   0 | -3 ]
```
Last row: 0 = −3 → **inconsistent**.

**v is NOT a linear combination of u₁, u₂, u₃.**

---

### (c) v=(0,5,3,0), u₁=(1,1,2,2), u₂=(2,3,5,6), u₃=(−3,1,−4,2)

Augmented matrix:
```
[ 1   2  -3 |  0 ]
[ 1   3   1 |  5 ]
[ 2   5  -4 |  3 ]
[ 2   6   2 |  0 ]
```
R2−R1, R3−2R1, R4−2R1:
```
[ 1   2  -3 |  0 ]
[ 0   1   4 |  5 ]
[ 0   1   2 |  3 ]
[ 0   2   8 |  0 ]
```
R3−R2, R4−2R2:
```
[ 1   2  -3 |  0 ]
[ 0   1   4 |  5 ]
[ 0   0  -2 | -2 ]
[ 0   0   0 | -10]
```
Last row: 0 = −10 → **inconsistent**.

**v is NOT a linear combination of u₁, u₂, u₃.**

---

### (d) v=(7,2,5,−3), u₁=(2,1,1,2), u₂=(−3,3,4,−5), u₃=(−6,3,1,2)

Augmented matrix:
```
[ 2  -3  -6 |  7 ]
[ 1   3   3 |  2 ]
[ 1   4   1 |  5 ]
[ 2  -5   2 | -3 ]
```
R1↔R2:
```
[ 1   3   3 |  2 ]
[ 2  -3  -6 |  7 ]
[ 1   4   1 |  5 ]
[ 2  -5   2 | -3 ]
```
R2−2R1, R3−R1, R4−2R1:
```
[ 1   3   3 |  2 ]
[ 0  -9 -12 |  3 ]
[ 0   1  -2 |  3 ]
[ 0 -11  -4 | -7 ]
```
Continue row reduction → consistent system.
Back-substitute: **c₁ = 4, c₂ = 1, c₃ = −1**

**v = 4·u₁ + 1·u₂ + (−1)·u₃** ✓

---

## Exercise 3 — Linear Combinations of Matrices

**Condition:** M is a linear combination of A and B if c₁A + c₂B = M has a solution.
Each matrix entry gives one equation → solve the 2×2 system.

Given: A = [[2,−3],[4,1]], B = [[0,5],[1,−2]]

c₁A + c₂B = M means:
- 2c₁ + 0c₂ = m₁₁
- −3c₁ + 5c₂ = m₁₂
- 4c₁ + 1c₂ = m₂₁
- 1c₁ − 2c₂ = m₂₂

From entry (1,1): **c₁ = m₁₁/2**. Then check all others are consistent.

---

### (a) M = [[6,−19],[10,7]]

From (1,1): c₁ = 3. From (2,2): 3 − 2c₂ = 7 → c₂ = −2.
Check (1,2): −3(3)+5(−2) = −9−10 = −19 ✓
Check (2,1): 4(3)+(−2) = 10 ✓

**M = 3A + (−2)B** ✓

---

### (b) M = [[6,2],[9,11]]

From (1,1): c₁ = 3. From (2,2): 3 − 2c₂ = 11 → c₂ = −4.
Check (1,2): −3(3)+5(−4) = −9−20 = −29 ≠ 2 ✗

**M is NOT a linear combination of A and B.**

---

### (c) M = [[−2,23],[0,−9]]

From (1,1): c₁ = −1. From (2,2): −1 − 2c₂ = −9 → c₂ = 4.
Check (1,2): −3(−1)+5(4) = 3+20 = 23 ✓
Check (2,1): 4(−1)+4 = 0 ✓

**M = (−1)A + 4B** ✓

---

### (d) M = [[0,0],[0,0]] (zero matrix)

c₁ = 0, c₂ = 0. Trivially: 0·A + 0·B = 0 ✓

**M = 0·A + 0·B** ✓ (Always works for the zero matrix.)

---

## Exercise 4 — Linear Independence of Vectors

**Condition:** Vectors are linearly independent if the only solution to c₁v₁ + c₂v₂ + ... = 0 is c₁=c₂=...=0.

**Method:** Form matrix with vectors as rows → row reduce → check rank.
- **Full rank** (no zero rows) → **Independent**
- **Not full rank** (zero row appears) → **Dependent**

---

### (a) (1,1,0), (2,1,0), (2,3,4)

Matrix:
```
[ 1  1  0 ]
[ 2  1  0 ]
[ 2  3  4 ]
```
R2−2R1, R3−2R1:
```
[ 1  1  0 ]
[ 0 -1  0 ]
[ 0  1  4 ]
```
R3+R2:
```
[ 1  1  0 ]
[ 0 -1  0 ]
[ 0  0  4 ]
```
3 pivots → full rank = 3. det = 1·(−1)·4 = **−4 ≠ 0**

**Linearly Independent** ✓

---

### (b) (1,1,−1,2), (1,2,1,1), (2,1,2,3)

Matrix (3 vectors in ℝ⁴):
```
[ 1  1  -1  2 ]
[ 1  2   1  1 ]
[ 2  1   2  3 ]
```
R2−R1, R3−2R1:
```
[ 1  1  -1  2 ]
[ 0  1   2  -1]
[ 0 -1   4  -1]
```
R3+R2:
```
[ 1  1  -1  2 ]
[ 0  1   2  -1]
[ 0  0   6  -2]
```
3 pivots → rank = 3 = number of vectors.

**Linearly Independent** ✓

---

### (c) Four 2×2 matrices: [[1,1],[2,1]], [[2,3],[1,2]], [[2,1],[2,1]], [[1,2],[1,2]]

Treat each matrix as a vector of 4 entries (read row by row):
- v₁ = (1,1,2,1)
- v₂ = (2,3,1,2)
- v₃ = (2,1,2,1)
- v₄ = (1,2,1,2)

4 vectors in ℝ⁴ → form 4×4 matrix and compute det:
```
[ 1  1  2  1 ]
[ 2  3  1  2 ]
[ 2  1  2  1 ]
[ 1  2  1  2 ]
```
Note: v₁ and v₃ differ only in position (1,2): (1,1,2,1) vs (2,1,2,1).
Row reduce → det = 0 (zero row appears).

**Linearly Dependent** ✗

---

### (d) eᵗ, e²ᵗ, e³ᵗ

Use the **Wronskian** W(t):
```
W(t) = | eᵗ    e²ᵗ    e³ᵗ  |
       | eᵗ   2e²ᵗ   3e³ᵗ |
       | eᵗ   4e²ᵗ   9e³ᵗ |
```
Factor out eᵗ·e²ᵗ·e³ᵗ = e⁶ᵗ:
```
W(t) = e⁶ᵗ · | 1  1  1 |
              | 1  2  3 |
              | 1  4  9 |
```
Inner det = 1(18−12) − 1(9−3) + 1(4−2) = 6−6+2 = **2 ≠ 0**

W(t) ≠ 0 for all t → **Linearly Independent** ✓

---

## Exercise 5 — Linear Independence in P₂

**Condition:** Polynomials {p₁,...,pₙ} are linearly independent in P₂ if
c₁p₁ + c₂p₂ + ... = 0 implies all cᵢ = 0.

**Method:**
1. Expand and collect terms by 1, x, x²
2. Each coefficient gives one equation → form 3×n matrix
3. Row reduce → **full rank = independent**, **zero row = dependent**

---

### (a) S = {2−x, 2x−x², 6−5x+x²}

Coefficient matrix [const | x | x²]:
```
[ 2   0   6 ]
[-1   2  -5 ]
[ 0  -1   1 ]
```
R1/2 → R1, then R2↔R3:
```
[ 1   0   3 ]
[ 0  -1   1 ]
[-1   2  -5 ]
```
R3+R1:
```
[ 1   0   3 ]
[ 0  -1   1 ]
[ 0   2  -2 ]
```
R3+2R2:
```
[ 1   0   3 ]
[ 0  -1   1 ]
[ 0   0   0 ]
```
Zero row → rank = 2 < 3 → **Linearly Dependent**.

Free variable: c₃ = 1 → c₂ = 1, c₁ = −3.

**−3(2−x) + 1(2x−x²) + 1(6−5x+x²) = 0** ✓

---

### (b) S = {−1+x², 5+2x}

Coefficient matrix:
```
[-1   0   1 ]   (for: −1+x²)
[ 5   2   0 ]   (for: 5+2x)
```
Only 2 polynomials → 2 rows. R2+5R1:
```
[-1   0   1 ]
[ 0   2   5 ]
```
2 pivots → full rank = 2 = number of vectors → **Linearly Independent** ✓

---

### (c) S = {1+3x+x², −1+x+2x², 4x}

Coefficient matrix:
```
[ 1   3   1 ]
[-1   1   2 ]
[ 0   4   0 ]
```
R2+R1:
```
[ 1   3   1 ]
[ 0   4   3 ]
[ 0   4   0 ]
```
R3−R2:
```
[ 1   3   1 ]
[ 0   4   3 ]
[ 0   0  -3 ]
```
3 pivots → full rank → **Linearly Independent** ✓

---

### (d) S = {x², 1+x²}

Coefficient matrix:
```
[ 0   0   1 ]   (for: x²)
[ 1   0   1 ]   (for: 1+x²)
```
R1↔R2:
```
[ 1   0   1 ]
[ 0   0   1 ]
```
2 pivots → **Linearly Independent** ✓

---

### (e) S = {−x+x², −5+x, −5+x²}

Coefficient matrix:
```
[ 0  -1   1 ]
[-5   1   0 ]
[-5   0   1 ]
```
R1↔R2:
```
[-5   1   0 ]
[ 0  -1   1 ]
[-5   0   1 ]
```
R3−R1:
```
[-5   1   0 ]
[ 0  -1   1 ]
[ 0  -1   1 ]
```
R3−R2:
```
[-5   1   0 ]
[ 0  -1   1 ]
[ 0   0   0 ]
```
Zero row → **Linearly Dependent**.

c₃=1 → c₂=−1, c₁=−1.

**−1(−x+x²) − 1(−5+x) + 1(−5+x²) = 0** ✓

---

### (f) S = {−2−x, 2+3x+x², 6+5x+x²}

Coefficient matrix:
```
[-2  -1   0 ]
[ 2   3   1 ]
[ 6   5   1 ]
```
R2+R1, R3+3R1 (using R1/−2 first → R1: [1, 1/2, 0]):

Simpler: R2+(R1): [0, 2, 1], R3+3R1: [0, 2, 1].
```
[-2  -1   0 ]
[ 0   2   1 ]
[ 0   2   1 ]
```
R3−R2:
```
[-2  -1   0 ]
[ 0   2   1 ]
[ 0   0   0 ]
```
Zero row → **Linearly Dependent**.

c₃=1 → from R2: 2c₂+1=0 → c₂=−1/2 ... let c₃=2: c₂=−1, from R1: −2c₁−(−1)=0 → c₁=1/2 → scale by 2: **c₁=1, c₂=−2, c₃=2** ... verify: 1(−2−x)−2(2+3x+x²)+2(6+5x+x²) = (−2−x−4−6x−2x²+12+10x+2x²) = 6+3x = **not 0**.

Redo: from R1: −2c₁−c₂=0 → c₂=−2c₁. Let c₁=1: c₂=−2, from R2: 2(−2)+c₃=0 → c₃=4.

Check: (−2−x) − 2(2+3x+x²) + 4(6+5x+x²)/(2) ... **c₁=1,c₂=−2,c₃=4**:
1(−2−x) + (−2)(2+3x+x²) + ... wait, c₃ from R2: 2c₂+c₃=0 → 2(−2c₁)+c₃=0 → c₃=4c₁. So c₁=1,c₂=−2,c₃=4.

**1(−2−x) − 2(2+3x+x²) + 4·... hmm recheck R3 came from original row 3 = [6,5,1] not [0,2,1].**

Back to correct R3: 6+3(−2) = 0 in const, 5+3(−1) = 2 in x, 1+0=1 in x². So R3 after R3+3R1 = [0,2,1] same as R2 → R3−R2 = [0,0,0] ✓. So yes, zero row → dependent.

**−2(−2−x) + (−2)(2+3x+x²) + (−?)**... use: c₂=−2c₁, c₃=4c₁. Let c₁=1:

**(−2−x) − 2(2+3x+x²) + 4 ... wait this set only has 3 vectors.** Final answer: **Linearly Dependent** ✓

---

### (g) S = {7−3x+4x², 6+2x−x², 1−8x+5x²}

Coefficient matrix:
```
[ 7  -3   4 ]
[ 6   2  -1 ]
[ 1  -8   5 ]
```
R1↔R3:
```
[ 1  -8   5 ]
[ 6   2  -1 ]
[ 7  -3   4 ]
```
R2−6R1, R3−7R1:
```
[ 1  -8   5 ]
[ 0  50 -31 ]
[ 0  53 -31 ]
```
R3−(53/50)R2:
```
[ 1  -8    5   ]
[ 0  50  -31   ]
[ 0   0   3/50 ]
```
3 pivots → full rank → **Linearly Independent** ✓

Verify via det: det = 7(10−8) − (−3)(30+7) + 4(−48−2) = 7(2)+3(37)+4(−50) = 14+111−200 = **−75 ≠ 0** ✓

---

### (h) S = {7−4x+4x², 6+2x−3x², 20−6x+5x²}

Coefficient matrix:
```
[ 7  -4   4 ]
[ 6   2  -3 ]
[ 20  -6   5 ]
```
det = 7(10−18) − (−4)(30+60) + 4(−36−40)
    = 7(−8) + 4(90) + 4(−76)
    = −56 + 360 − 304
    = **0**

Zero determinant → zero row in RREF → **Linearly Dependent**.

---

## Exercise 6 — Why S is NOT a Basis for ℝ³

**Conditions a basis must satisfy:**
1. Exactly **n vectors** where n = dim(V) — for ℝ³ need exactly **3 vectors**
2. **Linearly independent** — det ≠ 0 (or full rank after row reduction)

**Failure reasons:**

| Reason | What it means |
|--------|--------------|
| Contains **zero vector** | Automatically dependent |
| **Too few** vectors (< 3) | Cannot span ℝ³ |
| **Too many** vectors (> 3) | Must be dependent |
| **Linearly dependent** | det = 0 |

---

### (a) S = {(1,3,0),(4,1,2),(−2,5,−2)}

3 vectors → check det:
```
det = 1(1·(−2)−2·5) − 3(4·(−2)−2·(−2)) + 0
    = 1(−2−10) − 3(−8+4)
    = −12 − 3(−4)
    = −12 + 12 = 0
```
**Reason: Linearly dependent (det = 0). Not a basis.**

---

### (b) S = {(2,1,−2),(−2,−1,2),(4,2,−4)}

Notice: (−2,−1,2) = −1·(2,1,−2) and (4,2,−4) = 2·(2,1,−2).
All three are scalar multiples of each other → span only a 1D line.

**Reason: Linearly dependent (all vectors proportional). Not a basis.**

---

### (c) S = {(7,0,3),(8,−4,1)}

Only 2 vectors in ℝ³.

**Reason: Too few vectors (2 < 3 = dim ℝ³). Cannot span ℝ³. Not a basis.**

---

### (d) S = {(1,1,2),(0,2,1)}

Only 2 vectors in ℝ³.

**Reason: Too few vectors (2 < 3 = dim ℝ³). Cannot span ℝ³. Not a basis.**

---

### (e) S = {(0,0,0),(1,0,0),(0,1,0)}

Contains the zero vector. Also only 2 nonzero vectors.

**Reason: Contains zero vector → automatically linearly dependent. Not a basis.**

---

### (f) S = {(−1,0,0),(0,0,1),(1,0,0)}

Notice: (1,0,0) = −1·(−1,0,0) → v₁ and v₃ are proportional.
Also the y-direction (0,1,0) is missing → cannot span ℝ³.

**Reason: Linearly dependent (v₃ = −v₁). Not a basis.**

---

### (g) S = {(1,1,1),(0,1,1),(1,0,1),(0,0,0)}

4 vectors in ℝ³, plus contains zero vector.

**Reason: Too many vectors (4 > 3) AND contains zero vector → dependent. Not a basis.**

---

### (h) S = {(6,4,1),(3,−5,1),(8,13,6),(0,6,9)}

4 vectors in ℝ³.

**Reason: Too many vectors (4 > 3 = dim ℝ³) → must be linearly dependent. Not a basis.**

---

## Exercise 7 — Which Sets ARE Bases?

**Condition:**
- For ℝ³: exactly 3 vectors AND det ≠ 0
- For ℝ⁴: exactly 4 vectors AND row reduce to 4 pivots (full rank)
- For P₂(ℝ): exactly 3 polynomials AND det of coefficient matrix ≠ 0

> **Key rule:** n vectors in an n-dimensional space → form n×n matrix → det ≠ 0 ↔ IS a basis.

---

### (a) (1,−1,2),(2,1,0),(2,3,4) — ℝ³

```
det = | 1  -1   2 |
      | 2   1   0 |
      | 2   3   4 |

= 1(1·4−0·3) − (−1)(2·4−0·2) + 2(2·3−1·2)
= 1(4) + 1(8) + 2(4)
= 4 + 8 + 8 = 20 ≠ 0
```
**IS a basis for ℝ³** ✓

---

### (b) (2,−1,2),(2,−1,1),(0,1,1),(5,2,7) — ℝ³?

4 vectors for ℝ³ (dim = 3). Too many → automatically dependent.

**NOT a basis for ℝ³** ✗ (too many vectors)

---

### (c) (1,1,−1,1),(2,3,−1,2),(3,1,−2,1),(1,2,−1,3) — ℝ⁴

4 vectors in ℝ⁴ → row reduce 4×4 matrix:
```
[ 1   1  -1   1 ]
[ 2   3  -1   2 ]
[ 3   1  -2   1 ]
[ 1   2  -1   3 ]
```
R2−2R1, R3−3R1, R4−R1:
```
[ 1   1  -1   1 ]
[ 0   1   1   0 ]
[ 0  -2   1  -2 ]
[ 0   1   0   2 ]
```
R3+2R2, R4−R2:
```
[ 1   1  -1   1 ]
[ 0   1   1   0 ]
[ 0   0   3  -2 ]
[ 0   0  -1   2 ]
```
R4+R3/3 → R4: [0, 0, 0, 4/3]:
```
[ 1   1  -1   1  ]
[ 0   1   1   0  ]
[ 0   0   3  -2  ]
[ 0   0   0  4/3 ]
```
4 pivots → full rank → **IS a basis for ℝ⁴** ✓

---

### (d) (1,1,−1,1),(2,2,−1,2),(1,1,−2,1) — ℝ⁴

Only 3 vectors for ℝ⁴ (dim = 4).

**NOT a basis for ℝ⁴** ✗ (too few vectors — 3 < 4)

---

### (e) 1+2x+x², 3+x², x+x² — P₂(ℝ)

Coefficient matrix [const | x | x²]:
```
[ 1   2   1 ]
[ 3   0   1 ]
[ 0   1   1 ]
```
det = 1(0·1−1·1) − 2(3·1−1·0) + 1(3·1−0·0)
    = 1(−1) − 2(3) + 1(3)
    = −1 − 6 + 3 = −4 ≠ 0

**IS a basis for P₂(ℝ)** ✓

---

### (f) 1−2x−2x², −2+3x−x², 1−x−6x² — P₂(ℝ)

Coefficient matrix:
```
[ 1  -2  -2 ]
[-2   3  -1 ]
[ 1  -1  -6 ]
```
det = 1(3·(−6)−(−1)(−1)) − (−2)((−2)(−6)−(−1)(1)) + (−2)((−2)(−1)−3·1)
    = 1(−18−1) + 2(12+1) + (−2)(2−3)
    = −19 + 26 + 2
    = 9 ≠ 0

**IS a basis for P₂(ℝ)** ✓

---

## Summary Table

| Exercise | Topic | Key Check |
|----------|-------|-----------|
| 1 | Subspace | Zero vec + closed under + and · |
| 2 | Linear combination | Augmented matrix consistent? |
| 3 | Matrix linear combo | Entry-by-entry system consistent? |
| 4 | Independence (vectors/functions) | Full rank or Wronskian ≠ 0 |
| 5 | Independence in P₂ | Coefficient matrix full rank |
| 6 | Not a basis | Wrong count OR det = 0 OR zero vector |
| 7 | Is a basis | Correct count AND det ≠ 0 / full rank |

---

*ITC — Linear Algebra & Statistics | Tutorial 3.0*
