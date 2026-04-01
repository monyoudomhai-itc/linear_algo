# Chapter III: Vector Space — Study Guide

**Course**: Linear Algebra | **Institution**: Institute of Technology of Cambodia  
**Instructor**: PHAUK Sokkhey, Ph.D | **Department**: Applied Mathematics and Statistics

---

## Overview

This chapter builds the foundation of **vector space theory** — one of the most powerful frameworks in mathematics. The ideas here appear everywhere: computer graphics, machine learning, physics simulations, signal processing, and more. The chapter is divided into four major sections:

1. Definitions (what a vector space is)
2. Subspace (vector spaces living inside other vector spaces)
3. Spanning Sets and Linear Independence (what vectors can "reach" and whether they overlap)
4. Bases and Dimension (the minimal, efficient description of a space)

---

## Section 1: Definitions — What Is a Vector Space?

### The Core Idea

A **vector space** is any set $V$ where you can add elements together and multiply them by scalars (numbers), as long as those two operations follow 10 specific rules. The set of all real-valued vectors like $(1, 2, 3)$ is the classic example, but polynomials, matrices, and even functions can form vector spaces too.

### The 10 Axioms (Rules)

For $u, v, w \in V$ and scalars $\alpha, \beta \in \mathbb{K}$:

| # | Rule | What It Means |
|---|------|---------------|
| 1 | $u + v \in V$ | Adding two vectors stays inside $V$ (closure under addition) |
| 2 | $u + v = v + u$ | Order of addition doesn't matter (commutativity) |
| 3 | $(u+v)+w = u+(v+w)$ | Grouping doesn't matter (associativity) |
| 4 | $\exists\, 0 \in V: v + 0 = v$ | There is a zero vector (additive identity) |
| 5 | $\exists\, {-v} \in V: v+(-v)=0$ | Every vector has an opposite (additive inverse) |
| 6 | $\alpha \cdot v \in V$ | Scaling a vector stays inside $V$ (closure under scalar multiplication) |
| 7 | $\alpha(u+v)=\alpha u + \alpha v$ | Scalar distributes over vector addition |
| 8 | $(\alpha+\beta)v = \alpha v + \beta v$ | Vector distributes over scalar addition |
| 9 | $(\alpha\beta)v = \alpha(\beta v)$ | Scalar multiplication is associative |
| 10 | $1 \cdot v = v$ | Multiplying by 1 does nothing (scalar identity) |

### Worked Example

**Given**: $\mathbf{u} = (2,-1,5,0)$, $\mathbf{v} = (4,3,1,-1)$, $\mathbf{w} = (-6,2,0,3)$ in $\mathbb{R}^4$.

**Find** $\mathbf{x} = 2\mathbf{u} - (\mathbf{v} + 3\mathbf{w})$:

$$\mathbf{x} = 2(2,-1,5,0) - (4,3,1,-1) - 3(-6,2,0,3)$$
$$= (4,-2,10,0) - (4,3,1,-1) - (-18,6,0,9)$$
$$= (4-4+18,\ -2-3-6,\ 10-1-0,\ 0+1-9)$$
$$= (18,\ -11,\ 9,\ -8)$$

**Find** $\mathbf{x}$ from $3(\mathbf{x} + \mathbf{w}) = 2\mathbf{u} - \mathbf{v} + \mathbf{x}$:

$$3\mathbf{x} + 3\mathbf{w} = 2\mathbf{u} - \mathbf{v} + \mathbf{x}$$
$$2\mathbf{x} = 2\mathbf{u} - \mathbf{v} - 3\mathbf{w}$$
$$\mathbf{x} = \tfrac{1}{2}(18,-11,9,-8) = \left(9,\ -\tfrac{11}{2},\ \tfrac{9}{2},\ -4\right)$$

### Key Properties (Theorem 1)

These follow automatically from the 10 axioms — you do not need to prove them separately:

- The zero vector **0** is unique.
- The additive inverse $-v$ is unique.
- $0 \cdot v = \mathbf{0}$ (the scalar zero times any vector gives the zero vector)
- $c \cdot \mathbf{0} = \mathbf{0}$ (any scalar times the zero vector gives the zero vector)
- If $cv = \mathbf{0}$, then $c = 0$ **or** $v = \mathbf{0}$.
- $-(-v) = v$

---

## Section 2: Subspace

### The Core Idea

A **subspace** is a subset $S$ of a vector space $V$ that is itself a vector space, using the same operations. Instead of checking all 10 axioms, we only need 2 conditions (Theorem 2):

> $S$ is a subspace of $V$ **if and only if**:
> 1. $S \neq \emptyset$ (it is not empty — equivalently, it contains the zero vector)
> 2. For any $u, v \in S$ and scalar $\alpha$: both $u + v \in S$ and $\alpha v \in S$

In plain words: a subspace must be **closed under addition** and **closed under scalar multiplication**.

### Worked Example — Symmetric Matrices (Example 2)

**Claim**: The set $W$ of all $2 \times 2$ symmetric matrices (where $A = A^T$) is a subspace of the vector space $M_{2,2}$ of all $2 \times 2$ matrices.

**Proof**:

- $W$ is nonempty: the zero matrix $\begin{bmatrix}0&0\\0&0\end{bmatrix}$ is symmetric.
- Closed under addition: if $A_1 = A_1^T$ and $A_2 = A_2^T$, then $(A_1+A_2)^T = A_1^T+A_2^T = A_1+A_2$ ✓
- Closed under scalar multiplication: $(cA)^T = cA^T = cA$ ✓

Therefore $W$ is a subspace.

### Worked Example — The First Quadrant (Example 3)

**Claim**: $W = \{(x_1, x_2) : x_1 \geq 0,\ x_2 \geq 0\}$ is **not** a subspace of $\mathbb{R}^2$.

**Why it fails**: $(1,1) \in W$, but $(-1)(1,1) = (-1,-1) \notin W$.

It fails scalar multiplication closure — multiplying by a negative scalar sends you outside the set.

> **Key Insight**: Even though the first quadrant looks "natural," subspaces must be closed under *all* scalars, including negatives. Every subspace must contain the origin and be symmetric around it.

---

## Section 3: Spanning Sets and Linear Independence

### Linear Combination and Span

A vector $v$ is a **linear combination** of $v_1, v_2, \ldots, v_p$ if there exist scalars $k_1, \ldots, k_p$ such that:

$$v = k_1 v_1 + k_2 v_2 + \cdots + k_p v_p$$

The **span** of $\{v_1, \ldots, v_p\}$ is the set of *all* such linear combinations — every vector reachable by scaling and adding those vectors:

$$\text{Span}\{v_1, \ldots, v_p\} = \{k_1 v_1 + \cdots + k_p v_p : k_i \in \mathbb{K}\}$$

### Examples of Spanning Sets (Example 5)

**a.** The set $S = \{(1,0,0),\ (0,1,0),\ (0,0,1)\}$ spans $\mathbb{R}^3$ because any vector $(u_1, u_2, u_3)$ can be written as:
$$\mathbf{u} = u_1(1,0,0) + u_2(0,1,0) + u_3(0,0,1)$$

**b.** The set $S = \{1,\ x,\ x^2\}$ spans $P_2$ (polynomials of degree $\leq 2$) because:
$$p(x) = a(1) + b(x) + c(x^2) = a + bx + cx^2$$

### Redundancy in Spanning Sets (Example 6)

In $\mathbb{R}^3$, with $v_1 = (1,3,1)$, $v_2 = (0,1,2)$, $v_3 = (1,0,-5)$:

$$v_1 = 3v_2 + v_3 = 3(0,1,2)+(1,0,-5) = (1,3,1) \checkmark$$

So $v_1$ is redundant — it is already expressible from the others. This leads naturally to the concept of **linear independence**.

### Linear Independence (Definition 7)

Vectors $v_1, \ldots, v_p$ are **linearly independent** if the only solution to:

$$k_1 v_1 + k_2 v_2 + \cdots + k_p v_p = \mathbf{0}$$

is the **trivial solution** $k_1 = k_2 = \cdots = k_p = 0$.

If any non-trivial solution exists (some $k_i \neq 0$), the vectors are **linearly dependent** — meaning at least one vector is a linear combination of the others and carries no new information.

### Checking Independence — Gauss-Jordan Elimination (Example 8)

**Problem**: Are $\{(1,2,3),\ (0,1,2),\ (-2,0,1)\}$ linearly independent?

Set up $c_1 v_1 + c_2 v_2 + c_3 v_3 = \mathbf{0}$, which gives:

$$\begin{bmatrix}1 & 0 & -2 & 0\\ 2 & 1 & 0 & 0\\ 3 & 2 & 1 & 0\end{bmatrix} \longrightarrow \begin{bmatrix}1 & 0 & 0 & 0\\ 0 & 1 & 0 & 0\\ 0 & 0 & 1 & 0\end{bmatrix}$$

The only solution is $c_1 = c_2 = c_3 = 0$, so the vectors are **linearly independent** ✓.

> **Rule of Thumb**: Reduce the matrix formed by the vectors (as rows or columns). If the RREF has no free variables, the vectors are independent. A free variable means dependence.

---

## Section 4: Bases and Dimension

### What Is a Basis?

A **basis** for a vector space $V$ is a set $\{v_1, v_2, \ldots, v_n\} \subset V$ that is:

1. **Linearly independent** — no redundancy
2. **Spans** $V$ — every vector in $V$ can be reached

A basis is the most efficient description of a space: just enough vectors, with no overlap.

### Standard Bases

| Space | Standard Basis | Size |
|-------|---------------|------|
| $\mathbb{R}^2$ | $\{(1,0),\ (0,1)\}$ | 2 |
| $\mathbb{R}^3$ | $\{(1,0,0),\ (0,1,0),\ (0,0,1)\}$ | 3 |
| $\mathbb{R}^n$ | $\{e_1, e_2, \ldots, e_n\}$ (standard unit vectors) | $n$ |
| $P_2$ | $\{1,\ x,\ x^2\}$ | 3 |
| $P_3$ | $\{1,\ x,\ x^2,\ x^3\}$ | 4 |
| $M_{2,2}$ | Four elementary matrices | 4 |

### Verifying a Basis (Example 10)

**Show** $S = \{(1,0,0),\ (0,1,0),\ (0,0,1)\}$ is a basis for $\mathbb{R}^3$:

- *Spans $\mathbb{R}^3$*: Every vector $(u_1, u_2, u_3) = u_1 e_1 + u_2 e_2 + u_3 e_3$ ✓
- *Linearly independent*: $c_1(1,0,0)+c_2(0,1,0)+c_3(0,0,1)=\mathbf{0}$ gives $c_1=c_2=c_3=0$ ✓

### Unique Representation (Theorem 3)

If $S = \{v_1, \ldots, v_n\}$ is a basis for $V$, then every vector in $V$ can be written in **exactly one way** as a linear combination of vectors in $S$.

This is the coordinate system property — each basis gives a unique "address" to every vector.

### Dimension (Definition 12)

The **dimension** of a vector space $V$, written $\dim V$, is the number of vectors in any basis for $V$.

$$\dim \mathbb{R}^n = n \qquad \dim P_n = n+1 \qquad \dim M_{m,n} = mn$$

This is well-defined because of Theorem 5: all bases for the same vector space have the same number of vectors.

### Theorems on Dimension

**Theorem 4**: If $\dim V = n$, then any set of **more than $n$ vectors** in $V$ is linearly dependent.

**Theorem 5**: If $V$ has one basis with $n$ vectors, then **every** basis for $V$ has exactly $n$ vectors.

### Using Dimension to Rule Out Bases (Example 13)

**a.** Is $S_1 = \{(3,2,1),\ (7,-1,4)\}$ a basis for $\mathbb{R}^3$?

No — $\dim \mathbb{R}^3 = 3$, but $S_1$ has only 2 vectors. A set of 2 vectors cannot span a 3-dimensional space.

**b.** Is $S_2 = \{2+x,\ x^2,\ -1+x^3,\ 1+3x,\ 3-2x+x^2\}$ a basis for $P_3$?

No — $\dim P_3 = 4$, but $S_2$ has 5 vectors. By Theorem 4, any 5 vectors in a 4-dimensional space must be linearly dependent.

---

## Quick Reference: Rules Summary

| Concept | Key Condition | Fail Condition |
|---------|--------------|---------------|
| **Vector Space** | All 10 axioms hold | Any one axiom fails |
| **Subspace** | Non-empty + closed under $+$ and $\cdot$ | Not containing **0**, or escaping under negation |
| **Linear Independence** | Only trivial solution $c_i = 0$ | Some $c_i \neq 0$ solves the equation |
| **Spanning Set** | Every vector in $V$ is reachable | Some vector can't be expressed |
| **Basis** | Independent **and** spans $V$ | Too few (can't span), too many (dependent) |
| **Dimension** | Number of basis vectors | Not unique — guaranteed by Theorem 5 |

---

## Real-World Connections

**Digital Signal Processing**: Sampling from a union of vector subspaces (rather than a single space) improves signal reconstruction in radar and wireless communications.

**RGB Color Model**: Any visible color is a linear combination $c_1 \mathbf{r} + c_2 \mathbf{g} + c_3 \mathbf{b}$ of the standard basis vectors for $\mathbb{R}^3$, where $\mathbf{r} = (1,0,0)$, $\mathbf{g} = (0,1,0)$, $\mathbf{b} = (0,0,1)$. When $c_1 = c_2 = c_3$, the result is grayscale.

---

*Linear Algebra — Institute of Technology of Cambodia*