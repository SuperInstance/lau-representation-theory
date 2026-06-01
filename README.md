# lau-representation-theory

**Representation theory of finite groups** — a Rust library for computing group representations, character tables, irreducible decompositions, induced representations, Young tableaux, Clebsch-Gordan coefficients, and agent symmetry analysis.

## What This Does

This crate implements the core machinery of finite-group representation theory over the complex numbers:

- **Group algebra** — define finite groups via Cayley tables (with built-in generators for ℤ/nℤ, S₃, and the Klein four-group)
- **Representations** — construct complex matrix representations ρ: G → GL(n, ℂ) and verify homomorphism properties
- **Character theory** — compute characters χ(g) = Tr(ρ(g)), build character tables, and verify both orthogonality relations
- **Irreducible decomposition** — decompose any representation into irreducibles using inner products of characters (Maschke's theorem)
- **Induced representations** — compute induced characters and representations Indᴴᴳ(χ), verify Frobenius reciprocity
- **Tensor products** — Kronecker products of representations, Clebsch-Gordan coefficients for SU(2) (spin-½ × spin-½, spin-1 × spin-1, spin-1 × spin-½), and symmetric powers
- **Young tableaux** — partitions, hook-length formula, standard tableaux, and Murnaghan-Nakayama rule for Sₙ characters
- **Agent symmetry** — decompose multi-agent state spaces under group actions, find invariant subspaces, compute symmetrizers/antisymmetrizers

## Key Idea

Every representation of a finite group over ℂ decomposes uniquely into irreducible building blocks. This library computes that decomposition algorithmically: build the character table, then use the orthogonality relations to read off multiplicities. The number of irreducibles equals the number of conjugacy classes, so the problem is finite and tractable for small groups.

The Young tableau module implements the hook-length formula and Murnaghan-Nakayama rule, connecting symmetric group combinatorics to representation theory without constructing matrices at all.

## Install

```toml
[dependencies]
lau-representation-theory = "0.1.0"
```

Requires Rust 2021 edition. Depends on `nalgebra` (linear algebra), `num-complex` (complex numbers), and `serde` (serialization).

## Quick Start

```rust
use lau_representation_theory::*;

// Create the symmetric group S₃
let s3 = group::FiniteGroup::s3();

// Build the three irreducible representations: trivial, sign, and standard
let trivial = Representation::trivial(s3.order());
let sign = Representation::sign_s3();
let standard = Representation::standard_s3();

// Compute characters
let chi_triv = Character::from_representation(&trivial);
let chi_sign = Character::from_representation(&sign);
let chi_std = Character::from_representation(&standard);

// Verify orthogonality of the character table
let table = CharacterTable::new(&s3, vec![chi_triv, chi_sign, chi_std]);
assert!(table.verify_first_orthogonality());
assert!(table.verify_second_orthogonality());

// Decompose a tensor product into irreducibles
let tensor = tensor::tensor_product(&standard, &standard);
let irr_reps = vec![trivial, sign, standard];
let mults = decomposition::decompose(&tensor, &irr_reps, &s3);
println!("Tensor product decomposes as: {:?}", mults);
// → [1, 1, 1]  (trivial ⊕ sign ⊕ standard)

// Young tableaux: dimension of the S₆ irrep for partition [3,2,1]
let diagram = young::YoungDiagram::new(vec![3, 2, 1]);
assert_eq!(diagram.dimension(), 16);

// Clebsch-Gordan: spin-½ ⊗ spin-½ = spin-1 ⊕ spin-0
let cg = tensor::clebsch_gordan(0.5, 0.5);
assert_eq!(cg, vec![0.0, 1.0]);
```

## API Reference

### `group` — Finite Group Definitions

| Type / Fn | Description |
|-----------|-------------|
| `FiniteGroup` | A finite group stored as a Cayley table with labeled elements |
| `FiniteGroup::cyclic(n)` | Construct ℤ/nℤ |
| `FiniteGroup::s3()` | Construct the symmetric group S₃ |
| `FiniteGroup::klein_four()` | Construct ℤ/2ℤ × ℤ/2ℤ |
| `.order()`, `.identity()`, `.inverse(i)` | Basic group operations |
| `.conjugacy_classes()` | Partition elements into conjugacy classes |
| `GroupOps` | Trait abstracting group operations for generic code |

### `representation` — Group Representations

| Type / Fn | Description |
|-----------|-------------|
| `Representation` | A complex matrix representation ρ(g) ∈ GL(n, ℂ) for each group element |
| `Representation::trivial(n)` | 1D trivial representation (all elements → [1]) |
| `Representation::sign_s3()` | 1D sign representation of S₃ |
| `Representation::standard_s3()` | 2D standard representation of S₃ |
| `.direct_sum(&other)` | Direct sum of two representations |
| `.verify_homomorphism(group)` | Check ρ(gh) = ρ(g)ρ(h) and ρ(e) = I |
| `.trace(i)` | Character value χ(g_i) |

### `character` — Character Theory

| Type / Fn | Description |
|-----------|-------------|
| `Character` | Character χ(g) = Tr(ρ(g)), stored as values per element |
| `Character::from_representation(rep)` | Extract character from a representation |
| `.inner_product(other, group)` | Standard inner product ⟨χ₁, χ₂⟩ = (1/|G|) Σ χ₁(g) χ̄₂(g) |
| `.is_irreducible(group)` | Check ⟨χ, χ⟩ = 1 |
| `CharacterTable` | Full character table indexed by conjugacy classes |
| `.verify_first_orthogonality()` | Verify Σ_g χᵢ(g)χ̄ⱼ(g) = |G| δᵢⱼ |
| `.verify_second_orthogonality()` | Verify Σ_χ χ(gᵢ)χ̄(gⱼ) = |G|/|C(gᵢ)| δᵢⱼ |
| `.decompose_character(chi)` | Decompose a character into irreducible multiplicities |

### `irreducible` — Schur's Lemma and Isotypic Projections

| Fn | Description |
|----|-------------|
| `is_irreducible(rep, group)` | Test irreducibility via ⟨χ, χ⟩ = 1 |
| `schur_lemma(rep1, rep2, T, group)` | Verify an intertwining map; returns scalar λ if T = λI |
| `isotypic_projection(rep, chi, group)` | Compute Pᵢ = (dim χᵢ / |G|) Σ χ̄ᵢ(g) ρ(g) |

### `decomposition` — Maschke's Theorem

| Fn | Description |
|----|-------------|
| `decompose(rep, irr_reps, group)` | Returns multiplicities of each irreducible in rep |
| `verify_maschke(rep, irr_reps, group)` | Check Σ(mᵢ × dim) = dim(rep) |
| `block_diagonalize(rep, irr_reps, group)` | Compute block structure of the decomposition |

### `induced` — Induced Representations and Frobenius Reciprocity

| Type / Fn | Description |
|-----------|-------------|
| `SubgroupEmbedding` | Maps subgroup element indices → group element indices |
| `induced_character_embedded(group, χₕ, emb)` | Compute Indᴴᴳ(χ) |
| `restrict_character_embedded(χ_G, emb)` | Compute Resᴴᴳ(χ) |
| `induced_representation_embedded(...)` | Construct the induced representation matrix-by-matrix |
| `verify_frobenius_reciprocity_embedded(...)` | Check ⟨Ind χ, ψ⟩_G = ⟨χ, Res ψ⟩_H |

Also provides legacy label-based APIs: `induced_character`, `restrict_character`, `verify_frobenius_reciprocity`.

### `tensor` — Tensor Products and Clebsch-Gordan

| Fn | Description |
|----|-------------|
| `tensor_product(rep1, rep2)` | Kronecker product ρ₁ ⊗ ρ₂ |
| `tensor_character(chi1, chi2)` | Pointwise product χ₁ · χ₂ |
| `clebsch_gordan(j1, j2)` | SU(2) CG decomposition: returns list of j values |
| `clebsch_gordan_coefficient(j1, m1, j2, m2, j, m)` | Individual CG coefficient C(j1 m1 j2 m2 \| J M) |
| `verify_cg_dimension(j1, j2)` | Check dim conservation: (2j₁+1)(2j₂+1) = Σ(2j+1) |
| `symmetric_power(rep, n)` | Symmetric n-th power of a representation |

### `young` — Young Tableaux and Symmetric Group Characters

| Type / Fn | Description |
|-----------|-------------|
| `YoungDiagram` | A partition (non-increasing row lengths) |
| `.conjugate()` | Transpose (flip rows ↔ columns) |
| `.hook_length(i, j)`, `.hook_product()` | Hook lengths and their product |
| `.dimension()` | dim = n! / hook_product (hook-length formula) |
| `YoungDiagram::partitions(n)` | Enumerate all partitions of n |
| `YoungTableau` | A filling of a Young diagram |
| `.is_standard()` | Check strict increase in rows and columns |
| `YoungTableau::standard_tableaux(diagram)` | Enumerate all standard tableaux of a shape |
| `symmetric_group_character(partition, cycle_type)` | Sₙ character via Murnaghan-Nakayama rule |

### `agent_symmetry` — Multi-Agent Symmetry Analysis

| Type / Fn | Description |
|-----------|-------------|
| `AgentStateSpace` | Multi-agent state space with group action |
| `.permutation_representation(group)` | Natural permutation action on agent indices |
| `analyze_symmetry(space, group, irr_reps)` | Full symmetry decomposition of the state space |
| `invariant_subspaces(rep, irr_reps, group)` | List invariant subspaces with dimensions |
| `symmetrizer(rep, group)` | Projection onto the fully symmetric subspace |
| `antisymmetrizer(rep, signs, group)` | Projection onto the alternating subspace |

## How It Works

1. **Groups are Cayley tables.** A `FiniteGroup` stores the full multiplication table. This is simple and correct but O(n²) in group order — fine for the small groups (≤ dozens of elements) where you'd do explicit representation theory.

2. **Representations are vectors of matrices.** Each group element maps to a `DMatrix<Complex64>` via `nalgebra`. Characters are extracted by taking traces.

3. **Decomposition uses inner products.** The multiplicity of irreducible χᵢ in χ is nᵢ = ⟨χ, χᵢ⟩ = (1/|G|) Σ χ(g)χ̄ᵢ(g). This is always a non-negative integer (proven by the orthogonality theorem).

4. **Young tableaux bypass matrices.** The Murnaghan-Nakayama rule computes Sₙ characters directly from partitions and cycle types using border-strip combinatorics — no matrix construction needed.

5. **Clebsch-Gordan uses known tables.** For the most common cases (spin-½ × spin-½, spin-1 × spin-1, spin-1 × spin-½), the library uses explicit tabulated coefficients. The general case falls back to the Wigner 3j-symbol formula.

## The Math

**Maschke's Theorem:** Every finite-dimensional representation of a finite group over ℂ is completely reducible: ρ ≅ n₁ρ₁ ⊕ n₂ρ₂ ⊕ ⋯, where the ρᵢ are pairwise non-isomorphic irreducibles.

**Character Orthogonality (First Relation):** Σ_{g∈G} χᵢ(g)χ̄ⱼ(g) = |G| δᵢⱼ, where the sum runs over all group elements.

**Frobenius Reciprocity:** For a subgroup H ≤ G and characters χ of H, ψ of G: ⟨Indᴴᴺ(χ), ψ⟩_G = ⟨χ, Resᴴᴺ(ψ)⟩_H.

**Hook-Length Formula:** The dimension of the irreducible Sₙ representation indexed by partition λ is n! / ∏ h(i,j), where the product is over all boxes of the Young diagram and h(i,j) is the hook length.

**Murnaghan-Nakayama Rule:** χ^λ(μ) = Σ (-1)ˡ h(μ), summed over all border strips of size μ₁ that can be removed from λ, where l is the leg length of the strip.

**Clebsch-Gordan Series:** For SU(2), V_{j₁} ⊗ V_{j₂} = ⊕_{j=|j₁-j₂|}^{j₁+j₂} V_j. The CG coefficients C(j₁m₁ j₂m₂ | J M) are the expansion coefficients coupling the product basis to the total angular momentum basis.

## Test Suite

67 tests covering:
- Group structure (identity, inverse, conjugacy classes for cyclic, S₃, Klein four)
- Representation verification (homomorphism checks for trivial, sign, standard)
- Character orthogonality (both relations)
- Irreducibility testing and Schur's lemma
- Decomposition and Maschke verification
- Induced characters and Frobenius reciprocity
- Tensor products and CG coefficient values
- Young diagram/tableaux construction and hook-length formula
- Murnaghan-Nakayama rule for S₃ characters
- Agent symmetry decomposition

## License

MIT
