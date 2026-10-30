# A tree automaton with infinitely many rogue fixed points

Paper: *A tree automaton with infinitely many rogue fixed points under every supercritical offspring law*, by Yitzchak Shmalo.

Johnson, Podder and Skerman (Random Structures & Algorithms 56 (2020); arXiv:1808.03019, Section 6) asked whether a tree automaton with finitely many, but at least three, states can have infinitely many fixed points of its automaton distributional map. Under binary branching an example is implicit in Rozikov and Zada's work on quadratic stochastic operators (arXiv:0712.3979, Theorem 2(iv)). This paper extends it to every supercritical offspring law.

On states {a, b1, ..., bk}, let a vertex take state bi when exactly one of its children is in a state b and that child is in state bi, and state a otherwise. If the offspring law has generating function G and mean greater than one, the fixed points are the point mass at a and every law that gives a the mass s*, where G'(s*) = 1, a (k-1)-dimensional simplex. For Poisson offspring with mean λ > 1, s* = 1 - (log λ)/λ. Every point of the simplex is rogue, with no moment assumption, so the point mass at a is the only interpretable fixed point. For the points other than the vertices there is a short direct proof; the vertices use the spine argument of Johnson, Podder and Skerman, whose logarithmic moment assumption is not needed for this automaton. A remark shows that two states already suffice when the offspring law has infinite mean, and that three are needed for Poisson offspring.

Files: `paper.pdf`, and its source `paper.tex` with `references.bib` (amsart, cleveref; compile with pdflatex and bibtex).
