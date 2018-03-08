using SciCompDSL
using ReactionDSL
import ReactionDSL: Formula, Reactant, Reaction
using Base.Test

# Basic elements
@Species s1 s2 s3 = 1
@test s1 == Species(:s1)
@test s3 == Species(:s3, 1)

# Basic reaction formula
f1 = 2s1 + s2 ⟺ s3
f2 = Formula(:⟺, [Reactant(s1, Constant(2)), Reactant(s2, Constant(1))],
             [Reactant(s3, Constant(1))])
@test f1 == f2

# Use parameter as stoichiometry
@Param st
f1 = 2s1 + st*s2 ⟺ s3
f2 = Formula(:⟺, [Reactant(s1, Constant(2)), Reactant(s2, st)],
             [Reactant(s3, Constant(1))])
@test f1 == f2

# An actual reaction
@Param k Keq
r1 = 2s1 + st*s2 ⟺ s3 ~ k*(s1*s2 - s3/Keq)
op = k*(s1*s2 - s3/Keq)
r2 = Reaction(f1, op)
@test r1 == r2

# A simple system with sinks and sources
@Param gₓ xₒ kₓ Kₑ vₛ Kₘ
@IVar t
@Species x y
eq = [∅ ⟺ x ~ gₓ*(xₒ - x),
      x ⟺ y ~ kₓ*(x - y/Kₑ),
      y ⟺ ∅ ~ vₛ*y/(Kₘ + y)]
ode1 = ReactionODE(eq, t)
generate_ode_function(ode1)

# Like before, but adding a differential equation and intermediate variable
@Param vₛₘ Kₘₛ kₐ
@DVar vₛ # overrides previous definition
@Var vₛₛₛ
@Deriv D'~t
eq2 = [∅ ⟺ x ~ gₓ*(xₒ - x),
       x ⟺ y ~ kₓ*(x - y/Kₑ),
       vₛₛₛ ~ vₛₘ*xₒ/(vₛ + xₒ),
       D*vₛ ~ (vₛ - vₛₛₛ)*kₐ,
       y ⟺ ∅ ~ vₛ*y/(Kₘ + y)]
ode2 = ReactionODE(eq2, t)
generate_ode_function(ode2)
