export ⟺

# All possible operands of an arrow operator
const ReactionSide = Union{Environment, SciCompDSL.Operation, SciCompDSL.Variable}

# ⟺ No assumption of mass kinetics, assumes reaction reversible
⟺(lhs::ReactionSide, rhs::ReactionSide) = arrowop(:⟺, lhs, rhs)

# Determines whether lhs is product or substrate (reversible considered forward)
const forward = Dict(:⟺ => true)

# Determines whether arrow is open or closed
const open = Dict(:⟺ => true)

# Determines whether arrow is reversible or not
const reversible = Dict(:⟺ => true)

# Create a Formula from an arrow operator
function arrowop(arrow::Symbol, lhs::ReactionSide, rhs::ReactionSide)
    if forward[arrow]
        Formula(arrow, getreactants(lhs), getreactants(rhs))
    else
        Formula(arrow, getreactants(rhs), getreactants(lhs))
    end
end

# Parsing the AST of the formula (limited DSL)
getreactants(r::Environment) = Reactant[]
getreactants(r::SciCompDSL.Variable) = [Reactant(r, SciCompDSL.Constant(1))]
function getreactants(r::SciCompDSL.Operation)
    if r.op == +
        return [getreactants(arg)[1] for arg in r.args]
    elseif r.op == *
        return [Reactant(r.args[2], r.args[1])]
    else
        error("Only operators + and * are support in reaction formulae")
    end
end
