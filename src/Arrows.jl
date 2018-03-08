export ⟺

# All possible operands of an arrow operator
const ReactionSide = Union{Environment, SciCompDSL.Operation, Species}

# ⟺ No assumption of mass kinetics, assumes reaction reversible
⟺(lhs::ReactionSide, rhs::ReactionSide) = arrowop(:⟺, lhs, rhs)

# Determines whether lhs is product or substrate (reversible considered forward)
const forward = Dict(:⟺ => true)

# Determines whether arrow is open or closed
const open = Dict(:⟺ => true)

# Determines whether arrow is reversible or not
const reversible = Dict(:⟺ => true)

# Create a Formula from an arrow operator. Inherit compartment from first species
function arrowop(arrow::Symbol, lhs::ReactionSide, rhs::ReactionSide)
    reactlhs = getreactants(lhs)
    reactrhs = getreactants(rhs)
    compartment = length(reactlhs) > 0 ? reactlhs[1].species.compartment : reactrhs[1].species.compartment
    if forward[arrow]
        Formula(arrow, compartment, reactlhs, reactrhs)
    else
        Formula(arrow, compartment, reactrhs, reactlhs)
    end
end

# Parsing the AST of the formula (limited DSL)
getreactants(r::Environment) = Reactant[]
getreactants(r::Species) = [Reactant(r, SciCompDSL.Constant(1))]
function getreactants(r::SciCompDSL.Operation)
    if r.op == +
        return [getreactants(arg)[1] for arg in r.args]
    elseif r.op == *
        return [Reactant(r.args[2], r.args[1])]
    else
        error("Only operators + and * are support in reaction formulae")
    end
end
